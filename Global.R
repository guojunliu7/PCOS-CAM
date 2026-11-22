# PCOS-CAM cluster scoring, adapted from guojunliu7/FPCAM (Apache-2.0).
read_pcos_dictionary <- function(path) {
  d <- as.data.frame(readxl::read_excel(path, sheet = "Dictionary"))
  m <- as.data.frame(readxl::read_excel(path, sheet = "Mouse_to_human"))
  required <- c("cell_cluster", "marker_gene", "Number_of_votes_received", "votes", "percent")
  if (!all(required %in% names(d)) || !all(c("mouse_gene", "human_gene") %in% names(m)))
    stop("字典缺少标记、权重或人鼠映射列。")
  if (!nrow(d) || anyNA(d[required]) || any(!nzchar(d$cell_cluster)) ||
      any(!nzchar(d$marker_gene)) || anyDuplicated(d[c("cell_cluster", "marker_gene")]) ||
      any(!is.finite(d$percent) | d$percent <= 0) ||
      any(d$votes <= 0 | d$Number_of_votes_received <= 0) ||
      any(abs(d$percent - d$Number_of_votes_received / d$votes) > 1e-10) ||
      any(abs(tapply(d$percent, d$cell_cluster, sum) - 1) > 1e-10))
    stop("字典权重或标签–基因组合异常。")
  if (anyNA(m[c("mouse_gene", "human_gene")]) || anyDuplicated(m$mouse_gene) ||
      anyDuplicated(m$human_gene)) stop("小鼠映射须为唯一的一对一配对。")
  list(markers = d, ortholog = setNames(m$human_gene, m$mouse_gene))
}

read_pcos_expression <- function(path, filename = path) {
  ext <- tolower(tools::file_ext(filename))
  x <- switch(ext,
    xlsx = as.data.frame(readxl::read_excel(path, .name_repair = "minimal")),
    csv = read.csv(path, check.names = FALSE, colClasses = "character", na.strings = c("", "NA"), stringsAsFactors = FALSE),
    tsv = read.delim(path, check.names = FALSE, colClasses = "character", na.strings = c("", "NA"), stringsAsFactors = FALSE),
    stop("支持 XLSX、CSV 和 TSV。"))
  if (!nrow(x) || anyDuplicated(names(x)) || any(!nzchar(names(x))))
    stop("输入为空或存在空白、重复列名。")
  if (all(c("cluster", "gene", "expression") %in% names(x)))
    return(x[c("cluster", "gene", "expression")])
  if (any(grepl("_Markergene$", names(x), ignore.case = TRUE))) {
    genes <- seq.int(1L, ncol(x), by = 2L)
    ids <- sub("_Markergene$", "", names(x)[genes], ignore.case = TRUE)
    if (ncol(x) %% 2 || !all(grepl("_Markergene$", names(x)[genes], ignore.case = TRUE)) ||
        !all(tolower(names(x)[genes+1L]) == tolower(paste0(ids, "_Express"))))
      stop("每个簇需要相邻、同名的 _Markergene 和 _Express 列。")
    parts <- lapply(seq_along(ids), function(i) {
      g <- x[[genes[i]]]; v <- x[[genes[i]+1L]]
      keep <- !(is.na(g) & is.na(v))
      if (!any(keep)) stop("空簇：", ids[i])
      data.frame(cluster = ids[i], gene = g, expression = v)[keep,,drop = FALSE]
    })
    return(do.call(rbind, parts))
  }
  if (ncol(x) < 2 || tolower(names(x)[1]) != "gene")
    stop("宽表首列应为 gene，后续列为各簇的平均表达。")
  do.call(rbind, lapply(names(x)[-1], function(id)
    data.frame(cluster = id, gene = x[[1]], expression = x[[id]])))
}

pcos_annotate <- function(expression, dictionary, species = c("Human", "Mouse")) {
  species <- match.arg(species)
  x <- expression
  if (!all(c("cluster","gene","expression") %in% names(x)) || !nrow(x))
    stop("需提供 cluster、gene 和 expression 三列及至少一行数据。")
  x$gene <- trimws(as.character(x$gene))
  x$cluster <- as.character(x$cluster)
  x$expression <- suppressWarnings(as.numeric(as.character(x$expression)))
  if (anyNA(x[c("cluster","gene")]) || any(!nzchar(x$gene)) || any(!nzchar(trimws(x$cluster))) ||
      any(!is.finite(x$expression) | x$expression < 0))
    stop("簇名和基因名须非空，表达量须为有限非负值。")
  if (anyDuplicated(x[c("cluster","gene")])) stop("同簇同基因重复，请先核实并汇总表达量。")
  x$symbol <- if (species == "Mouse") unname(dictionary$ortholog[x$gene]) else x$gene
  d <- dictionary$markers
  ids <- unique(x$cluster); labels <- unique(d$cell_cluster)
  score <- matrix(0, length(labels), length(ids), dimnames = list(labels, ids))
  evidence <- vector("list", length(ids)); calls <- evidence
  for (j in seq_along(ids)) {
    z <- x[x$cluster == ids[j],,drop = FALSE]
    total <- sum(z$expression)
    if (!is.finite(total)) stop("表达总量超出数值范围。")
    idx <- match(d$marker_gene, z$symbol)
    keep <- !is.na(idx) & z$expression[idx] > 0
    hit <- d[keep, c("cell_cluster","marker_gene","percent"),drop = FALSE]
    hit$cluster <- rep(ids[j], nrow(hit))
    hit$input_gene <- z$gene[idx[keep]]
    hit$contribution <- z$expression[idx[keep]] / total * hit$percent
    evidence[[j]] <- hit
    if (nrow(hit)) {
      sums <- tapply(hit$contribution, hit$cell_cluster, sum)
      score[names(sums),j] <- sums
    }
    v <- score[,j]; top <- max(v)
    # Relative tolerance only: very small positive scores must not tie with zero.
    tied <- abs(v-top) <= 32*.Machine$double.eps*top
    winners <- if (top > 0) labels[tied] else character()
    rest <- v[!tied]
    second <- if (sum(tied)>1) top else if (length(rest)) max(rest) else NA_real_
    runners <- if (is.finite(second) && second > 0) labels[abs(v-second) <= 32*.Machine$double.eps*second] else character()
    lower <- if (length(rest)) max(rest) else NA_real_
    delta2 <- if (length(rest) && top-lower >= 1e-5) lower-median(v[v != lower]) else NA_real_
    calls[[j]] <- data.frame(cluster = ids[j],
      annotation = if (!length(winners)) "Unassigned" else if (length(winners)>1) "Ambiguous" else winners,
      top_candidates = paste(winners, collapse = " | "), score = top,
      second_candidates = paste(runners, collapse = " | "), second_score = second, margin = top-second,
      delta1 = if (length(rest)) top-median(rest) else NA_real_, delta2 = delta2,
      input_genes = nrow(z), unmapped_genes = sum(is.na(z$symbol)),
      matched_genes = length(unique(hit$marker_gene)), expression_total = total)
  }
  list(annotations = do.call(rbind,calls), scores = score,
       evidence = do.call(rbind,evidence), unmapped = x[is.na(x$symbol),c("cluster","gene","expression"),drop=FALSE])
}

pcos_cluster_labels <- function(cluster, result) {
  idx <- match(as.character(cluster), result$annotations$cluster)
  if (anyNA(idx)) stop("坐标或对象中的簇编号与评分输入不一致。")
  result$annotations$annotation[idx]
}
