# PCOS-CAM hierarchical annotation (B). Version A remains in Global.R.
read_pcos_hierarchy <- function(path, dictionary, mapping_path) {
  map <- read.csv(mapping_path, stringsAsFactors = FALSE)
  labels <- unique(dictionary$markers$cell_cluster)
  if (!all(c("cell_cluster","parent") %in% names(map)) || anyNA(map) ||
      any(!nzchar(map$parent)) || anyDuplicated(map$cell_cluster) ||
      !setequal(map$cell_cluster, labels)) stop("父类映射与字典标签不一致。")
  parent <- setNames(map$parent,map$cell_cluster)
  e <- as.data.frame(readxl::read_excel(path, sheet="Evidence",
                                       col_types="text"))
  e <- e[!is.na(e$integrated_gene),]
  e$parent <- unname(parent[e$cell_cluster])
  e$rank <- as.numeric(e$source_rank)
  if (anyNA(e[c("parent","vote_unit","rank")]) || any(e$rank < 1))
    stop("来源证据缺少标签、投票单位或排名。")
  # One best rank per independent voting unit, parent and human-coordinate gene.
  best <- aggregate(rank ~ parent + integrated_gene + vote_unit, e, min)
  keys <- c("parent","integrated_gene")
  n <- aggregate(best["rank"],best[keys],length)
  names(n)[3] <- "votes"
  ranks <- aggregate(best["rank"],best[keys],mean)
  genes <- merge(n,ranks,by=keys)
  genes <- genes[order(genes$parent,-genes$votes,genes$rank,genes$integrated_gene,
                       method="radix"),]
  genes <- do.call(rbind,lapply(split(genes,genes$parent),head,10L))
  markers <- data.frame(cell_cluster=genes$parent,marker_gene=genes$integrated_gene,
    Number_of_votes_received=genes$votes,votes=ave(genes$votes,genes$parent,FUN=sum))
  markers$percent <- markers$Number_of_votes_received/markers$votes
  list(parent=parent,markers=markers)
}

# Summarize one competition using only candidates at the current level.
pcos_evidence_summary <- function(scores,hits) {
  top <- max(scores)
  winners <- function(s) {
    m <- max(s)
    if (m <= 0) return(character())
    names(s)[abs(s-m) <= 32*.Machine$double.eps*m]
  }
  best <- winners(scores)
  runner <- if (length(scores)>1L) sort(scores,decreasing=TRUE)[2] else 0
  margin <- if (length(best)>1L) 0 else unname(top-runner)
  out <- data.frame(winner=if (length(best)==1L) best else NA_character_,
    candidates=paste(best,collapse=" | "),score=top,
    runner_up=paste(names(scores)[runner>0 & abs(scores-runner)<=32*.Machine$double.eps*runner],collapse=" | "),
    runner_score=unname(runner),margin=margin,
    relative_margin=if (top>0) margin/top else 0,
    matched=0L,weight_coverage=0,leading_gene="",leading_share=NA_real_,
    stable=NA,after_removal="")
  if (length(best)!=1L) return(out)
  relevant <- hits[hits$cell_cluster==best,,drop=FALSE]
  out$matched <- nrow(relevant)
  out$weight_coverage <- sum(relevant$percent)
  maximum <- max(relevant$contribution)
  leading <- relevant$marker_gene[abs(relevant$contribution-maximum)<=32*.Machine$double.eps*maximum]
  out$leading_gene <- paste(leading,collapse=" | ")
  out$leading_share <- maximum/top
  # Remove each co-leading gene separately across every competitor. The common
  # expression denominator does not affect the order of the remaining scores.
  alternatives <- lapply(leading,function(g) {
    h <- hits[hits$marker_gene!=g & hits$cell_cluster %in% names(scores),,drop=FALSE]
    updated <- setNames(numeric(length(scores)),names(scores))
    if (nrow(h)) {
      summed <- tapply(h$contribution,h$cell_cluster,sum)
      updated[names(summed)] <- summed
    }
    winners(updated)
  })
  out$stable <- all(vapply(alternatives,function(a) identical(a,best),logical(1)))
  out$after_removal <- paste(unique(vapply(alternatives,function(a)
    if (length(a)) paste(a,collapse=" | ") else "Unassigned",character(1))),collapse=" / ")
  out
}

pcos_granulosa_candidate <- function(annotation, baseline, parent) {
  eligible <- annotation == "Granulosa" & baseline != "Granulosa" &
    baseline %in% names(parent)[parent == "Granulosa"]
  ifelse(eligible, baseline, NA_character_)
}

pcos_annotate_b <- function(expression,dictionary,hierarchy,species=c("Human","Mouse"),
                            min_markers=2L,min_relative_margin=0.10) {
  species <- match.arg(species)
  if (length(min_markers)!=1L || !is.finite(min_markers) || min_markers<1 ||
      min_markers!=floor(min_markers) || length(min_relative_margin)!=1L ||
      !is.finite(min_relative_margin) || min_relative_margin<0 || min_relative_margin>1)
    stop("请检查命中数和相对分差阈值。")
  fine <- pcos_annotate(expression,dictionary,species)
  coarse <- pcos_annotate(expression,list(markers=hierarchy$markers,
                                          ortholog=dictionary$ortholog),species)
  calls <- list(); details <- list()
  assess <- function(z) {
    if (z$score<=0) return("no_hit")
    if (is.na(z$winner)) return("tie")
    if (z$matched<min_markers) return("few_markers")
    if (z$relative_margin<min_relative_margin) return("small_margin")
    "pass"
  }
  for (id in colnames(fine$scores)) {
    c_hits <- subset(coarse$evidence,cluster==id)
    p <- pcos_evidence_summary(coarse$scores[,id],c_hits)
    state <- assess(p)
    annotation <- if (state %in% c("tie","small_margin")) "Ambiguous" else "Unassigned"
    level <- "unresolved"
    reason <- paste0("parent_",state)
    details[[length(details)+1L]] <- data.frame(cluster=id,stage="parent",p)
    if (state=="pass") {
      annotation <- p$winner;level <- "parent"
      members <- names(hierarchy$parent)[hierarchy$parent==p$winner]
      # Include the general label when available: a broad signature can win.
      s <- fine$scores[members,id];names(s) <- members
      f <- pcos_evidence_summary(s,subset(fine$evidence,cluster==id))
      details[[length(details)+1L]] <- data.frame(cluster=id,stage="within_parent",f)
      decision <- assess(f)
      reason <- paste0("within_parent_",decision)
      if (decision=="pass" && f$winner!=p$winner) {
        annotation <- f$winner;level <- "fine"
      } else if (decision=="pass") reason <- "broad_label_wins"
    }
    calls[[id]] <- data.frame(cluster=id,annotation=annotation,level=level,reason=reason,
      baseline_A=fine$annotations$annotation[match(id,fine$annotations$cluster)],
      parent_candidates=p$candidates,parent_score=p$score,parent_margin=p$margin,
      parent_stable=p$stable)
  }
  annotations <- do.call(rbind,calls)
  annotations$candidate_subtype <- pcos_granulosa_candidate(annotations$annotation,
    annotations$baseline_A, hierarchy$parent)
  list(annotations=annotations,diagnostics=do.call(rbind,details),
       scores=fine$scores,parent_scores=coarse$scores,evidence=fine$evidence,
       parent_evidence=coarse$evidence,baseline_A=fine$annotations,unmapped=fine$unmapped,
       parent_dictionary=hierarchy$markers,
       settings=data.frame(species=species,min_markers=min_markers,min_relative_margin=min_relative_margin))
}

# Public workflow. Keep the original scoring outputs available for reproducibility.
pcos_cam <- function(expression, dictionary, hierarchy, species=c("Human","Mouse"),
                     min_markers=2L, min_relative_margin=0.10) {
  result <- pcos_annotate_b(expression, dictionary, hierarchy, match.arg(species),
                            min_markers, min_relative_margin)
  names(result$annotations)[names(result$annotations)=="baseline_A"] <- "global_top_label"
  names(result)[names(result)=="baseline_A"] <- "global_ranking"
  result
}
