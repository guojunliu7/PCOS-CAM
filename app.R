library(shiny)
source("Global.R", local = TRUE)
source("hierarchy.R", local = TRUE)
pcos_dict <- read_pcos_dictionary("data/File_1.xlsx")
hierarchy <- read_pcos_hierarchy("data/File_1.xlsx",pcos_dict,"data/label_hierarchy.csv")
ui <- fluidPage(titlePanel("PCOS-CAM 卵巢细胞簇注释"),
  sidebarLayout(sidebarPanel(selectInput("species","输入物种",c("Human","Mouse")),
    numericInput("min_hits","每级最少命中基因",2,min=1,step=1),
    numericInput("min_gap","每级最小相对分差",0.10,min=0,max=1,step=.01),
    fileInput("expr","标记表达表（XLSX、CSV、TSV）",accept=c(".xlsx",".csv",".tsv")),
    helpText("宽表首列为 gene，其余列为各簇的非负平均表达；也支持同簇基因与表达量成对列或 cluster、gene、expression 长表。"),
    tags$a("下载 PCOS 研究人源输入示例", href="PCOS_GSE240688_cluster_expression.csv", download="PCOS_GSE240688_cluster_expression.csv"),
    fileInput("coords","可选 UMAP 坐标（CSV）",accept=".csv"),
    helpText("每个细胞一行，列名为 cluster、UMAP_1、UMAP_2，簇编号与表达表完全一致。"),
    downloadButton("download","下载评分与证据")),
    mainPanel(tabsetPanel(tabPanel("簇注释",tableOutput("calls")),
      tabPanel("分层证据",tableOutput("diagnostics")),
      tabPanel("得分热图",plotOutput("heat")),tabPanel("UMAP",plotOutput("umap")),
      tabPanel("分差",plotOutput("margin"))))))
server <- function(input, output, session) {
  result <- reactive({
    req(input$expr)
    tryCatch({
      x <- read_pcos_expression(input$expr$datapath,input$expr$name)
      pcos_cam(x,pcos_dict,hierarchy,input$species,input$min_hits,input$min_gap)
    },error=function(e) validate(need(FALSE,conditionMessage(e))))
  })
  output$calls <- renderTable(result()$annotations, digits=6)
  output$diagnostics <- renderTable(result()$diagnostics,digits=4)
  output$heat <- renderPlot({
    s <- result()$parent_scores
    top <- unique(unlist(lapply(seq_len(ncol(s)),function(i) {
      pos <- which(s[,i]>0); head(pos[order(s[pos,i],decreasing=TRUE)],3)
    })))
    validate(need(length(top)>0,"没有匹配到表达大于零的字典标记。"))
    z <- as.data.frame(as.table(s[top,,drop=FALSE])); names(z) <- c("label","cluster","score")
    ggplot2::ggplot(z,ggplot2::aes(cluster,label,fill=score))+ggplot2::geom_tile()+
      ggplot2::scale_fill_gradient(low="white",high="#00665d")+ggplot2::labs(x="Cluster",y=NULL)
  })
  output$margin <- renderPlot({
    z <- result()$annotations
    z$display_margin <- z$parent_margin
    ggplot2::ggplot(z,ggplot2::aes(cluster,display_margin))+ggplot2::geom_col(fill="#00665d")+
      ggplot2::labs(x="Cluster",y="Top score minus runner-up")
  })
  output$umap <- renderPlot({
    req(input$coords); z <- read.csv(input$coords$datapath,colClasses="character",check.names=FALSE)
    validate(need(all(c("cluster","UMAP_1","UMAP_2") %in% names(z)),"缺少坐标列。"))
    for (col in c("UMAP_1","UMAP_2")) z[[col]] <- suppressWarnings(as.numeric(z[[col]]))
    validate(need(nrow(z)>0 && all(is.finite(z$UMAP_1) & is.finite(z$UMAP_2)),"坐标应为有限数值。"))
    z$annotation <- pcos_cluster_labels(z$cluster,result())
    ggplot2::ggplot(z,ggplot2::aes(UMAP_1,UMAP_2,color=annotation))+ggplot2::geom_point(size=.4)+ggplot2::coord_equal()
  })
  output$download <- downloadHandler(filename=function() "PCOS-CAM_results.xlsx",content=function(file) {
    r <- result()
    r$scores <- data.frame(cell_cluster=rownames(r$scores),r$scores,check.names=FALSE)
    if (!is.null(r$parent_scores)) r$parent_scores <- data.frame(parent=rownames(r$parent_scores),r$parent_scores,check.names=FALSE)
    openxlsx::write.xlsx(r,file,overwrite=TRUE)
  })
}
shinyApp(ui,server)
