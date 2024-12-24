# UI Fonksiyonu
DifferentialExpressionUI <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(ns("countsFile"), "Upload Gene Expression Data (CSV)", accept = c(".csv")),
    actionButton(ns("runAnalysis"), "Run Differential Expression Analysis"),
    tableOutput(ns("differentialResults")),
    plotOutput(ns("volcanoPlot"))
  )
}
# Server Fonksiyonu
DifferentialExpressionServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # Diferansiyel gen ekspresyonu fonksiyonu
    findDifferentialExpression <- function(counts_file) {
      counts <- read.csv(counts_file, row.names = 1)
      
      # S??tun say??s??n?? al
      num_samples <- ncol(counts)
      
      # Grup etiketlerini dinamik olarak olu??tur
      # Bu ??rnekte her iki grup i??in e??it say??da ??rnek al??yoruz (bu say??y?? de??i??tirebilirsiniz)
      group <- factor(rep(c("Group1", "Group2"), length.out = num_samples))
      
      # DGEList nesnesi olu??tur
      dge <- DGEList(counts = counts, group = group)
      dge <- calcNormFactors(dge)
      
      # Tahmini dispersiyon hesapla
      dge <- estimateDisp(dge)
      
      # GLM fit modelini olu??tur
      fit <- glmQLFit(dge)
      
      # Diferansiyel gen ekspresyonu testi
      results <- glmQLFTest(fit)
      
      # En y??ksek farkl??la??an genler
      top_genes <- topTags(results, n = 10)
      
      # ????kt??y?? daha eri??ilebilir hale getirme
      top_genes_df <- as.data.frame(top_genes)
      
      return(top_genes_df)
    }
    
    # Kullan??c?? "Run Analysis" butonuna t??klad??????nda ??al????acak
    observeEvent(input$runAnalysis, {
      req(input$countsFile)  # Dosyan??n y??klendi??inden emin ol
      
      # Diferansiyel gen ekspresyonu analizini ba??lat
      differential_results <- findDifferentialExpression(input$countsFile$datapath)
      
      # Sonu??lar?? ekranda g??ster
      output$differentialResults <- renderTable({
        differential_results
      })
      
      # Volcano grafi??i olu??turma
      output$volcanoPlot <- renderPlot({
        plot(differential_results$logFC, -log10(differential_results$PValue), 
             xlab = "Log Fold Change", ylab = "-Log10(p-value)", 
             main = "Volcano Plot", pch = 20, col = "blue")
        abline(h = -log10(0.05), col = "red", lty = 2)  # 0.05 p-de??eri i??in k??rm??z?? ??izgi
      })
    })
  })
}
