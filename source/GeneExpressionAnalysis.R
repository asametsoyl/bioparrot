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
      
      # E??er gruplar belirtilmemi??se, her iki grup i??in e??it say??da ??rnek varsayal??m
      group_labels <- factor(rep(c("Group1", "Group2"), length.out = num_samples))
      
      # DGEList nesnesi olu??turma
      dge <- DGEList(counts = counts, group = group_labels)
      dge <- calcNormFactors(dge)
      
      # CPM de??erlerini hesapla
      cpm_values <- cpm(dge, normalized.lib.sizes = TRUE)
      
      # ??statistiksel test ve sonu??lar?? hesapla
      design <- model.matrix(~ group_labels)
      fit <- glmFit(dge, design)
      results <- glmLRT(fit)
      
      # En fazla farkl??la??an 10 genin ????kt??s??
      top_genes <- topTags(results, n = 10)
      
      # Geriye CPM ve Volcano plotu d??nd??r??yoruz
      cpm_plot <- ggplot(data = as.data.frame(t(cpm_values)), aes(x = 1:nrow(cpm_values), y = log10(V1 + 1))) +
        geom_boxplot() +
        labs(title = "Gene Expression CPM Values", x = "Genes", y = "Log10(CPM + 1)") +
        theme_minimal()
      
      # Volcano plotu olu??tur
      volcano_plot <- ggplot(data = top_genes, aes(x = logFC, y = -log10(PValue))) +
        geom_point(aes(color = PValue < 0.05), size = 3) +
        labs(title = "Volcano Plot", x = "Log Fold Change", y = "-Log10(p-value)") +
        theme_minimal() +
        scale_color_manual(values = c("black", "red"))
      
      # Sonu??lar?? d??nd??r??yoruz
      return(list(cpm_values = cpm_values, top_genes = top_genes, cpm_plot = cpm_plot, volcano_plot = volcano_plot))
    }
    
    # Kullan??c?? "Run Analysis" butonuna t??klad??????nda ??al????acak
    observeEvent(input$runAnalysis, {
      req(input$countsFile)  # Dosyan??n y??klendi??inden emin ol
      
      # Diferansiyel gen ekspresyonu analizini ba??lat
      differential_results <- findDifferentialExpression(input$countsFile$datapath)
      
      # Sonu??lar?? ekranda g??ster
      output$differentialResults <- renderTable({
        differential_results$top_genes
      })
      
      # Volcano grafi??i olu??turma
      output$volcanoPlot <- renderPlot({
        differential_results$volcano_plot
      })
    })
  })
}
