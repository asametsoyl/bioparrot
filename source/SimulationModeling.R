library(shiny)
library(ggplot2)

# Sim??lasyon Fonksiyonu
simulateGeneExpression <- function(n_genes = 100, n_samples = 10, 
                                   noise_level = 0.1, 
                                   n_groups = 2, 
                                   group_sizes = c(5, 5), 
                                   mean_diff = 2) {
  # Grup isimleri
  group_names <- paste0("Group_", 1:n_groups)
  
  # Sim??le edilen gen ekspresyon verisi olu??turuluyor
  expression_data <- matrix(NA, nrow = n_genes, ncol = n_samples)
  colnames(expression_data) <- paste0("Sample_", 1:n_samples)
  rownames(expression_data) <- paste0("Gene_", 1:n_genes)
  
  # Gruplar i??in gen ekspresyon de??erlerini olu??turma
  for (i in 1:n_genes) {
    group1_expr <- rnorm(group_sizes[1], mean = 10, sd = noise_level)  # Control grubu
    group2_expr <- rnorm(group_sizes[2], mean = 10 + mean_diff, sd = noise_level)  # Treatment grubu
    
    expression_data[i, 1:group_sizes[1]] <- group1_expr
    expression_data[i, (group_sizes[1]+1):n_samples] <- group2_expr
  }
  
  # Veri g??rselle??tirmeyi iyile??tirme: Gruplar aras??ndaki farklar?? g??stermek
  data_for_plot <- data.frame(
    Gene = rep(rownames(expression_data), each = n_samples),
    Expression = as.vector(expression_data),
    Sample = rep(colnames(expression_data), times = n_genes),
    Group = rep(c(rep("Control", group_sizes[1]), rep("Treatment", group_sizes[2])), each = n_genes)
  )
  
  # Grafik olu??turma
  expression_plot <- ggplot(data_for_plot, aes(x = Gene, y = Expression, color = Group)) +
    geom_point() +
    theme_minimal() +
    labs(title = "Simulated Gene Expression", x = "Genes", y = "Expression Level") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_manual(values = c("blue", "red"))  # Control ve Treatment i??in farkl?? renkler
  
  return(list(simulated_data = expression_data, plot = expression_plot))
}

# Mod??l UI
RNASeqModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        numericInput(ns("n_genes"), "Number of Genes:", 100, min = 1),
        numericInput(ns("n_samples"), "Number of Samples:", 10, min = 1),
        numericInput(ns("noise_level"), "Noise Level:", 0.1, min = 0.01, step = 0.01),
        actionButton(ns("simulate"), "Start Simulation")
      ),
      mainPanel(
        plotOutput(ns("geneExpressionPlot")),
        verbatimTextOutput(ns("simulationSummary"))
      )
    )
  )
}

# Mod??l Server
RNASeqModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Sim??lasyon Fonksiyonu Sonucunu Saklamak i??in reaktif de??er
    simulation_result <- reactiveVal(NULL)
    
    observeEvent(input$simulate, {
      result <- simulateGeneExpression(
        n_genes = input$n_genes, 
        n_samples = input$n_samples, 
        noise_level = input$noise_level
      )
      simulation_result(result)
    })
    
    # Sim??lasyon Grafi??i
    output$geneExpressionPlot <- renderPlot({
      result <- simulation_result()
      if (!is.null(result)) {
        print(result$plot)
      }
    })
    
    # Sim??lasyon ??zeti
    output$simulationSummary <- renderPrint({
      result <- simulation_result()
      if (!is.null(result)) {
        summary(result$simulated_data)
      }
    })
  })
}
