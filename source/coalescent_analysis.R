# Gerekli K??t??phaneler
library(ape)
library(pegas)
library(ggplot2)
library(dplyr)
library(reshape2)

# === FONKS??YONLAR ===

# SNP Frekans Analizi
analyzeSNPFrequencies <- function(data) {
  if (!"SNP" %in% colnames(data) || !"Frequency" %in% colnames(data)) {
    stop("Data must include 'SNP' and 'Frequency' columns.")
  }
  
  diversity <- data %>%
    summarize(Mean_Frequency = mean(Frequency), SD_Frequency = sd(Frequency))
  
  return(diversity)
}

# Koalesans Zaman?? Tahmini
estimateCoalescentTime <- function(sequence_data) {
  # Genetik mesafe matrisini hesapla
  dist_matrix <- dist.dna(sequence_data, model = "K80")
  
  # Neighbor Joining A??a??
  nj_tree <- nj(dist_matrix)
  
  # Koalesans Zaman??
  coalescent_time <- sum(branching.times(nj_tree))
  
  return(list(tree = nj_tree, time = coalescent_time))
}

# Pop??lasyon Ge??mi??ini G??rselle??tirme
plotPopulationHistory <- function(tree) {
  ggtree(tree) +
    geom_tiplab() +
    labs(title = "Population History via Coalescent Analysis")
}

# === UI VE SERVER ===

coalescentUI <- function(id) {
  ns <- NS(id)
  
  tabsetPanel(
    tabPanel("SNP Frequency Analysis",
             fileInput(ns("snpFile"), "Upload SNP Data (e.g., .csv)"),
             actionButton(ns("analyzeSNP"), "Analyze SNP Frequencies"),
             tableOutput(ns("snpResults"))
    ),
    tabPanel("Coalescent Time Estimation",
             fileInput(ns("sequenceFile"), "Upload Sequence Data (FASTA format)"),
             actionButton(ns("estimateCoalescent"), "Estimate Coalescent Time"),
             tableOutput(ns("coalescentResults")),
             plotOutput(ns("coalescentTree"))
    )
  )
}

coalescentServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # SNP Frekans Analizi
    observeEvent(input$analyzeSNP, {
      req(input$snpFile)
      snp_data <- read.csv(input$snpFile$datapath)
      
      results <- tryCatch({
        analyzeSNPFrequencies(snp_data)
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        NULL
      })
      
      output$snpResults <- renderTable({
        results
      })
    })
    
    # Koalesans Zaman?? Tahmini
    observeEvent(input$estimateCoalescent, {
      req(input$sequenceFile)
      sequence_data <- read.dna(input$sequenceFile$datapath, format = "fasta")
      
      results <- tryCatch({
        estimateCoalescentTime(sequence_data)
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        NULL
      })
      
      if (!is.null(results)) {
        output$coalescentResults <- renderTable({
          data.frame(Coalescent_Time = results$time)
        })
        output$coalescentTree <- renderPlot({
          plotPopulationHistory(results$tree)
        })
      }
    })
  })
}
