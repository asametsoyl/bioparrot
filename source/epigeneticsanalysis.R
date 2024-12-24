# Gerekli K??t??phaneler
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(rmarkdown)
library(randomForest)

# === DNA METHYLATION ===

# CpG Adas?? Tespiti
findCpGIslands <- function(sequence) {
  matches <- gregexpr("CG", sequence)
  positions <- unlist(matches)
  return(data.frame(Position = positions))
}

# Diferansiyel Metilasyon Analizi
differentialMethylation <- function(data1, data2) {
  if (!all(c("Position", "Methylation_Level") %in% colnames(data1)) || 
      !all(c("Position", "Methylation_Level") %in% colnames(data2))) {
    stop("Data must include 'Position' and 'Methylation_Level' columns.")
  }
  
  combined_data <- merge(data1, data2, by = "Position", suffixes = c("_1", "_2"))
  combined_data$Difference <- combined_data$Methylation_Level_2 - combined_data$Methylation_Level_1
  
  return(combined_data)
}

# CpG ve Diferansiyel Analiz G??rselle??tirme
plotDifferentialMethylation <- function(data) {
  ggplot(data, aes(x = Position, y = Difference, color = Difference > 0)) +
    geom_point() +
    scale_color_manual(values = c("red", "blue")) +
    labs(title = "Differential Methylation", x = "Position", y = "Difference") +
    theme_minimal()
}

# === HISTONE MODIFICATIONS ===

# Yo??unluk Grafi??i
plotHistoneDensity <- function(data) {
  ggplot(data, aes(x = Position, fill = Modification)) +
    geom_density(alpha = 0.5) +
    labs(title = "Histone Modification Density", x = "Position", y = "Density") +
    theme_minimal()
}

# Modifikasyon ve Gen Kesi??imi
intersectHistoneRegions <- function(histone_data, gene_regions) {
  overlaps <- histone_data %>%
    filter(Position %in% gene_regions$Position)
  
  return(overlaps)
}

# === ATAC-SEQ VE ER??????LEB??L??RL??K ===

# ATAC-Seq Analizi
analyzeATACSeq <- function(data) {
  accessible_regions <- data %>% filter(Accessibility_Score > 0.5)
  return(accessible_regions)
}

# === MAK??NE ????REN??M?? ===

# CpG Metilasyon Tahmini
predictCpGMethylation <- function(data) {
  model <- randomForest(Methylation_Level ~ Position, data = data, ntree = 100)
  return(model)
}

# Hastal??k Risk Tahmini
predictDiseaseRisk <- function(data, model) {
  predictions <- predict(model, newdata = data)
  return(data.frame(Position = data$Position, Predicted_Risk = predictions))
}

# === RAPORLAMA ===

# Analiz Sonu??lar??n?? PDF Olarak Kaydetme
generateReport <- function(data, output_file) {
  rmarkdown::render(
    input = "report_template.Rmd",
    output_file = output_file,
    params = list(data = data)
  )
}

# === UI VE SERVER ===

epigeneticsUI <- function(id) {
  ns <- NS(id)
  
  tabsetPanel(
    tabPanel("DNA Methylation",
             fileInput(ns("methFile1"), "Upload Condition 1 (e.g., .csv)"),
             fileInput(ns("methFile2"), "Upload Condition 2 (e.g., .csv)"),
             actionButton(ns("analyzeMeth"), "Analyze Methylation"),
             tableOutput(ns("methResults")),
             plotOutput(ns("methPlot"))
    ),
    tabPanel("Histone Modifications",
             fileInput(ns("histoneFile"), "Upload Histone Modification Data (e.g., .csv)"),
             actionButton(ns("analyzeHistone"), "Analyze Histone Modifications"),
             plotOutput(ns("histonePlot"))
    ),
    tabPanel("ATAC-Seq",
             fileInput(ns("atacFile"), "Upload ATAC-Seq Data"),
             actionButton(ns("analyzeATAC"), "Analyze Accessibility"),
             tableOutput(ns("atacResults"))
    ),
    tabPanel("Machine Learning",
             fileInput(ns("mlFile"), "Upload CpG Data for Prediction"),
             actionButton(ns("trainModel"), "Train Model"),
             actionButton(ns("predictRisk"), "Predict Disease Risk"),
             tableOutput(ns("mlResults"))
    ),
    tabPanel("Reports",
             downloadButton(ns("downloadReport"), "Download PDF Report")
    )
  )
}

epigeneticsServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # DNA Metilasyon
    observeEvent(input$analyzeMeth, {
      req(input$methFile1, input$methFile2)
      meth_data1 <- read.csv(input$methFile1$datapath)
      meth_data2 <- read.csv(input$methFile2$datapath)
      diff_results <- differentialMethylation(meth_data1, meth_data2)
      
      output$methResults <- renderTable({ diff_results })
      output$methPlot <- renderPlot({ plotDifferentialMethylation(diff_results) })
    })
    
    # Histon Modifikasyonu
    observeEvent(input$analyzeHistone, {
      req(input$histoneFile)
      histone_data <- read.csv(input$histoneFile$datapath)
      output$histonePlot <- renderPlot({ plotHistoneDensity(histone_data) })
    })
    
    # ATAC-Seq
    observeEvent(input$analyzeATAC, {
      req(input$atacFile)
      atac_data <- read.csv(input$atacFile$datapath)
      accessible_regions <- analyzeATACSeq(atac_data)
      output$atacResults <- renderTable({ accessible_regions })
    })
    
    # Makine ????renimi
    observeEvent(input$trainModel, {
      req(input$mlFile)
      ml_data <- read.csv(input$mlFile$datapath)
      model <- predictCpGMethylation(ml_data)
      session$userData$model <- model
    })
    
    observeEvent(input$predictRisk, {
      req(input$mlFile, session$userData$model)
      ml_data <- read.csv(input$mlFile$datapath)
      predictions <- predictDiseaseRisk(ml_data, session$userData$model)
      output$mlResults <- renderTable({ predictions })
    })
  })
}
