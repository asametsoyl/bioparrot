# Gerekli K??t??phaneler
library(ggplot2)
library(dplyr)
library(Biostrings)
library(ComplexHeatmap)

# === SOMAT??K MUTASYON ANAL??Z?? ===

# Somatik Mutasyon Tespiti
analyzeSomaticMutations <- function(data) {
  if (!"Position" %in% colnames(data) || !"Type" %in% colnames(data)) {
    stop("Data must include 'Position' and 'Type' columns.")
  }
  
  mutation_freq <- data %>%
    group_by(Type) %>%
    summarize(Frequency = n())
  
  return(mutation_freq)
}

# Mutasyon Frekans Grafi??i
plotSomaticMutations <- function(data) {
  ggplot(data, aes(x = Type, y = Frequency, fill = Type)) +
    geom_bar(stat = "identity") +
    labs(title = "Somatic Mutation Frequencies", x = "Mutation Type", y = "Frequency") +
    theme_minimal()
}

# === CRISPR TASARIMI ===

# sgRNA Hedef B??lgesi Belirleme
designCRISPR <- function(sequence) {
  sgRNA_sites <- matchPattern("NGG", sequence)  # Cas9 PAM b??lgesi
  positions <- start(sgRNA_sites) - 20  # sgRNA hedefleme b??lgesi (20 n??kleotid ??ncesi)
  sgRNAs <- sapply(positions, function(pos) substr(sequence, pos, pos + 19))
  
  return(data.frame(Position = positions, sgRNA = sgRNAs))
}

# sgRNA ??zellik Analizi
analyzeCRISPR <- function(sgRNA_data) {
  sgRNA_data$GC_Content <- sapply(sgRNA_data$sgRNA, function(sg) {
    gc <- sum(str_count(sg, c("G", "C")))
    round(gc / nchar(sg) * 100, 2)
  })
  
  return(sgRNA_data)
}

# === KANSER MUTASYON ANAL??Z?? ===

# Kanser Mutasyonlar??n?? S??n??fland??rma
classifyCancerMutations <- function(data) {
  if (!"Gene" %in% colnames(data) || !"Clinical_Significance" %in% colnames(data)) {
    stop("Data must include 'Gene' and 'Clinical_Significance' columns.")
  }
  
  classification <- data %>%
    group_by(Clinical_Significance) %>%
    summarize(Frequency = n())
  
  return(classification)
}

# Kanser Mutasyonlar??n?? G??rselle??tirme
plotCancerMutations <- function(data) {
  ggplot(data, aes(x = Clinical_Significance, y = Frequency, fill = Clinical_Significance)) +
    geom_bar(stat = "identity") +
    labs(title = "Cancer Mutation Clinical Significance", x = "Clinical Significance", y = "Frequency") +
    theme_minimal()
}

# === UI VE SERVER ===

cancerCRISPRUI <- function(id) {
  ns <- NS(id)
  
  tabsetPanel(
    tabPanel("Somatic Mutations",
             fileInput(ns("somaticFile"), "Upload Somatic Mutation Data (e.g., .csv)"),
             actionButton(ns("analyzeSomatic"), "Analyze Somatic Mutations"),
             tableOutput(ns("somaticResults")),
             plotOutput(ns("somaticPlot"))
    ),
    tabPanel("CRISPR Design",
             textInput(ns("dnaSeq"), "Enter DNA Sequence (FASTA format)"),
             actionButton(ns("designCRISPR"), "Design sgRNAs"),
             tableOutput(ns("crisprResults"))
    ),
    tabPanel("Cancer Mutation Analysis",
             fileInput(ns("cancerFile"), "Upload Cancer Mutation Data (e.g., .csv)"),
             actionButton(ns("analyzeCancer"), "Analyze Cancer Mutations"),
             tableOutput(ns("cancerResults")),
             plotOutput(ns("cancerPlot"))
    )
  )
}

cancerCRISPRServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Somatic Mutations
    observeEvent(input$analyzeSomatic, {
      req(input$somaticFile)
      somatic_data <- read.csv(input$somaticFile$datapath)
      somatic_results <- analyzeSomaticMutations(somatic_data)
      
      output$somaticResults <- renderTable({ somatic_results })
      output$somaticPlot <- renderPlot({ plotSomaticMutations(somatic_results) })
    })
    
    # CRISPR Design
    observeEvent(input$designCRISPR, {
      req(input$dnaSeq)
      sequence <- DNAString(input$dnaSeq)
      crispr_results <- designCRISPR(sequence)
      
      output$crisprResults <- renderTable({ analyzeCRISPR(crispr_results) })
    })
    
    # Cancer Mutation Analysis
    observeEvent(input$analyzeCancer, {
      req(input$cancerFile)
      cancer_data <- read.csv(input$cancerFile$datapath)
      cancer_results <- classifyCancerMutations(cancer_data)
      
      output$cancerResults <- renderTable({ cancer_results })
      output$cancerPlot <- renderPlot({ plotCancerMutations(cancer_results) })
    })
  })
}
