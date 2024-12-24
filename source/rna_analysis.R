# Gerekli K??t??phaneler
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(seqinr)
library(ape)

# === miRNA ve lncRNA Analizleri ===

# miRNA ??fade Analizi
analyzeMiRNA <- function(data) {
  if (!"miRNA" %in% colnames(data) || !"Expression" %in% colnames(data)) {
    stop("Data must include 'miRNA' and 'Expression' columns.")
  }
  
  summary_stats <- data %>%
    group_by(miRNA) %>%
    summarize(Average_Expression = mean(Expression))
  
  return(summary_stats)
}

# miRNA G??rselle??tirme
plotMiRNA <- function(data) {
  ggplot(data, aes(x = reorder(miRNA, -Average_Expression), y = Average_Expression, fill = miRNA)) +
    geom_bar(stat = "identity") +
    labs(title = "miRNA Expression Levels", x = "miRNA", y = "Average Expression") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# === RNA D??zenleme B??lgeleri ===

# RNA D??zenleme B??lgelerini Tespit Et
findRNAEditingRegions <- function(data) {
  if (!"Position" %in% colnames(data) || !"Base_Change" %in% colnames(data)) {
    stop("Data must include 'Position' and 'Base_Change' columns.")
  }
  
  editing_regions <- data %>%
    filter(Base_Change %in% c("A-to-I", "C-to-U"))
  
  return(editing_regions)
}

# D??zenleme B??lgelerini G??rselle??tirme
plotRNAEditingRegions <- function(data) {
  ggplot(data, aes(x = Position, y = Base_Change, color = Base_Change)) +
    geom_point() +
    labs(title = "RNA Editing Regions", x = "Position", y = "Base Change") +
    theme_minimal()
}

# === Mutasyon ve Genom Tespiti ===

# SNP Analizi
analyzeSNPs <- function(data) {
  if (!"Position" %in% colnames(data) || !"Mutation_Type" %in% colnames(data)) {
    stop("Data must include 'Position' and 'Mutation_Type' columns.")
  }
  
  mutation_freq <- data %>%
    group_by(Mutation_Type) %>%
    summarize(Frequency = n())
  
  return(mutation_freq)
}

# SNP G??rselle??tirme
plotSNPs <- function(data) {
  ggplot(data, aes(x = Mutation_Type, y = Frequency, fill = Mutation_Type)) +
    geom_bar(stat = "identity") +
    labs(title = "SNP Mutation Types", x = "Mutation Type", y = "Frequency") +
    theme_minimal()
}

# Filogenetik A??a??
buildPhylogeneticTree <- function(data) {
  dist_matrix <- dist.dna(as.DNAbin(data), model = "K80")
  tree <- nj(dist_matrix)
  plot(tree, main = "Phylogenetic Tree")
}

# === UI VE SERVER ===

rnaAnalysisUI <- function(id) {
  ns <- NS(id)
  
  tabsetPanel(
    tabPanel("miRNA & lncRNA Analysis",
             fileInput(ns("mirnaFile"), "Upload miRNA Data (e.g., .csv)"),
             actionButton(ns("analyzeMiRNA"), "Analyze miRNA"),
             tableOutput(ns("mirnaResults")),
             plotOutput(ns("mirnaPlot"))
    ),
    tabPanel("RNA Editing",
             fileInput(ns("rnaEditFile"), "Upload RNA Editing Data (e.g., .csv)"),
             actionButton(ns("analyzeRNAEditing"), "Analyze RNA Editing"),
             tableOutput(ns("rnaEditResults")),
             plotOutput(ns("rnaEditPlot"))
    ),
    tabPanel("Mutation & Genome Analysis",
             fileInput(ns("snpFile"), "Upload SNP Data (e.g., .csv)"),
             actionButton(ns("analyzeSNPs"), "Analyze SNPs"),
             plotOutput(ns("snpPlot")),
             fileInput(ns("fastaFile"), "Upload Genome Data (FASTA format)"),
             actionButton(ns("buildTree"), "Build Phylogenetic Tree"),
             plotOutput(ns("phyloTree"))
    )
  )
}

rnaAnalysisServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # miRNA Analizi
    observeEvent(input$analyzeMiRNA, {
      req(input$mirnaFile)
      mirna_data <- read.csv(input$mirnaFile$datapath)
      mirna_results <- analyzeMiRNA(mirna_data)
      
      output$mirnaResults <- renderTable({ mirna_results })
      output$mirnaPlot <- renderPlot({ plotMiRNA(mirna_results) })
    })
    
    # RNA D??zenleme Analizi
    observeEvent(input$analyzeRNAEditing, {
      req(input$rnaEditFile)
      rna_edit_data <- read.csv(input$rnaEditFile$datapath)
      editing_results <- findRNAEditingRegions(rna_edit_data)
      
      output$rnaEditResults <- renderTable({ editing_results })
      output$rnaEditPlot <- renderPlot({ plotRNAEditingRegions(editing_results) })
    })
    
    # SNP Analizi
    observeEvent(input$analyzeSNPs, {
      req(input$snpFile)
      snp_data <- read.csv(input$snpFile$datapath)
      snp_results <- analyzeSNPs(snp_data)
      
      output$snpPlot <- renderPlot({ plotSNPs(snp_results) })
    })
    
    # Filogenetik A??a??
    observeEvent(input$buildTree, {
      req(input$fastaFile)
      genome_data <- readDNAStringSet(input$fastaFile$datapath)
      output$phyloTree <- renderPlot({ buildPhylogeneticTree(genome_data) })
    })
  })
}
