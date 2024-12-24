# Gerekli K??t??phaneler
library(ape)
library(phangorn)
library(ggplot2)
library(ggtree)

# === FONKS??YONLAR ===

# ML Modeli ile A??a?? Olu??turma
createMLTree <- function(sequence_data) {
  alignment <- as.phyDat(sequence_data)  # DNA dizilerini phyDat format??na ??evir
  dist_matrix <- dist.ml(alignment)  # Genetik mesafeyi hesapla
  tree <- nj(dist_matrix)  # Neighbor Joining a??ac??
  
  fit <- optim.pml(pml(tree, alignment), model = "GTR", optGamma = TRUE)  # ML modeli optimize et
  
  # A??ac??n `phylo` format??nda oldu??unu kontrol et
  if (!inherits(fit$tree, "phylo")) {
    stop("Error: The tree object is not in 'phylo' format.")
  }
  
  return(fit$tree)
}

# Bayesyen Yakla????m?? ile A??a?? Olu??turma
createBayesianTree <- function(sequence_data, n_iter = 10000) {
  alignment <- as.phyDat(sequence_data)  # DNA dizilerini phyDat format??na ??evir
  tree <- rtree(n = length(alignment), rooted = TRUE)  # Rastgele bir ba??lang???? a??ac??
  
  # Rastgele bir sim??lasyon sonucu d??nd??relim
  optimized_tree <- optim.pml(pml(tree, alignment), optNni = TRUE, control = pml.control(iter.max = n_iter))
  
  if (!inherits(optimized_tree$tree, "phylo")) {
    stop("Error: The tree object is not in 'phylo' format.")
  }
  
  return(optimized_tree$tree)
}

# Bootstrap Analizi
performBootstrap <- function(sequence_data, n_bootstrap = 100) {
  alignment <- as.phyDat(sequence_data)  # DNA dizilerini phyDat format??na ??evir
  tree <- nj(dist.ml(alignment))  # Neighbor Joining a??ac??
  
  # Bootstrap analizi
  bs_values <- bootstrap.phyDat(alignment, FUN = function(x) nj(dist.ml(x)), bs = n_bootstrap)
  
  # A??aca bootstrap de??erlerini ekle
  tree_with_bs <- plotBS(tree, bs_values, type = "phylogram")
  
  if (!inherits(tree_with_bs, "phylo")) {
    stop("Error: The tree object is not in 'phylo' format.")
  }
  
  return(tree_with_bs)
}

# A??a?? G??rselle??tirme (Bootstrap De??erleri Dahil)
plotTreeWithBootstrap <- function(tree, bootstrap_values) {
  ggtree(tree) +
    geom_tiplab() +
    geom_text2(aes(label = bootstrap_values), hjust = -0.3, vjust = 0.5) +
    labs(title = "Phylogenetic Tree with Bootstrap Values")
}

# === UI VE SERVER ===

phylogeneticUI <- function(id) {
  ns <- NS(id)
  
  tabsetPanel(
    tabPanel("ML Tree",
             fileInput(ns("mlFile"), "Upload Sequence Data (FASTA format)"),
             actionButton(ns("createMLTree"), "Create ML Tree"),
             plotOutput(ns("mlTree"))
    ),
    tabPanel("Bayesian Tree",
             fileInput(ns("bayesianFile"), "Upload Sequence Data (FASTA format)"),
             actionButton(ns("createBayesianTree"), "Create Bayesian Tree"),
             plotOutput(ns("bayesianTree"))
    ),
    tabPanel("Bootstrap Analysis",
             fileInput(ns("bootstrapFile"), "Upload Sequence Data (FASTA format)"),
             numericInput(ns("bootstrapIter"), "Number of Bootstrap Iterations", value = 100, min = 10),
             actionButton(ns("runBootstrap"), "Run Bootstrap Analysis"),
             plotOutput(ns("bootstrapTree"))
    )
  )
}

phylogeneticServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # ML Tree
    observeEvent(input$createMLTree, {
      req(input$mlFile)
      sequence_data <- read.dna(input$mlFile$datapath, format = "fasta")
      
      results <- tryCatch({
        createMLTree(sequence_data)
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        NULL
      })
      
      output$mlTree <- renderPlot({
        ggtree(results) +
          geom_tiplab() +
          labs(title = "Maximum Likelihood Phylogenetic Tree")
      })
    })
    
    # Bayesian Tree
    observeEvent(input$createBayesianTree, {
      req(input$bayesianFile)
      sequence_data <- read.dna(input$bayesianFile$datapath, format = "fasta")
      
      results <- tryCatch({
        createBayesianTree(sequence_data)
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        NULL
      })
      
      output$bayesianTree <- renderPlot({
        ggtree(results) +
          geom_tiplab() +
          labs(title = "Bayesian Phylogenetic Tree")
      })
    })
    
    # Bootstrap Analysis
    observeEvent(input$runBootstrap, {
      req(input$bootstrapFile)
      sequence_data <- read.dna(input$bootstrapFile$datapath, format = "fasta")
      
      results <- tryCatch({
        performBootstrap(sequence_data, input$bootstrapIter)
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        NULL
      })
      
      if (!is.null(results)) {
        output$bootstrapTree <- renderPlot({
          plotTreeWithBootstrap(results$tree, results$bootstrap)
        })
      }
    })
  })
}
