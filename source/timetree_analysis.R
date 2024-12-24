# Gerekli K??t??phaneler
library(httr)
library(jsonlite)
library(ggtree)
library(ape)

# === TIME TREE FONKS??YONLARI ===

# T??rler aras?? evrimsel zaman sorgusu
queryTimeTree <- function(species) {
  base_url <- "http://timetree.org/api/v5/timetree"
  
  # T??rleri birle??tir
  species_query <- paste(species, collapse = ",")
  url <- paste0(base_url, "?taxon=", URLencode(species_query))
  
  # API iste??i
  response <- GET(url)
  
  if (status_code(response) != 200) {
    stop("API request failed. Please check your species names or connection.")
  }
  
  # JSON verilerini i??leme
  data <- content(response, "parsed", simplifyVector = TRUE)
  
  # T??rler aras?? zamanlar?? d??nd??r
  time_data <- data.frame(
    Species1 = sapply(data$speciesPairs, function(x) x$taxon1),
    Species2 = sapply(data$speciesPairs, function(x) x$taxon2),
    DivergenceTime = sapply(data$speciesPairs, function(x) x$time)
  )
  
  return(time_data)
}

# Filogenetik a??ac?? olu??turma
buildPhylogeneticTree <- function(species, divergence_times) {
  tree <- nj(as.dist(divergence_times))
  tree$tip.label <- species
  return(tree)
}

# === UI VE SERVER ===

timeTreeUI <- function(id) {
  ns <- NS(id)
  
  tabsetPanel(
    tabPanel("Evolutionary Time Estimate",
             textInput(ns("speciesInput"), "Enter Species Names (comma-separated)", 
                       placeholder = "e.g., Homo sapiens, Pan troglodytes, Mus musculus"),
             actionButton(ns("queryTimeTree"), "Query TimeTree"),
             tableOutput(ns("timeResults")),
             plotOutput(ns("phyloTree")),
             downloadButton(ns("downloadResults"), "Download Results (CSV)")
    )
  )
}

timeTreeServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    time_results <- reactiveVal(NULL)
    
    # TimeTree sorgusu
    observeEvent(input$queryTimeTree, {
      req(input$speciesInput)
      species <- strsplit(input$speciesInput, ",")[[1]]
      species <- trimws(species)  # T??r adlar??ndaki gereksiz bo??luklar?? kald??r
      
      # T??rler aras?? evrimsel zaman?? sorgula
      results <- tryCatch({
        queryTimeTree(species)
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        NULL
      })
      
      if (is.null(results) || nrow(results) == 0) {
        showNotification("No results returned from TimeTree API.", type = "error")
        return()
      }
      
      time_results(results)
      output$timeResults <- renderTable({ results })
      
      # Filogenetik a??a?? olu??turma
      if (!is.null(results)) {
        # Mesafe matrisi olu??turma
        all_species <- unique(c(results$Species1, results$Species2))
        divergence_matrix <- matrix(NA, nrow = length(all_species), ncol = length(all_species))
        rownames(divergence_matrix) <- all_species
        colnames(divergence_matrix) <- all_species
        
        # Divergence Time de??erlerini matrise yerle??tir
        for (i in seq_len(nrow(results))) {
          sp1 <- results$Species1[i]
          sp2 <- results$Species2[i]
          time <- results$DivergenceTime[i]
          
          divergence_matrix[sp1, sp2] <- time
          divergence_matrix[sp2, sp1] <- time
        }
        
        divergence_matrix[is.na(divergence_matrix)] <- 0  # Bo?? de??erleri s??f??rla
        
        # Filogenetik a??a?? olu??tur
        tree <- nj(as.dist(divergence_matrix))
        
        output$phyloTree <- renderPlot({
          ggtree(tree) +
            geom_tiplab() +
            labs(title = "Phylogenetic Tree", subtitle = "TimeTree Estimated Evolutionary Relationships")
        })
      }
    })
    
    # Sonu??lar?? indirilebilir hale getirme
    output$downloadResults <- downloadHandler(
      filename = function() {
        paste0("TimeTree_Results_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(time_results(), file, row.names = FALSE)
      }
    )
  })
}
