# Gerekli K??t??phaneler
library(httr)
library(jsonlite)
library(r3dmol)
library(DT)

# === FONKS??YONLAR ===

# PDB Veri ??ndirme ve G??rselle??tirme
fetchPDBStructure <- function(pdb_id) {
  base_url <- "https://files.rcsb.org/download/"
  url <- paste0(base_url, pdb_id, ".pdb")
  response <- httr::GET(url)
  
  if (httr::status_code(response) != 200) {
    stop(paste("Error fetching PDB structure for ID:", pdb_id))
  }
  
  pdb_content <- httr::content(response, as = "text")
  return(pdb_content)
}

# OMIM Veritaban??ndan Veri ??ekme
fetchOMIMData <- function(omim_id, api_key) {
  base_url <- "https://api.omim.org/api/entry"
  url <- paste0(base_url, "?mimNumber=", omim_id, "&format=json&apiKey=", api_key)
  
  response <- httr::GET(url)
  
  if (httr::status_code(response) != 200) {
    stop(paste("Error fetching OMIM data for ID:", omim_id))
  }
  
  omim_data <- httr::content(response, as = "parsed", type = "application/json")
  return(omim_data)
}

# GenBank veya RefSeq Verisi ??ekme
fetchGenBankData <- function(accession_id) {
  base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  url <- paste0(base_url, "?db=nucleotide&id=", accession_id, "&rettype=fasta&retmode=text")
  
  response <- httr::GET(url)
  
  if (httr::status_code(response) != 200) {
    stop(paste("Error fetching data for accession ID:", accession_id))
  }
  
  sequence_data <- httr::content(response, as = "text")
  return(sequence_data)
}

# === UI VE SERVER ===

dataIntegrationUI <- function(id) {
  ns <- NS(id)
  
  tabsetPanel(
    tabPanel("PDB Viewer",
             textInput(ns("pdbID"), "Enter PDB ID:", value = "1CRN"),
             actionButton(ns("fetchPDB"), "Fetch PDB Structure"),
             verbatimTextOutput(ns("pdbStructure")),
             r3dmolOutput(ns("pdbViewer")),
             downloadButton(ns("downloadPDB"), "Download PDB Structure")
    ),
    tabPanel("OMIM Search",
             textInput(ns("omimID"), "Enter OMIM ID:", value = "100100"),
             textInput(ns("apiKey"), "Enter OMIM API Key:", placeholder = "Your OMIM API Key"),
             actionButton(ns("fetchOMIM"), "Fetch OMIM Data"),
             verbatimTextOutput(ns("omimData")),
             downloadButton(ns("downloadOMIM"), "Download OMIM Data")
    ),
    tabPanel("GenBank/RefSeq Fetch",
             textInput(ns("accessionID"), "Enter Accession ID:", value = "NC_000852"),
             actionButton(ns("fetchGenBank"), "Fetch GenBank/RefSeq Data"),
             verbatimTextOutput(ns("genBankData")),
             downloadButton(ns("downloadGenBank"), "Download GenBank Data")
    )
  )
}

dataIntegrationServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # PDB Viewer
    pdb_structure <- reactiveVal(NULL)
    
    observeEvent(input$fetchPDB, {
      req(input$pdbID)
      
      result <- tryCatch({
        fetchPDBStructure(input$pdbID)
      }, error = function(e) {
        return(paste("Error:", e$message))
      })
      
      pdb_structure(result)
      output$pdbStructure <- renderText({ result })
      
      output$pdbViewer <- renderR3dmol({
        r3dmol() %>%
          m_add_model(data = result, format = "pdb") %>%
          m_set_style(style = m_style_cartoon()) %>%
          m_zoom_to()
      })
    })
    
    output$downloadPDB <- downloadHandler(
      filename = function() {
        paste0(input$pdbID, ".pdb")
      },
      content = function(file) {
        writeLines(pdb_structure(), file)
      }
    )
    
    # OMIM Search
    omim_data <- reactiveVal(NULL)
    
    observeEvent(input$fetchOMIM, {
      req(input$omimID, input$apiKey)
      
      result <- tryCatch({
        fetchOMIMData(input$omimID, input$apiKey)
      }, error = function(e) {
        return(paste("Error:", e$message))
      })
      
      omim_data(result)
      output$omimData <- renderText({
        jsonlite::toJSON(result, pretty = TRUE)
      })
    })
    
    output$downloadOMIM <- downloadHandler(
      filename = function() {
        paste0(input$omimID, "_OMIM.json")
      },
      content = function(file) {
        write(jsonlite::toJSON(omim_data(), pretty = TRUE), file)
      }
    )
    
    # GenBank/RefSeq Fetch
    genbank_data <- reactiveVal(NULL)
    
    observeEvent(input$fetchGenBank, {
      req(input$accessionID)
      
      result <- tryCatch({
        fetchGenBankData(input$accessionID)
      }, error = function(e) {
        return(paste("Error:", e$message))
      })
      
      genbank_data(result)
      output$genBankData <- renderText({ result })
    })
    
    output$downloadGenBank <- downloadHandler(
      filename = function() {
        paste0(input$accessionID, ".fasta")
      },
      content = function(file) {
        writeLines(genbank_data(), file)
      }
    )
  })
}
