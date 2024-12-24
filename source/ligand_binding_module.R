# Gerekli K??t??phaneler
library(shiny)
library(r3dmol)
library(DT)

# === FONKS??YONLAR ===

# Protein ve Ligand Y??kleme
loadProteinLigand <- function(protein_file, ligand_file) {
  protein <- readLines(protein_file)
  ligand <- readLines(ligand_file)
  list(protein = protein, ligand = ligand)
}

# Ligand Ba??lanma Sim??lasyonu (Basit bir model)
simulateLigandBinding <- function(protein, ligand) {
  # Basit ba??lanma sim??lasyonu ????kt??s??
  paste("Simulated binding interaction between protein and ligand.\n",
        "Protein: ", substr(protein, 1, 30), "...\n",
        "Ligand: ", substr(ligand, 1, 30), "...\n")
}

# Molek??ler Dinamik Sim??lasyonu i??in Haz??rl??k
prepareMDSimulation <- function(protein_file, ligand_file, output_path) {
  # Sim??lasyon i??in gerekli girdiler
  system(paste("cp", protein_file, output_path))
  system(paste("cp", ligand_file, output_path))
  paste("Simulation files prepared at:", output_path)
}

# === UI VE SERVER ===

ligandBindingUI <- function(id) {
  ns <- NS(id)
  
  tabsetPanel(
    tabPanel("Load Protein and Ligand",
             fileInput(ns("proteinFile"), "Upload Protein (.pdb)"),
             fileInput(ns("ligandFile"), "Upload Ligand (.pdb)"),
             actionButton(ns("loadFiles"), "Load Files"),
             verbatimTextOutput(ns("fileStatus")),
             r3dmolOutput(ns("structureViewer"))
    ),
    tabPanel("Ligand Binding Simulation",
             actionButton(ns("runBindingSim"), "Run Binding Simulation"),
             verbatimTextOutput(ns("bindingSimResult"))
    ),
    tabPanel("Molecular Dynamics Simulation",
             textInput(ns("outputPath"), "Enter Output Path for Simulation Files", value = "./simulation"),
             actionButton(ns("prepareMDSim"), "Prepare MD Simulation"),
             verbatimTextOutput(ns("mdSimResult"))
    )
  )
}

ligandBindingServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Protein ve Ligand Dosyalar??n?? Y??kleme
    protein_ligand <- reactiveVal(NULL)
    
    observeEvent(input$loadFiles, {
      req(input$proteinFile, input$ligandFile)
      protein_ligand_data <- tryCatch({
        loadProteinLigand(input$proteinFile$datapath, input$ligandFile$datapath)
      }, error = function(e) {
        return(paste("Error loading files:", e$message))
      })
      
      protein_ligand(protein_ligand_data)
      output$fileStatus <- renderText({
        if (is.null(protein_ligand_data)) {
          "Failed to load files."
        } else {
          paste("Files loaded successfully:\nProtein:", input$proteinFile$name, "\nLigand:", input$ligandFile$name)
        }
      })
      
      output$structureViewer <- renderR3dmol({
        r3dmol() %>%
          m_add_model(data = protein_ligand_data$protein, format = "pdb") %>%
          m_add_model(data = protein_ligand_data$ligand, format = "pdb") %>%
          m_set_style(style = m_style_cartoon()) %>%
          m_zoom_to()
      })
    })
    
    # Ligand Ba??lanma Sim??lasyonu
    observeEvent(input$runBindingSim, {
      req(protein_ligand())
      result <- tryCatch({
        simulateLigandBinding(protein_ligand()$protein, protein_ligand()$ligand)
      }, error = function(e) {
        return(paste("Error during simulation:", e$message))
      })
      
      output$bindingSimResult <- renderText({ result })
    })
    
    # Molek??ler Dinamik Sim??lasyonu Haz??rl??k
    observeEvent(input$prepareMDSim, {
      req(protein_ligand())
      result <- tryCatch({
        prepareMDSimulation(input$proteinFile$datapath, input$ligandFile$datapath, input$outputPath)
      }, error = function(e) {
        return(paste("Error preparing simulation:", e$message))
      })
      
      output$mdSimResult <- renderText({ result })
    })
  })
}
