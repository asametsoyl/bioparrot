library(shiny)
library(Biostrings)
library(msa)
library(ape)
library(plotly)
library(manipulate)
library(ggtree)
library(rentrez)
library(randomForest)
library(stringr)
library(httr)
library(jsonlite)
library(KEGGREST)
library(pathview)
library(reshape2)
library(phangorn)
library(DT)
library(ggmsa)
library(r3dmol)
library(rsconnect)
library(dplyr)
library(seqinr)
library(genetics)
library(popbio)
source("CRISPRModule.R")
source("MetilasyonModule.R")
source("HistoneModule.R")
source("RNASeqModule.R")
source("GeneExpressionAnalysis.R")
source("DifferentialExpression.R")
source("SimulationModeling.R")
source("GeneticNetworkModule.R")
source("MutationSimulationModule.R")
source("ReportingModule.R")
source("3DStructureViewerModule.R")
source("AllelFrequencyModule.R")
source("HardyWeinbergModule.R")
source("auth_ui.R")
source("auth_server.R")
source("epigeneticsanalysis.R")
source("rna_analysis.R")
source("cancer_crispr_analysis.R")
source("timetree_analysis.R")
source("coalescent_analysis.R")
source("phylogenetic_tree_analysis.R")
source("data_integration_module.R")
source("ligand_binding_module.R")
source("disease_prediction_ppi_module.R")
source("de_novo_gene_annotation.R")


ui <- fluidPage(
  # Include Head Elements
  tags$head(
    tags$link(rel = "icon", type = "image/png", href = "https://i.hizliresim.com/orjupw0.png"),
    tags$style(HTML("      
      body {
        background: linear-gradient(45deg, #007bff, #ff0000, #ffa500, #ffcc00);
        font-family: 'Poppins', sans-serif;
        color: #333;
      }
      .title-panel {
        background-color: #ffcc00;
        color: #004c9c;
        padding: 15px;
        border-radius: 8px;
        display: flex;
        align-items: center;
        box-shadow: 0px 4px 12px rgba(0, 0, 0, 0.15);
      }
      .title-panel img {
        margin-right: 15px;
        border-radius: 50%;
        width: 60px;
      }
      .sidebar {
        background-color: #ffa500;
        border: 1px solid #e0e0e0;
        padding: 15px;
        border-radius: 8px;
        box-shadow: 0px 4px 12px rgba(0, 0, 0, 0.1);
        color: white;
      }
      .btn-primary {
        background-color: #ff0000;
        border-color: #e60000;
        color: white;
        border-radius: 20px;
        box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.2);
      }
      .btn-primary:hover {
        background-color: #e60000;
        border-color: #cc0000;
      }
      .tab-content {
        background-color: #ffffff;
        border: 1px solid #cccccc;
        padding: 15px;
        border-radius: 8px;
        box-shadow: 0px 4px 12px rgba(0, 0, 0, 0.1);
        transition: all 0.3s ease;
      }
      .file-input-container {
        background-color: #007bff;
        padding: 10px;
        border-radius: 8px;
        color: white;
      }
       .tab-content {
    animation: fadein 0.5s;
  }
  @keyframes fadein {
    from { opacity: 0; }
    to { opacity: 1; }
  }
    "))
  ),
  
  # Welcome Panel
  titlePanel(
    tags$div(
      class = "title-panel",
      tags$img(src = "https://i.hizliresim.com/orjupw0.png", height = "60px"),
      "BioParrot: DNA, RNA and Protein Analysis Application"
    )
  ),
  
  # Sidebar and Main Layout
  sidebarLayout(
    sidebarPanel(
      class = "sidebar",
      div(
        class = "file-input-container",
        fileInput("fileInput", "Drag and Drop or Upload Sequence Files", multiple = TRUE, accept = c(".fasta", ".txt", ".genbank", ".pdb"))
      ),
      textInput("manualInput", "Manual Sequence Input", placeholder = "Enter the sequence here"),
      selectInput("analysisType", "Select Analysis Type", choices = c("DNA", "RNA", "Protein")),
      textInput("motifInput", "Enter Motif to Search", placeholder = "Enter Motif"),
      actionButton("submit", "Start Analysis", class = "btn-primary"),
      textInput("searchQuery", "Search NCBI Database", placeholder = "Enter organism or sequence"),
      actionButton("searchNCBI", "Search NCBI", class = "btn-primary"),
      textInput("downloadID", "Enter NCBI ID to Download", placeholder = "Enter NCBI ID"),
      actionButton("downloadNCBI", "Download Sequence", class = "btn-primary"),
      actionButton("findSNPs", "Find SNPs", class = "btn-primary"),
      actionButton("compareButton", "Compare Sequences", class = "btn-primary"),
      actionButton("compareHeatmap", "Generate Difference Heatmap", class = "btn-primary"),
      actionButton("runBootstrap", "Run Bootstrap Analysis", class = "btn-primary"),
      actionButton("runNeighborNet", "Run NeighborNet Analysis", class = "btn-primary"),
      fileInput("snpFile", "Upload SNP Data File (.csv or .tsv):"),  # class arg??man?? kald??r??ld??
      actionButton("predictMutation", "Predict Mutation Impact", class = "btn-primary"),
      textInput("proteinID", "Enter Protein ID:", value = ""),
      actionButton("fetchProtein", "Fetch Protein Data", class = "btn-primary"),
      downloadButton("downloadProtein", "Download Protein Data", class = "btn-primary"),
      textInput("geneID", "Enter Gene ID:", value = ""),
      actionButton("fetchGene", "Fetch Gene Data", class = "btn-primary"),
      downloadButton("downloadGene", "Download Gene Data", class = "btn-primary"),
    ),
    mainPanel(
      tabsetPanel(
        # Welcome Tab
        tabPanel("Welcome",
                 tags$div(
                   class = "welcome-content",
                   tags$img(src = "https://i.hizliresim.com/orjupw0.png", height = "80px"),
                   tags$h3("Welcome to BioParrot!"),
                   tags$ul(
                     tags$li("Upload your sequence file or input a manual sequence."),
                     tags$li("Select the desired analysis type (DNA, RNA, Protein)."),
                     tags$li("Explore various genetic and bioinformatics tools available.")
                   )
                 )
        ),
        
        # Input Sequence
        tabPanel("Input Sequence", verbatimTextOutput("inputSequence"), class = "tab-content"),
        
        # Genetic Analysis
        tabPanel("Genetic Analysis",
                 tabsetPanel(
                   tabPanel("Transcription", verbatimTextOutput("transcription"), class = "tab-content"),
                   tabPanel("Translation", verbatimTextOutput("translation"), class = "tab-content"),
                   tabPanel("GC/AT Ratio", verbatimTextOutput("gcAtRatio"), class = "tab-content"),
                   tabPanel("Summary", h3("Sequence Summary"),textOutput("summary"), class = "tab-content"),
                   tabPanel("ORF Analysis", verbatimTextOutput("orfResults")),
                   tabPanel("Bootstrap Analysis", plotOutput("phyloPlot")),
                   tabPanel("NeighborNet Analysis", plotOutput("neighborNetPlot")),
                   tabPanel("CpG Islands",
                            fileInput("fileInputCpG", "Upload DNA Sequence (FASTA format)", accept = c(".fasta", ".txt")),
                            actionButton("analyzeCpG", "Analyze CpG Islands"),
                            tableOutput("cpgTable"),
                            plotOutput("cpgPlot")
                   ),
                   tabPanel("Histon Modification",
                     histoneModificationsUI("histone1"),
                   ),
                   tabPanel("CRISPR Analysis", 
                            CRISPRModuleUI("crispr1")),
                   tabPanel("Phylogenetic Tree", plotOutput("phyloTree", height = "800px"), class = "tab-content"))
        ),
        #biyokimya gruplama
        tabPanel("Biochemistry",
                 tabsetPanel(
                   tabPanel("Nucleotide Composition", verbatimTextOutput("nucleotideComp"), class = "tab-content"),
                   tabPanel("Sequence Quality Control", verbatimTextOutput("seqQC"), class = "tab-content"),
                   tabPanel("Motif Analysis", verbatimTextOutput("motifAnalysis"), class = "tab-content"),
                   tabPanel("Motif Visualization", plotlyOutput("motifPlot"), class = "tab-content"),
                   tabPanel("GC Content Visualization", plotlyOutput("gcPlot"), class = "tab-content"),
                   tabPanel("GC Content Heatmap",
                            actionButton("analyzeGCContent", "Generate GC Heatmap"),
                            plotOutput("gcHeatmap")
                   ),
                   tabPanel("Mutation Impact Prediction",
                            fileInput("snpFile", "Upload SNP Data (CSV format)", accept = c(".csv")),
                            actionButton("predictMutation", "Predict Mutation Impact", class = "btn-primary"),
                            tableOutput("mutationImpactTable"),
                            plotOutput("mutationImpactPlot")
                   ),
                   
                   tabPanel("Motif Heatmap",
                            actionButton("analyzeMotifHeatmap", "Generate Motif Heatmap"),
                            plotOutput("motifHeatmap")
                   ),
                   tabPanel("Allel Frequency Analysis", 
                            AllelFrequencyModuleUI("alleleFreq1")),
                   tabPanel("UniProt Fetch",
                            tabPanel("Protein Details", htmlOutput("proteinData")),
                            tabPanel("Amino Acid Composition", plotOutput("aaComposition"))),
                   tabPanel("Differential Expression", 
                            DifferentialExpressionUI("DE1")),
                   tabPanel("Ligand Binding Tools", ligandBindingUI("ligandBindingModule")),
                   tabPanel("Difference Heatmap", plotOutput("heatmapPlot"))
                 ), class = "tab-content"),
        
        tabPanel("Metagenomic Analysis",
                 fileInput("metagenomicFile", "Upload Sequence Abundance Data", accept = c(".csv", ".tsv")),
                 actionButton("analyzeAlpha", "Calculate Alpha Diversity", class = "btn-primary"),
                 verbatimTextOutput("alphaDiversity"),
                 actionButton("analyzeBeta", "Calculate Beta Diversity", class = "btn-primary"),
                 verbatimTextOutput("betaDiversity"),
                 plotOutput("betaDiversityPlot")
        ),
        tabPanel("Bioinformatics",
                 tabsetPanel(
                   tabPanel("Multiple Sequence Alignment",
                            h3("Alignment Results"),
                            verbatimTextOutput("msaResult"),
                            plotOutput("msaPlot"),
                            plotOutput("msaDendrogram")),
                   tabPanel("Base Pair Counting", verbatimTextOutput("bpCount"), class = "tab-content"),
                   tabPanel("NCBI Search Results", verbatimTextOutput("ncbiResults")),
                   tabPanel("Sequence Classification",
                            fileInput("sequenceFile", "Upload Sequence File (FASTA format)", accept = c(".fasta")),
                            actionButton("classifySequence", "Classify Sequences", class = "btn-primary"),
                            tableOutput("classificationResults"),
                            plotOutput("classificationPlot")
                   ),
                   tabPanel(textInput("sequenceInput", "Enter Sequence (DNA or Protein)", placeholder = "e.g., ATGCGATAC or MVLSPADKTNVKAAW"),
                            actionButton("analyzeFrequency", "Analyze Frequency"),
                            plotOutput("frequencyPlot"),
                            tableOutput("frequencyTable")
                   ),
                   tabPanel("Gene Details", tableOutput("geneData")),
                   tabPanel("Gene Description", uiOutput("geneDescription")),
                   
                
                     tabPanel("Mutation Table", tableOutput("mutationImpactTable")),
                     tabPanel("Mutation Plot", plotOutput("mutationImpactPlot")),
                   tabPanel("Download Status", verbatimTextOutput("downloadStatus")),
                   tabPanel("Reverse Complement", verbatimTextOutput("revComplement"), class = "tab-content"),
                 ), class = "tab-content"),
        # KEGG Analysis
        tabPanel("KEGG Analysis",
                 tabsetPanel(
                   tabPanel("KEGG ID Info", tableOutput("keggInfo")),
                   tabPanel("Pathway Visualization", plotOutput("pathwayPlot")),
                   tabPanel("Dynamic Pathway", plotOutput("livePathwayPlot"))
                 )
        ),
      
        
        
        
        # K-mer Analysis
        tabPanel("K-mer Analysis",
                 numericInput("kValue", "Enter k-mer length:", value = 3, min = 1, max = 10),
                 actionButton("analyzeKMer", "Analyze K-mer Frequency", class = "btn-primary"),
                 tableOutput("kMerTable"),
                 plotOutput("kMerBarPlot")
        ),
        tabPanel(
          "SNP Finder",
          h3("SNP Results"),
          verbatimTextOutput("snpResults")  # SNP sonu??lar??n?? burada g??sterece??iz
        ),
        
        tabPanel("Hardy-Weinberg Testi", 
                 HardyWeinbergUI("hwTestModule")),
        
        # Genetic Tools
        tabPanel("Genetic Tools",
                 tabsetPanel(
                   tabPanel("Genetic Network", 
                            GeneticNetworkModuleUI("geneticNetwork1")),
                   tabPanel("RNA-Seq Analysis", 
                            RNASeqModuleUI("rnaSeq1")
                   ),
                   tabPanel("Gene Expression Simulation", 
                            RNASeqModuleUI("geneExpressionSimulation")),
                   tabPanel("Epigenetics Analysis",
                            epigeneticsUI("epigenetics1")),
                   tabPanel("RNA Analysis", rnaAnalysisUI("rnaAnalysisModule")),
                   tabPanel("Mutation Simulation", MutationSimulationModuleUI("mutationSimulation")),
                   tabPanel("Cancer & CRISOR Tools", cancerCRISPRUI("cancerCRISPRModule")),
                   tabPanel("TimeTree Integration", timeTreeUI("timeTreeModule")),
                   tabPanel("Coalescent Analysis Tools", coalescentUI("coalescentModule")),
                   tabPanel("Phylogenetic Tree Tools", phylogeneticUI("phylogeneticModule")),
                   tabPanel("Data Integration Tools", dataIntegrationUI("dataIntegrationModule")),
                   tabPanel("Analysis Tool", diseasePredictionPPIUI("diseasePredictionPPIModule")),
                   tabPanel("Gen Annotation Tools", geneAnnotationUI("geneAnnotationModule")),
                   tabPanel("Reporting", ReportingModuleUI("reporting"))),
    
                 
        ),
        
        # Help
        tabPanel("Help",
                 HTML("<h4>How to Use:</h4>
            <ul>
              <li>Upload sequences in .fasta, .txt, or .genbank format.</li>
              <li>Enter DNA, RNA, or Protein sequences manually if preferred.</li>
              <li>Use 'Motif Analysis' to search for specific patterns in sequences.</li>
              <li>Generate visualizations for GC content and motifs in sequences.</li>
              <li>Perform multiple sequence alignments and phylogenetic analyses.</li>
            </ul>")
        )
      )
    )
  )
)

# Server: Server Logic
server <- function(input, output, session) {
  # Login functionality
  observeEvent(input$login, {
    showModal(modalDialog(
      title = "User Login",
      textInput("username", "Username", placeholder = "Enter your username"),
      passwordInput("password", "Password", placeholder = "Enter your password"),
      footer = tagList(
        actionButton("login_confirm", "Login", class = "btn-primary"),
        modalButton("Cancel")
      )
    ))
  })
  
  observeEvent(input$login_confirm, {
    if (input$username == "samesoyl" && input$password == "samet2003") {
      removeModal()
      showModal(modalDialog(title = "Success", "Login Successful!", easyClose = TRUE))
    } else {
      showModal(modalDialog(title = "Error", "Invalid Credentials!", easyClose = TRUE))
    }
  })
  # Validate input file format and size
  validateInputFile <- reactive({
    # Check if a file is uploaded
    if (is.null(input$fileInput)) return(TRUE)
    
    # Extract file extension
    ext <- tools::file_ext(input$fileInput$name)
    
    # Supported file extensions
    supported_formats <- c("fasta", "txt", "genbank")
    
    # Check file format
    if (!ext %in% supported_formats) {
      showModal(modalDialog(
        title = "Input Error",
        paste0(
          "Unsupported file format: .", ext, ". ",
          "Please upload one of the following formats: .fasta, .txt, .genbank."
        ),
        easyClose = TRUE,
        footer = NULL
      ))
      return(FALSE)
    }
    
    # Check file size (limit: 10 MB)
    file_info <- file.info(input$fileInput$datapath)
    max_size <- 10 * 1024 * 1024  # 10 MB in bytes
    if (file_info$size > max_size) {
      showModal(modalDialog(
        title = "Input Error",
        "The uploaded file exceeds the size limit of 10 MB. Please upload a smaller file.",
        easyClose = TRUE,
        footer = NULL
      ))
      return(FALSE)
    }
    
    # File validation passed
    TRUE
  })
  
  
  # Retrieve and validate user input
  # Retrieve and validate user input, ignoring unknown bases (e.g., N, R, Y)
  inputSequence <- reactive({
    # Validate input file
    if (!validateInputFile()) {
      return(DNAStringSet())
    }
    
    sequences <- DNAStringSet()
    
    # Process file input
    if (!is.null(input$fileInput)) {
      seqData <- tryCatch(
        {
          lapply(input$fileInput$datapath, readDNAStringSet)
        },
        error = function(e) {
          showModal(modalDialog(
            title = "File Read Error",
            "An error occurred while reading the file. Please check the file format and content.",
            easyClose = TRUE,
            footer = NULL
          ))
          NULL
        }
      )
      
      if (!is.null(seqData)) {
        sequences <- do.call(c, seqData)
      }
    }
    
    # Process manual input
    if (input$manualInput != "") {
      tryCatch(
        {
          sequences <- DNAStringSet(input$manualInput)
        },
        error = function(e) {
          showModal(modalDialog(
            title = "Input Error",
            "The manual input contains invalid characters. Please enter a valid DNA sequence (A, T, G, C, N).",
            easyClose = TRUE,
            footer = NULL
          ))
        }
      )
    }
    
    # Filter sequences to ignore unknown bases (e.g., N, R, Y)
    cleanSequences <- endoapply(sequences, function(seq) {
      # Replace unknown bases with empty strings
      DNAString(gsub("[^ATGCatgc]", "", as.character(seq)))
    })
    
    # Provide feedback if sequences were modified
    if (any(sapply(sequences, nchar) != sapply(cleanSequences, nchar))) {
      showModal(modalDialog(
        title = "Sequence Cleaned",
        "Some unknown bases (e.g., N, R, Y) were removed from the sequences.",
        easyClose = TRUE,
        footer = NULL
      ))
    }
    
    # Return cleaned sequences
    cleanSequences
  })
  

  
# Display sequence
output$inputSequence <- renderText({
  seq <- inputSequence()
  
  if (length(seq) == 0) {
    "No sequence input provided."
  } else {
    paste(sapply(names(seq), function(name) {
      paste0(">", name, "\n", as.character(seq[[name]]))
    }), collapse = "\n")
  }
})

# Sequence Quality Control
output$seqQC <- renderText({
  seq <- inputSequence()
  
  if (length(seq) == 0) {
    "No sequence provided for quality control."
  } else {
    # Count sequences with invalid characters
    invalidChars <- vapply(seq, function(x) {
      sum(nchar(gsub("[ATGCatgc]", "", as.character(x))))
    }, FUN.VALUE = numeric(1))
    
    # Identify invalid sequences
    invalidSequences <- names(seq)[invalidChars > 0]
    
    # Summary message
    totalSequences <- length(seq)
    validSequences <- sum(invalidChars == 0)
    invalidSequencesCount <- totalSequences - validSequences
    
    paste0(
      "Total sequences: ", totalSequences, "\n",
      "Valid sequences: ", validSequences, "\n",
      "Sequences with invalid characters: ", invalidSequencesCount, "\n",
      if (invalidSequencesCount > 0) {
        paste0("Invalid sequences: ", paste(invalidSequences, collapse = ", "))
      } else {
        ""
      }
    )
  }
})

# Highlight invalid bases in sequences
output$highlightInvalid <- renderText({
  seq <- inputSequence()
  
  if (length(seq) == 0) {
    "No sequence input provided."
  } else {
    paste(sapply(names(seq), function(name) {
      originalSeq <- as.character(seq[[name]])
      highlightedSeq <- gsub("[^ATGCatgc]", "<span style='color: red;'>\\0</span>", originalSeq)
      paste0(">", name, "\n", highlightedSeq)
    }), collapse = "\n")
  }
})

# Auto-correct invalid sequences
output$correctedSequences <- renderText({
  seq <- inputSequence()
  
  if (length(seq) == 0) {
    "No sequence input provided."
  } else {
    cleanedSeq <- endoapply(seq, function(s) {
      DNAString(gsub("[^ATGCatgc]", "", as.character(s)))
    })
    
    paste(sapply(names(cleanedSeq), function(name) {
      paste0(">", name, "\n", as.character(cleanedSeq[[name]]))
    }), collapse = "\n")
  }
})

# Base Pair Counting
output$bpCount <- renderText({
  seq <- inputSequence()
  
  if (length(seq) == 0) {
    "No sequence provided for base pair counting."
  } else {
    # Count individual bases
    baseCounts <- lapply(seq, function(x) {
      counts <- table(strsplit(as.character(x), split = "")[[1]])
      as.list(counts)
    })
    
    # Format output
    baseCountsFormatted <- sapply(names(seq), function(name) {
      counts <- baseCounts[[name]]
      paste0(name, ": ", paste(names(counts), counts, sep = "=", collapse = ", "))
    })
    
    paste("Base Pair Counts:\n", paste(baseCountsFormatted, collapse = "\n"))
  }
})

# Reverse Complement
output$revComplement <- renderText({
  seq <- inputSequence()
  
  if (length(seq) == 0) {
    "No sequence input provided."
  } else {
    # Calculate reverse complement
    revComp <- reverseComplement(seq)
    
    # Format output
    paste(sapply(names(revComp), function(name) {
      paste0(">", name, "\n", as.character(revComp[[name]]))
    }), collapse = "\n")
  }
})

  
# Nucleotide Composition
output$nucleotideComp <- renderText({
  seq <- inputSequence()
  
  if (length(seq) == 0) {
    "No sequence provided for nucleotide composition analysis."
  } else {
    # Calculate nucleotide composition for each sequence
    composition <- lapply(seq, function(s) {
      comp <- letterFrequency(DNAStringSet(s), letters = c("A", "T", "G", "C"), as.prob = TRUE) * 100
      round(comp, 2)
    })
    
    # Format output
    compositionFormatted <- sapply(names(seq), function(name) {
      paste0(name, ": ", paste(names(composition[[name]]), composition[[name]], sep = "=", collapse = ", "))
    })
    
    paste("Nucleotide Composition (%):\n", paste(compositionFormatted, collapse = "\n"))
  }
})

  
  # BLAST functionality placeholder (system call to local BLAST)
  observeEvent(input$submit, {
    seq <- inputSequence()
    if (length(seq) > 0) {
      writeXStringSet(seq, file = "query.fasta")
      system("blastn -query query.fasta -db nt -out results.txt")
      showModal(modalDialog(title = "BLAST Analysis", "BLAST analysis completed. Results saved to results.txt.", easyClose = TRUE))
    }
  })
  
  # Transcription process
  output$transcription <- renderText({
    seq <- inputSequence()
    
    if (is.null(seq) || length(seq) == 0) {
      return("No sequence input provided.")
    }
    
    tryCatch({
      # Ensure input is DNA before transcription
      dnaSeq <- DNAStringSet(seq)
      transcribed <- RNAStringSet(dnaSeq)
      paste(sapply(names(transcribed), function(name) {
        paste0(">", name, "\n", as.character(transcribed[[name]]))
      }), collapse = "\n")
    }, error = function(e) {
      paste("Error during transcription: ", e$message)
    })
  })
  
  output$translation <- renderText({
    seq <- inputSequence()
    
    if (is.null(seq) || length(seq) == 0) {
      return("No sequence input provided.")
    }
    
    tryCatch({
      # Convert input to uppercase and retain valid RNA characters
      seq <- toupper(seq)
      seq <- gsub("[^AUCG]", "", seq)
      
      # Check if the sequence is valid
      if (nchar(seq) == 0) {
        return("Invalid sequence: No valid RNA nucleotides (A, U, C, G) found.")
      }
      
      # Convert to RNAStringSet
      seq <- RNAStringSet(seq)
      
      # Perform translation
      translated <- translate(seq)
      translated <- as.character(translated)
      translated <- gsub("X", "*", translated)  # Replace invalid codons with '*'
      
      # Format output
      paste(">Translated Sequence\n", translated, collapse = "\n")
    }, error = function(e) {
      paste("Error during translation: ", e$message)
    })
  })
  
  
  output$summary <- renderText({
    seq <- inputSequence()
    
    if (is.null(seq) || length(seq) == 0) {
      return("No sequence input provided.")
    }
    
    tryCatch({
      # Convert sequence to RNA format and remove invalid characters
      seq <- toupper(seq)
      seq <- gsub("[^AUCG]", "", seq)
      
      if (!inherits(seq, "RNAStringSet")) {
        seq <- RNAStringSet(seq)
      }
      
      # Transcription and translation
      transcribed <- seq
      translated <- translate(seq)
      
      # Stop codon detection
      stopCodons <- lapply(as.character(transcribed), function(s) {
        positions <- gregexpr("UAA|UAG|UGA", s)[[1]]
        if (positions[1] != -1) {
          return(positions)
        } else {
          return("No stop codons")
        }
      })
      
      stopCodonSummary <- paste(sapply(seq_along(stopCodons), function(i) {
        positions <- stopCodons[[i]]
        name <- names(transcribed)[i]
        if (is.character(positions)) {
          paste0(">", name, ": ", positions)
        } else {
          paste0(">", name, ": Positions: ", paste(positions, collapse = ", "))
        }
      }), collapse = "\n")
      
      paste(
        "Summary:",
        "\n- Total input sequences: ", length(transcribed),
        "\n- Transcribed RNA sequences: ", length(transcribed),
        "\n- Translated protein sequences: ", length(translated),
        "\n- Stop Codons:\n", stopCodonSummary
      )
    }, error = function(e) {
      paste("Error during processing: ", e$message)
    })
  })
  
  
  
  # GC/AT ratio calculation
  output$gcAtRatio <- renderText({
    seq <- inputSequence()
    if (length(seq) == 0) {
      "No sequence input provided."
    } else {
      gcContent <- letterFrequency(seq, "GC")
      atContent <- letterFrequency(seq, "AT")
      gcRatio <- gcContent / (gcContent + atContent)
      atRatio <- atContent / (gcContent + atContent)
      paste("GC Ratio: ", gcRatio, "\nAT Ratio: ", atRatio)
    }
  })
  
  # Motif Visualization
output$motifPlot <- renderPlotly({
  seq <- inputSequence()
  motif <- input$motifInput
  
  # Input validation
  if (length(seq) == 0 || is.null(motif) || motif == "") {
    return(NULL)
  }
  
  tryCatch({
    # Data frame to store motif positions
    df <- data.frame(Sequence = character(), Position = numeric())
    
    # Find motif positions in each sequence
    for (i in seq_along(seq)) {
      matches <- matchPattern(motif, seq[[i]])
      if (length(matches) > 0) {
        df <- rbind(df, data.frame(
          Sequence = paste("Sequence", names(seq)[i]),
          Position = start(matches)
        ))
      }
    }
    
    # Handle no matches found
    if (nrow(df) == 0) {
      showModal(modalDialog(
        title = "No Matches Found",
        paste("The motif", motif, "was not found in any of the sequences."),
        easyClose = TRUE
      ))
      return(NULL)
    }
    
    # Create an interactive plot with Plotly
    plot_ly(df, x = ~Position, y = ~Sequence, type = 'scatter', mode = 'markers',
            marker = list(color = 'rgba(255, 99, 71, 0.8)', size = 10)) %>%
      layout(title = paste("Motif Positions: ", motif),
             xaxis = list(title = "Position"),
             yaxis = list(title = "Sequences"))
  }, error = function(e) {
    showModal(modalDialog(
      title = "Error",
      paste("An error occurred during motif visualization:", e$message),
      easyClose = TRUE
    ))
    NULL
  })
})

  
  # GC Content Visualization
output$gcPlot <- renderPlotly({
  seq <- inputSequence()
  
  if (length(seq) == 0) {
    return(NULL)
  }
  
  tryCatch({
    # Calculate GC content for each sequence
    gcData <- letterFrequency(seq, c("G", "C"), as.prob = FALSE)
    gcContent <- rowSums(gcData) / width(seq)
    
    # Create a data frame for plotting
    df <- data.frame(Sequence = names(seq), GCContent = round(gcContent, 2))
    
    # Create interactive bar plot
    plot_ly(data = df, x = ~Sequence, y = ~GCContent, type = 'bar',
            marker = list(color = df$GCContent, colorscale = 'Viridis')) %>%
      layout(title = "GC Content per Sequence",
             xaxis = list(title = "Sequence"),
             yaxis = list(title = "GC Content (Proportion)"))
  }, error = function(e) {
    showModal(modalDialog(
      title = "Error",
      paste("An error occurred during GC content visualization:", e$message),
      easyClose = TRUE
    ))
    NULL
  })
})

  
  output$phyloTree <- renderPlot({
  seq <- inputSequence()
  
  if (length(seq) == 0) {
    # No sequence data - Placeholder plot
    plot(1:10, rnorm(10), main = "No sequence data to create a phylogenetic tree")
  } else {
    tryCatch({
      # Multiple sequence alignment
      msaResult <- msa(seq)
      
      # Compute distance matrix
      distMatrix <- dist.dna(as.DNAbin(msaResult), model = "raw")
      
      # Construct phylogenetic tree
      phyloTree <- nj(distMatrix)
      
      # Update tip labels with shorter sequence names
      shortNames <- substr(names(seq), 33, 41)  # Adjust length if necessary
      phyloTree$tip.label <- shortNames
      
      # Plot the phylogenetic tree
      plot(phyloTree, main = "Phylogenetic Tree", cex = 0.7)
    }, error = function(e) {
      # Handle errors gracefully
      plot(1:10, rnorm(10), main = "Error in creating phylogenetic tree")
      text(5, 0, labels = paste("Error:", e$message), cex = 0.8, col = "red")
    })
  }
})

  
  # Motif Analysis
output$motifAnalysis <- renderText({
  seq <- inputSequence()
  motif <- input$motifInput
  
  # Validation: Check if sequences or motif are missing
  if (length(seq) == 0 || is.null(motif) || motif == "") {
    return("No valid sequences or motif provided.")
  }
  
  tryCatch({
    # Find motif locations in each sequence
    motifLocations <- lapply(seq, function(dna) {
      matchPattern(DNAString(motif), dna)
    })
    
    # Format results
    result <- sapply(seq_along(motifLocations), function(i) {
      matches <- motifLocations[[i]]
      if (length(matches) > 0) {
        positions <- start(matches)
        paste0("Sequence ", names(seq)[i], ": Motif found at positions: ", paste(positions, collapse = ", "))
      } else {
        paste0("Sequence ", names(seq)[i], ": Motif not found.")
      }
    })
    
    paste(result, collapse = "\n")
  }, error = function(e) {
    # Error handling
    paste("An error occurred during motif analysis:", e$message)
  })
})

  
# Multiple Sequence Alignment with Visualization
output$msaResult <- renderText({
  seq <- inputSequence()
  
  if (length(seq) < 2) {
    return("At least two sequences are required for alignment.")
  }
  
  tryCatch({
    # Perform multiple sequence alignment
    msaAligned <- msa(seq, method = "ClustalW")
    # Format aligned sequences as text
    paste(as.character(msaAligned), collapse = "\n")
  }, error = function(e) {
    paste("Error during alignment:", e$message)
  })
})

  
# Colorful Visualization of MSA
output$msaPlot <- renderPlot({
  seq <- inputSequence()
  
  if (length(seq) < 2) {
    ggplot() + labs(title = "At least two sequences are required for alignment") + theme_void()
  } else {
    tryCatch({
      msaAligned <- msa(seq, method = "ClustalW")
      alignedSequences <- as(msaAligned, "DNAStringSet")
      ggmsa::ggmsa(alignedSequences, font = NULL, color = "Clustal") +
        labs(title = "Sequence Alignment Visualization") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))  # X ekseni etiketlerini optimize et
    }, error = function(e) {
      ggplot() + labs(title = paste("Error during alignment visualization:", e$message)) + theme_void()
    })
  }
})


  
# Dendrogram for Multiple Sequence Alignment
output$msaDendrogram <- renderPlot({
  seq <- inputSequence()
  
  if (length(seq) < 2) {
    ggplot() + labs(title = "At least two sequences are required for dendrogram") + theme_void()
  } else {
    tryCatch({
      # Perform multiple sequence alignment
      msaAligned <- msa(seq, method = "ClustalW")
      # Compute distance matrix
      distMatrix <- dist.ml(as.phyDat(msaAligned, type = "DNA"))
      # Perform hierarchical clustering
      hc <- hclust(distMatrix, method = "average")
      # Plot dendrogram
      plot(hc, main = "Dendrogram for Sequence Alignment", xlab = "Sequences", sub = "")
    }, error = function(e) {
      ggplot() + labs(title = paste("Error during dendrogram creation:", e$message)) + theme_void()
    })
  }
})

  
  
  #ORF Bulma Fonksiyonu
  findORFs <- function(seq) {
    startCodon <- "ATG"
    stopCodons <- c("TAA", "TAG", "TGA")
    
    seqString <- as.character(seq)
    starts <- unlist(gregexpr(startCodon, seqString))
    stops <- unlist(lapply(stopCodons, function(codon) gregexpr(codon, seqString)))
    
    orfs <- list()
    for (start in starts) {
      stop <- min(stops[stops > start & (stops - start) %% 3 == 0], na.rm = TRUE)
      if (!is.na(stop)) {
        orfs[[paste0("ORF from ", start, " to ", stop)]] <- substr(seqString, start, stop + 2)
      }
    }
    
    if (length(orfs) == 0) {
      return("No ORFs found.")
    } else {
      return(orfs)
    }
  }
  # ORF sonuC'larD1nD1 gC6ster
  output$orfResults <- renderText({
    seq <- inputSequence()
    if (length(seq) == 0) {
      "No sequence input provided."
    } else {
      orfs <- findORFs(seq[[1]])
      if (is.character(orfs)) {
        orfs
      } else {
        paste(sapply(names(orfs), function(name) paste0(name, ":\n", orfs[[name]])), collapse = "\n\n")
      }
    }
  })
  # NCBI Arama Fonksiyonu
  searchNCBI <- function(query) {
    if (nzchar(query)) {
      # Sorgunun geC'erli olup olmadD1DD1nD1 kontrol et
      if (!grepl("^\\w+", query)) {
        return("Please enter a valid query!")
      }
      # NCBI veritabanD1nda sorgu iC'in arama yap
      search_result <- tryCatch(
        entrez_search(db = "nucleotide", term = query, retmax = 10),
        error = function(e) return(paste("Hata:", e$message))
      )
      if (!is.null(search_result) && length(search_result$ids) > 0) {
        # DetaylD1 bilgi almak iC'in her bir ID iC'in summary alD1nD1r
        summaries <- sapply(search_result$ids, function(id) {
          summary <- entrez_summary(db = "nucleotide", id = id)
          paste0("ID: ", id, "\nTitle: ", summary$title, 
                 "\nDescription: ", summary$source,
                 "\nOrganism: ", summary$organism, 
                 "\nGenomic Context: ", summary$genomic_context)
        })
        return(paste(summaries, collapse = "\n\n"))
      } else {
        return("No query found in NCBI database.")
      }
    } else {
      return("Please enter a query!")
    }
  }
  
  
  # NCBI DC6ndC<rme Fonksiyonu
  downloadNCBI <- function(id) {
    if (nzchar(id)) {
      tryCatch({
        fasta_data <- entrez_fetch(db = "nucleotide", id = id, rettype = "fasta")
        file_name <- paste0(id, ".fasta")
        write(fasta_data, file_name)
        return(paste("File downloaded successfully:", file_name))
      }, error = function(e) {
        return("Error: Unable to download the file. Please check the ID.")
      })
    } else {
      return("Please enter a valid NCBI ID.")
    }
  }
  # NCBI arama sonuC'larD1nD1 gC6ster
  observeEvent(input$searchNCBI, {
    output$ncbiResults <- renderText({
      searchNCBI(input$searchQuery)
    })
  })
  
  # NCBI'dan veri indirme iElemini gerC'ekleEtirir
  observeEvent(input$downloadNCBI, {
    output$downloadStatus <- renderText({
      downloadNCBI(input$downloadID)
    })
  })
  # SNP Bulma Fonksiyonu
  findSNPs <- function(seq) {
    snps <- list()
    seqString <- as.character(seq)
    for (i in 1:(nchar(seqString) - 1)) {
      if (substr(seqString, i, i) != substr(seqString, i + 1, i + 1)) {
        snps[[paste0("Position ", i)]] <- paste0(substr(seqString, i, i), " -> ", substr(seqString, i + 1, i + 1))
      }
    }
    
    if (length(snps) == 0) {
      return("No SNPs found.")
    } else {
      return(snps)
    }
  }
  # SNP sonuC'larD1nD1 gC6ster
  observeEvent(input$findSNPs, {
    output$snpResults <- renderText({
      seq <- inputSequence()
      if (length(seq) == 0) {
        "No sequence input provided."
      } else {
        snps <- findSNPs(seq[[1]])
        if (is.character(snps)) {
          snps
        } else {
          paste(sapply(names(snps), function(name) paste0(name, ": ", snps[[name]])), collapse = "\n")
        }
      }
    })
  })
  
  # Hata YC6netimi ve KullanD1cD1 Geri Bildirimi
  observe({
    if (input$manualInput != "" && !grepl("^[ATGCatgc]+$", input$manualInput)) {
      showModal(modalDialog(
        title = "Input Error",
        "Invalid characters detected in manual sequence input. Please enter valid DNA sequences only.",
        easyClose = TRUE,
        footer = NULL
      ))
    }
  })
  
  calculateGCRatio <- function(seq) {
    seq <- DNAStringSet(seq)
    gcContent <- rowSums(letterFrequency(seq, c("G", "C"), as.prob = TRUE)) * 100
    return(gcContent)
  }
  
  
  
  observeEvent(input$analyzeGCContent, {
  seq <- tryCatch({
    inputSequence()  # Kullan??c??n??n DNA dizilerini al
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(seq) || length(seq) == 0) {
    # E??er dizi yoksa bo?? bir grafik g??ster
    output$gcHeatmap <- renderPlot({
      ggplot() + 
        labs(title = "No sequences provided for GC content analysis") +
        theme_void()
    })
    return()
  }
  
  gcContent <- tryCatch({
    calculateGCRatio(seq)  # GC oran??n?? hesapla
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(gcContent)) {
    # E??er GC i??eri??i hesaplanam??yorsa hata mesaj?? g??ster
    output$gcHeatmap <- renderPlot({
      ggplot() +
        labs(title = "Error calculating GC content") +
        theme_void()
    })
    return()
  }
  
  # Veriyi bir veri ??er??evesine d??n????t??r
  gcData <- data.frame(
    Sequence = names(seq),
    GC_Content = gcContent
  )
  
  # Is?? haritas??n?? olu??tur ve g??ster
  output$gcHeatmap <- renderPlot({
    ggplot(gcData, aes(x = Sequence, y = 1, fill = GC_Content)) +
      geom_tile(color = "white") +
      scale_fill_gradient(low = "blue", high = "red") +
      labs(title = "GC Content Heatmap", x = "Sequence", y = "", fill = "GC (%)") +
      theme_minimal() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank())
  })
})

  
  observeEvent(input$analyzeMotifHeatmap, {
    seq <- inputSequence()
    motif <- input$motifInput
    
    if (length(seq) == 0 || motif == "") {
      output$motifHeatmap <- renderPlot({
        ggplot() +
          labs(title = "No sequences or motif provided for analysis") +
          theme_void()
      })
      return()
    }
    
    # Motif e??le??melerini bul
    motifMatches <- lapply(seq, function(dna) {
      match <- matchPattern(DNAString(motif), dna)
      if (length(match) > 0) {
        return(start(match))
      } else {
        return(integer(0))
      }
    })
    
    # Veriyi g??rselle??tirme i??in d??zenle
    motifData <- data.frame(
      Sequence = rep(names(seq), sapply(motifMatches, length)),
      Position = unlist(motifMatches)
    )
    
    if (nrow(motifData) == 0) {
      output$motifHeatmap <- renderPlot({
        ggplot() +
          labs(title = "No motif matches found in the sequences") +
          theme_void()
      })
      return()
    }
    
    # Motif ??s?? haritas??
    output$motifHeatmap <- renderPlot({
      ggplot(motifData, aes(x = Position, y = Sequence, fill = Position)) +
        geom_tile(color = "white") +
        scale_fill_gradient(low = "green", high = "yellow") +
        labs(title = "Motif Positions Heatmap", x = "Position", y = "Sequence", fill = "Position") +
        theme_minimal() +
        theme(axis.text.y = element_text(angle = 0, hjust = 1))
    })
  })
  #kmer analysis
  calculateKMerFrequency <- function(seq, k) {
    kmers <- oligonucleotideFrequency(seq, width = k)
    kmers_df <- as.data.frame(kmers)
    kmers_df$Kmer <- rownames(kmers_df)  # K-mer isimlerini ekle
    return(kmers_df)
  }
  
  
  
  # K-mer Analizi Event'i
  observeEvent(input$analyzeKMer, {
    seq <- inputSequence()  # Kullan??c??n??n DNA dizisi girdisi
    k <- input$kValue       # Kullan??c??n??n K-mer uzunlu??u girdisi
    
    if (length(seq) == 0) {
      output$kMerTable <- renderTable({
        data.frame(Message = "No sequence provided.")
      })
      output$kMerBarPlot <- renderPlot({
        ggplot() + labs(title = "No sequence provided for K-mer analysis") + theme_void()
      })
      return()
    }
    
    tryCatch({
      # K-mer frekanslar??n?? hesaplama
      kmerData <- calculateKMerFrequency(seq, k)
      colnames(kmerData)[1] <- "Frequency"  # S??tun ismini d??zelt
      
      # Tablo olarak g??ster
      output$kMerTable <- renderTable({
        kmerData
      }, rownames = TRUE)
      
      # ??ubuk grafikle g??rselle??tir
      output$kMerBarPlot <- renderPlot({
        ggplot(kmerData, aes(x = reorder(Kmer, -Frequency), y = Frequency)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          theme_minimal() +
          labs(title = paste0("K-mer Frequency (k=", k, ")"), x = "K-mer", y = "Frequency") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
      })
    }, error = function(e) {
      # Hata durumunda mesaj g??ster
      output$kMerTable <- renderTable({
        data.frame(Message = paste("Error:", e$message))
      })
      output$kMerBarPlot <- renderPlot({
        ggplot() + labs(title = "Error in K-mer analysis") + theme_void()
      })
    })
  })
  
  loadGenomeData <- function(filepath) {
    file_ext <- tools::file_ext(filepath)
    switch(file_ext,
           "bed" = rtracklayer::import(filepath, format = "BED"),
           "bam" = rtracklayer::import(filepath, format = "BAM"),
           "vcf" = VariantAnnotation::readVcf(filepath),
           "gff" = rtracklayer::import(filepath, format = "GFF"),
           "gtf" = rtracklayer::import(filepath, format = "GTF"),
           stop("Unsupported file format. Please upload .bed, .bam, .vcf, .gff, or .gtf files.")
    )
  }
  
  renderGenomeBrowser <- function(genomeData, region) {
    region <- tryCatch({
      rtracklayer::GRanges(region)
    }, error = function(e) {
      stop("Invalid genomic region format. Use 'chr:start-end'.")
    })
    
    autoplot(genomeData, which = region)
  }
  
  observeEvent(input$viewGenome, {
    req(input$genomeFile)
    
    genomeData <- tryCatch({
      loadGenomeData(input$genomeFile$datapath)
    }, error = function(e) {
      output$genomePlot <- renderPlot({
        ggplot() + labs(title = paste("Error:", e$message)) + theme_void()
      })
      return(NULL)
    })
    
    region <- input$regionInput
    output$genomePlot <- renderPlot({
      tryCatch({
        renderGenomeBrowser(genomeData, region)
      }, error = function(e) {
        ggplot() + labs(title = paste("Error:", e$message)) + theme_void()
      })
    })
  })
  
  calculateAlphaDiversity <- function(data) {
  vegan::diversity(data, index = "shannon")
}

calculateBetaDiversity <- function(data) {
  vegan::vegdist(data, method = "bray")
}

observeEvent(input$analyzeAlpha, {
  req(input$metagenomicFile)
  
  data <- tryCatch({
    read.csv(input$metagenomicFile$datapath, row.names = 1)
  }, error = function(e) {
    output$alphaDiversity <- renderText({
      paste("Error loading data:", e$message)
    })
    return(NULL)
  })
  
  output$alphaDiversity <- renderText({
    tryCatch({
      alpha <- calculateAlphaDiversity(data)
      paste("Alpha Diversity (Shannon Index):", paste(alpha, collapse = ", "))
    }, error = function(e) {
      paste("Error calculating alpha diversity:", e$message)
    })
  })
})

observeEvent(input$analyzeBeta, {
  req(input$metagenomicFile)
  
  data <- tryCatch({
    read.csv(input$metagenomicFile$datapath, row.names = 1)
  }, error = function(e) {
    output$betaDiversity <- renderText({
      paste("Error loading data:", e$message)
    })
    return(NULL)
  })
  
  output$betaDiversity <- renderText({
    tryCatch({
      beta <- calculateBetaDiversity(data)
      paste("Beta Diversity (Bray-Curtis Distance):", beta)
    }, error = function(e) {
      paste("Error calculating beta diversity:", e$message)
    })
  })
  
  output$betaDiversityPlot <- renderPlot({
    tryCatch({
      beta <- calculateBetaDiversity(data)
      heatmap(as.matrix(beta), main = "Beta Diversity (Bray-Curtis)", col = colorRampPalette(c("blue", "white", "red"))(100))
    }, error = function(e) {
      ggplot() + labs(title = "Error in beta diversity visualization") + theme_void()
    })
  })
})

# SNP Verisi Y??kleme Fonksiyonu
loadSNPData <- function(filepath) {
  file_ext <- tools::file_ext(filepath)
  if (file_ext == "csv") {
    read.csv(filepath, header = TRUE)
  } else if (file_ext == "tsv") {
    read.delim(filepath, header = TRUE)
  } else {
    stop("Unsupported file format. Please upload a .csv or .tsv file.")
  }
}

# Mutasyon Etki Tahmini Fonksiyonu
predictMutationImpact <- function(snpData) {
  if (!"SNP" %in% colnames(snpData)) {
    stop("The SNP column is missing in the input data.")
  }
  if (!"Effect" %in% colnames(snpData)) {
    stop("The Effect column is missing in the input data.")
  }
  
  snpData$Impact <- ifelse(grepl("nonsense|missense", snpData$Effect, ignore.case = TRUE), 
                           "High",
                           ifelse(grepl("synonymous", snpData$Effect, ignore.case = TRUE), 
                                  "Low", "Moderate"))
  return(snpData)
}

# SNP Etki Tahmini Event'i
observeEvent(input$predictMutation, {
  req(input$snpFile)
  
  # SNP Verisini Y??kle
  snpData <- tryCatch({
    loadSNPData(input$snpFile$datapath)
  }, error = function(e) {
    output$mutationImpactTable <- renderTable({
      data.frame(Message = paste("Error loading SNP data:", e$message))
    })
    return(NULL)
  })
  
  if (is.null(snpData)) return()
  
  # Etki Tahmini Yap
  predictedData <- tryCatch({
    predictMutationImpact(snpData)
  }, error = function(e) {
    output$mutationImpactTable <- renderTable({
      data.frame(Message = paste("Error in mutation impact prediction:", e$message))
    })
    return(NULL)
  })
  
  if (is.null(predictedData)) return()
  
  # Sonu??lar?? Tablo Olarak G??ster
  output$mutationImpactTable <- renderTable({
    predictedData
  })
  
  # Sonu??lar?? G??rselle??tir
  output$mutationImpactPlot <- renderPlot({
    ggplot(predictedData, aes(x = Impact, fill = Impact)) +
      geom_bar() +
      theme_minimal() +
      labs(title = "Mutation Impact Distribution", x = "Impact", y = "Count") +
      scale_fill_brewer(palette = "Set3")
  })
})

  # Sekanslar?? Y??kleme Fonksiyonu
  loadSequences <- function(filepath) {
    sequences <- readDNAStringSet(filepath)
    data.frame(Name = names(sequences), Sequence = as.character(sequences))
  }
  
  # ??zellik ????kar??m?? Fonksiyonu
  extractFeatures <- function(sequenceData) {
    features <- data.frame(
      GC_Content = sapply(sequenceData$Sequence, function(seq) {
        gc <- sum(str_count(seq, c("G", "C")))
        gc / nchar(seq)
      }),
      Length = nchar(sequenceData$Sequence)
    )
    return(features)
  }
  
  # Model E??itimi Fonksiyonu
  trainModel <- function(trainingData, labels) {
    model <- randomForest(trainingData, as.factor(labels), ntree = 100)
    return(model)
  }
  
  # S??n??fland??rma Fonksiyonu
  classifySequences <- function(features, model) {
    predictions <- predict(model, newdata = features)
    return(predictions)
  }
  
  # Server Taraf??nda ????lem
  observeEvent(input$classifySequence, {
    req(input$sequenceFile)
    
    # Sekanslar?? Y??kle
    sequenceData <- tryCatch({
      loadSequences(input$sequenceFile$datapath)
    }, error = function(e) {
      output$classificationResults <- renderTable({
        data.frame(Message = paste("Error loading sequences:", e$message))
      })
      return(NULL)
    })
    
    if (is.null(sequenceData)) return()
    
    # ??zellik ????kar??m??
    features <- tryCatch({
      extractFeatures(sequenceData)
    }, error = function(e) {
      output$classificationResults <- renderTable({
        data.frame(Message = paste("Error extracting features:", e$message))
      })
      return(NULL)
    })
    
    if (is.null(features)) return()
    
    # ??rnek Model E??itimi
    model <- tryCatch({
      trainModel(features, sample(c("Class1", "Class2"), nrow(features), replace = TRUE))
    }, error = function(e) {
      output$classificationResults <- renderTable({
        data.frame(Message = paste("Error training/loading model:", e$message))
      })
      return(NULL)
    })
    
    if (is.null(model)) return()
    
    # S??n??fland??rma
    predictions <- tryCatch({
      classifySequences(features, model)
    }, error = function(e) {
      output$classificationResults <- renderTable({
        data.frame(Message = paste("Error classifying sequences:", e$message))
      })
      return(NULL)
    })
    
    if (is.null(predictions)) {
      output$classificationPlot <- renderPlot({
        ggplot() +
          labs(title = "No predictions available") +
          theme_void()
      })
      return()
    }
    
    sequenceData$Prediction <- predictions
    
    # Sonu??lar?? G??ster
    output$classificationResults <- renderTable({
      sequenceData
    })
    
    # G??rselle??tirme
    output$classificationPlot <- renderPlot({
      ggplot(sequenceData, aes(x = Prediction, fill = Prediction)) +
        geom_bar() +
        theme_minimal() +
        labs(title = "Sequence Classification Results", x = "Class", y = "Count") +
        scale_fill_brewer(palette = "Set3")
    })
  })
  # UniProt Veri ??ekme Fonksiyonu
  fetchUniProtData <- function(protein_id) {
    base_url <- "https://rest.uniprot.org/uniprotkb"
    url <- paste0(base_url, "/", URLencode(protein_id), "?format=json")
    
    response <- httr::GET(url)
    
    if (httr::status_code(response) != 200) {
      return(data.frame(Message = paste("Error fetching data for protein ID:", protein_id)))
    }
    
    protein_data <- httr::content(response, as = "parsed", type = "application/json")
    
    # ??nemli bilgileri ay??kla
    name <- ifelse(!is.null(protein_data$proteinDescription$recommendedName$fullName$value),
                   protein_data$proteinDescription$recommendedName$fullName$value,
                   "N/A")
    sequence <- ifelse(!is.null(protein_data$sequence$value),
                       as.character(protein_data$sequence$value),
                       "N/A")
    organism <- ifelse(!is.null(protein_data$organism$scientificName),
                       protein_data$organism$scientificName,
                       "N/A")
    
    # Sekans?? her 15 karakterde bir alt sat??ra b??l
    formatted_sequence <- paste(strwrap(sequence, width = 15), collapse = "<br>")
    
    return(data.frame(
      Protein_ID = protein_id,
      Name = name,
      Organism = organism,
      Sequence = formatted_sequence
    ))
  }
  
  # Amino Asit Analizi
  analyzeSequence <- function(sequence) {
    aa_counts <- seqinr::count(s2c(sequence), wordsize = 1)  # Amino asit frekanslar??
    aa_data <- data.frame(AminoAcid = names(aa_counts), Count = unlist(aa_counts))
    return(aa_data)
  }
  
  # UniProt Verisi ve Analiz
  observeEvent(input$fetchProtein, {
    req(input$proteinID)
    
    # Veri ??ekme
    protein_data <- tryCatch({
      fetchUniProtData(input$proteinID)
    }, error = function(e) {
      return(data.frame(Message = paste("Error:", e$message)))
    })
    
    # Protein Detaylar??n?? G??ster
    output$proteinData <- renderUI({
      if ("Message" %in% names(protein_data)) {
        HTML(paste("<strong>Error:</strong>", protein_data$Message))
      } else {
        HTML(paste(
          "<strong>Protein ID:</strong>", protein_data$Protein_ID[1], "<br>",
          "<strong>Name:</strong>", protein_data$Name[1], "<br>",
          "<strong>Organism:</strong>", protein_data$Organism[1], "<br>",
          "<strong>Sequence:</strong><br><pre>", protein_data$Sequence[1], "</pre>"
        ))
      }
    })
    
    # Amino Asit Analizi ve G??rselle??tirme
    if (!is.null(protein_data$Sequence) && protein_data$Sequence != "N/A") {
      sequence <- gsub("<br>", "", protein_data$Sequence[1])  # HTML'den temizlenmi?? dizi
      aa_analysis <- analyzeSequence(sequence)
      
      # Amino Asit Kompozisyon Grafi??i
      output$aaComposition <- renderPlot({
        ggplot(aa_analysis, aes(x = AminoAcid, y = Count, fill = AminoAcid)) +
          geom_bar(stat = "identity") +
          theme_minimal() +
          labs(title = "Amino Acid Composition", x = "Amino Acid", y = "Count") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
      })
      
      # Veri ??ndirilebilirlik
      output$downloadProtein <- downloadHandler(
        filename = function() {
          paste0("protein_", input$proteinID, "_aa_composition.csv")
        },
        content = function(file) {
          write.csv(aa_analysis, file, row.names = FALSE)
        }
      )
    }
  })
  # Ensembl Veri ??ekme Fonksiyonu
  fetchEnsemblData <- function(gene_id) {
    base_url <- "https://rest.ensembl.org"
    url <- paste0(base_url, "/lookup/id/", URLencode(gene_id), "?content-type=application/json")
    
    response <- GET(url)
    
    if (status_code(response) != 200) {
      return(data.frame(
        Message = paste("Error fetching data for gene ID:", gene_id, 
                        "- HTTP Status Code:", status_code(response))
      ))
    }
    
    gene_data <- content(response, as = "parsed", type = "application/json")
    
    # Gerekli bilgileri ay??kla
    return(data.frame(
      Gene_ID = gene_data$id,
      Name = gene_data$display_name,
      Species = gene_data$species,
      Start = gene_data$start,
      End = gene_data$end,
      Strand = ifelse(gene_data$strand == 1, "Forward", "Reverse"),
      Description = gene_data$description
    ))
  }
  
  # Ensembl Veri Sorgusu
  observeEvent(input$fetchGene, {
    req(input$geneID)
    
    gene_data <- tryCatch({
      fetchEnsemblData(input$geneID)
    }, error = function(e) {
      data.frame(Message = paste("Error:", e$message))
    })
    
    # Gen Detaylar??n?? G??ster
    output$geneData <- renderTable({
      gene_data
    })
    
    # Gen A????klamas??n?? G??ster
    output$geneDescription <- renderUI({
      if ("Message" %in% names(gene_data)) {
        HTML(paste("<strong>Error:</strong>", gene_data$Message))
      } else {
        HTML(paste(
          "<strong>Gene ID:</strong>", gene_data$Gene_ID[1], "<br>",
          "<strong>Name:</strong>", gene_data$Name[1], "<br>",
          "<strong>Species:</strong>", gene_data$Species[1], "<br>",
          "<strong>Start Position:</strong>", gene_data$Start[1], "<br>",
          "<strong>End Position:</strong>", gene_data$End[1], "<br>",
          "<strong>Strand:</strong>", gene_data$Strand[1], "<br>",
          "<strong>Description:</strong><br>", gene_data$Description[1]
        ))
      }
    })
    
    # Veri ??ndirme
    output$downloadGene <- downloadHandler(
      filename = function() {
        paste0("gene_", input$geneID, ".csv")
      },
      content = function(file) {
        write.csv(gene_data, file, row.names = FALSE)
      }
    )
  })
  # KEGG ID Sorgulama Fonksiyonu
  fetchKEGGID <- function(gene_name) {
    tryCatch({
      kegg_result <- keggFind("genes", gene_name)
      
      if (length(kegg_result) == 0) {
        return(data.frame(Message = paste("No KEGG ID found for gene:", gene_name)))
      }
      
      return(data.frame(
        Gene = gene_name,
        KEGG_ID = names(kegg_result),
        Description = kegg_result
      ))
    }, error = function(e) {
      return(data.frame(Message = paste("Error fetching KEGG ID:", e$message)))
    })
  }
  
  # KEGG Yol Haritas?? G??rselle??tirme Fonksiyonu
  visualizePathway <- function(kegg_id, pathway_id) {
    tryCatch({
      pathview(gene.data = kegg_id, pathway.id = pathway_id, species = "hsa", 
               out.suffix = "pathway", kegg.native = TRUE)
    }, error = function(e) {
      stop("Error visualizing pathway: ", e$message)
    })
  }
  
  # KEGG ID Sorgusu
  observeEvent(input$fetchPathway, {
    req(input$geneName)
    
    kegg_data <- fetchKEGGID(input$geneName)
    
    output$keggInfo <- renderTable({
      kegg_data
    })
    
    if (!is.null(kegg_data$KEGG_ID) && nrow(kegg_data) > 0) {
      output$pathwayPlot <- renderPlot({
        visualizePathway(kegg_id = kegg_data$KEGG_ID[1], pathway_id = "hsa04110") # ??rnek pathway ID
      })
    } else {
      output$pathwayPlot <- renderPlot({
        plot(1, type = "n", main = "No pathway visualization available.")
      })
    }
  })
  
  # Y??klenen Dosyay?? ??nizleme
  output$previewTable <- renderDataTable({
    req(input$fileInput)
    uploaded_data <- read.csv(input$fileInput$datapath)
    datatable(uploaded_data)
  })
  
  # Dinamik Yol Haritas?? G??rselle??tirme
  observeEvent(input$viewPathway, {
    req(input$pathwayID)
    
    output$livePathwayPlot <- renderPlot({
      tryCatch({
        pathview(gene.data = NULL, pathway.id = input$pathwayID, species = "hsa", kegg.native = TRUE)
      }, error = function(e) {
        plot(1, type = "n", main = paste("Error:", e$message))
      })
    })
  })
  
  # Motif Arama Fonksiyonu
  observeEvent(input$searchMotif, {
    req(input$fileInput, input$motifInput)
    
    # Sekanslar?? y??kle
    sequences <- readDNAStringSet(input$fileInput$datapath)
    
    # Motifleri aray??n
    motif_matches <- lapply(sequences, function(seq) {
      matchPattern(input$motifInput, seq)
    })
    
    # Sonu??lar?? haz??rlay??n
    result_table <- data.frame(
      Sequence_Name = names(sequences),
      Matches = sapply(motif_matches, length)
    )
    
    # Sonu??lar?? g??ster
    output$motifSearchResults <- renderDataTable({
      datatable(result_table, options = list(pageLength = 10))
    })
  })
  
  # Minimum uzunlu??a g??re sekanslar?? k??rpma
  truncateToMinLength <- function(sequences) {
    min_length <- min(nchar(as.character(sequences)))
    DNAStringSet(lapply(sequences, function(seq) substr(seq, 1, min_length)))
  }
  
  # Y??klenen dosyadan minimum uzunluklu dizileri hesapla ve k??rp
  loadTruncatedSequences <- reactive({
    req(input$fileInput)
    sequences <- readDNAStringSet(input$fileInput$datapath)
    truncateToMinLength(sequences)
  })
  
  # Frekans Analizi
  observeEvent(input$analyzeFrequency, {
    sequences <- loadTruncatedSequences()
    req(length(sequences) >= 1)
    
    freq_data <- lapply(as.character(sequences), nucleotideFrequency)
    
    # Frekans tablosu
    output$frequencyTable <- renderTable({
      freq_table <- do.call(rbind, freq_data)
      rownames(freq_table) <- names(sequences)
      freq_table
    }, rownames = TRUE)
  })
  
  # ??kili Kar????la??t??rma
  observeEvent(input$compareButton, {
    sequences <- loadTruncatedSequences()
    req(length(sequences) >= 2)
    
    # ??lk iki diziyi al
    seq1 <- as.character(sequences[[1]])
    seq2 <- as.character(sequences[[2]])
    
    # Dizileri hizala
    alignment <- pairwiseAlignment(seq1, seq2, type = "global")
    matches <- alignedPattern(alignment) == alignedSubject(alignment)
    
    comparison <- data.frame(
      Position = seq_along(matches),
      Sequence1 = unlist(strsplit(as.character(alignedPattern(alignment)), "")),
      Sequence2 = unlist(strsplit(as.character(alignedSubject(alignment)), "")),
      Match = matches
    )
    
    # Kar????la??t??rma grafi??i
    output$differencePlot <- renderPlot({
      ggplot(comparison, aes(x = Position, y = as.numeric(Match), fill = Match)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
        theme_minimal() +
        labs(title = "Sequence Comparison", x = "Position", y = "Match (1/0)")
    })
    
    # Kar????la??t??rma tablosu
    output$comparisonTable <- renderTable({
      comparison
    }, rownames = FALSE)
  })
  
  # Minimum uzunlu??a g??re sekanslar?? k??rpma
  truncateToMinLength <- function(sequences) {
    min_length <- min(nchar(as.character(sequences)))
    DNAStringSet(lapply(sequences, function(seq) substr(seq, 1, min_length)))
  }
  
  # Y??klenen dosyadan minimum uzunluklu dizileri hesapla ve k??rp
  loadTruncatedSequences <- reactive({
    req(input$fileInput)
    sequences <- readDNAStringSet(input$fileInput$datapath)
    truncateToMinLength(sequences)
  })
  
  # Is?? Haritas??
  observeEvent(input$compareHeatmap, {
    sequences <- loadTruncatedSequences()
    req(length(sequences) >= 2)
    
    n <- length(sequences)
    mat <- matrix(0, nrow = n, ncol = n, dimnames = list(names(sequences), names(sequences)))
    
    # T??m dizileri birbirleriyle kar????la??t??r
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        mat[i, j] <- sum(compareStrings(as.character(sequences[[i]]), as.character(sequences[[j]])) != "|")
      }
    }
    
    # Is?? haritas?? i??in veri ??er??evesi olu??tur
    heatmap_data <- melt(mat)
    colnames(heatmap_data) <- c("Sequence1", "Sequence2", "Differences")
    
    # Is?? Haritas?? ??izimi
    output$heatmapPlot <- renderPlot({
      ggplot(heatmap_data, aes(x = Sequence1, y = Sequence2, fill = Differences)) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "red") +
        labs(title = "Difference Heatmap", x = "Sequence 1", y = "Sequence 2", fill = "Differences") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
  })
  # Bootstrap Analizi ve Filogenetik A??a G??venilirlik De??erlerini Eklemek
  observeEvent(input$runBootstrap, {
    req(input$fileInput)
    
    # Y??klenen dizileri oku
    sequences <- readDNAStringSet(input$fileInput$datapath)
    
    # ??oklu dizi hizalamas?? (ClustalW ile)
    alignment <- msa(sequences, method = "ClustalW")
    alignment <- as.phyDat(alignment, type = "DNA")
    
    # Filogenetik A??a i??in Mesafe Matrisi Hesaplama
    dist_matrix <- dist.ml(alignment)
    tree <- nj(dist_matrix)  # Neighbor Joining a??ac?? olu??tur
    
    # Bootstrap analizi
    bs_results <- bootstrap.phyDat(alignment, FUN = function(x) nj(dist.ml(x)), bs = 100)
    
    # Bootstrap de??erlerini a??aca ekleme
    tree_with_bootstrap <- plotBS(tree, bs_results, p = 50, type = "phylogram")
    
    # Filogenetik a??ac?? ve bootstrap de??erlerini g??rselle??tirme
    output$phyloPlot <- renderPlot({
      plot(tree_with_bootstrap, main = "Phylogenetic Tree with Bootstrap Values")
    })
  })
  
  # NeighborNet A??a?? G??rselle??tirme
  observeEvent(input$runNeighborNet, {
    req(input$fileInput)
    
    # Y??klenen dizileri oku
    sequences <- readDNAStringSet(input$fileInput$datapath)
    
    # ??oklu dizi hizalamas?? (ClustalW ile)
    alignment <- msa(sequences, method = "ClustalW")
    alignment <- as.phyDat(alignment, type = "DNA")
    
    # Mesafe Matrisi Hesaplama
    dist_matrix <- dist.ml(alignment)
    
    # NeighborNet A??a?? ??emas?? olu??turma
    neighbornet_network <- neighborNet(dist_matrix)
    
    # G??rselle??tirme (etiketlerin konumunu daha d??zg??n hale getirme)
    output$neighborNetPlot <- renderPlot({
      plot(neighbornet_network, 
           main = "NeighborNet Phylogenetic Network", 
           type = "2D", 
           show.node.label = TRUE, 
           cex = 0.7, 
           font = 2,    # Etiket fontunu kal??n yap
           no.margin = TRUE)
    })
  })
  # CpG Adalar??n?? Bulan Fonksiyon
  find_CpG <- function(dna_sequence) {
    # CpG adas?? "CG" dizisini kullan??yoruz
    cpg_pattern <- "CG"
    positions <- gregexpr(cpg_pattern, dna_sequence)  # CpG adalar??n?? arar
    matches <- unlist(positions)
    
    if (length(matches) > 0) {
      result <- data.frame(Position = matches)
    } else {
      result <- data.frame(Position = NA)
    }
    
    return(result)
  }
  
  # Dosya Y??kleme ve CpG Adas?? Analizi
  observeEvent(input$analyzeCpG, {
    req(input$fileInputCpG)  # Dosyan??n y??klenmi?? olmas??n?? sa??la
    
    # Dosyay?? oku
    dna_sequence <- readDNAStringSet(input$fileInputCpG$datapath)
    dna_sequence <- as.character(dna_sequence[[1]])  # FASTA dosyas??ndaki ilk diziyi al
    
    # CpG adalar??n?? bulma
    cpg_results <- find_CpG(dna_sequence)
    
    # CpG sonu??lar??n?? tablo olarak g??stermek
    output$cpgTable <- renderTable({
      cpg_results
    })
    
    # CpG adalar?? ile ilgili basit bir g??rselle??tirme
    output$cpgPlot <- renderPlot({
      # CpG adalar??n??n bulundu??u pozisyonlar?? g??rselle??tir
      ggplot(cpg_results, aes(x = Position)) +
        geom_bar(stat = "count", fill = "blue") +
        labs(title = "CpG Islands Distribution", x = "Position", y = "Frequency") +
        theme_minimal()
    })
  })
  
  observeEvent(input$visualizeHistones, {
    req(input$histoneFile)
    
    histone_data <- read.csv(input$histoneFile$datapath)  # Histon verisini y??kle
    output$histonePlot <- renderPlot({
      plotHistoneModifications(histone_data)  # Histon modifikasyonlar??n?? g??rselle??tir
    })
  })
  observeEvent(input$runRNASeq, {
    req(input$countFile, input$conditionFile)
    
    count_data <- read.csv(input$countFile$datapath, row.names = 1)  # Say??m verilerini y??kle
    condition_data <- read.csv(input$conditionFile$datapath)  # Durum verilerini y??kle
    
    results <- analyzeRNASeq(count_data, condition_data)  # RNA-Seq analizi yap
    output$rnaSeqTable <- renderTable({
      results$table  # Sonu??lar?? tablo olarak g??ster
    })
    
    output$rnaSeqPlot <- renderPlot({
      plotRNASeq(results$table)  # RNA-Seq sonu??lar??n?? g??rselle??tir
    })
  })
  observeEvent(input$analyzeExpression, {
    req(input$geneExpFile)
    results <- analyzeGeneExpression(input$geneExpFile$datapath)  # Gen ifadesi analizi
    output$geneExpPlot <- renderPlot({ results$plot })  # Gen ifadesi grafi??ini g??ster
    output$cpmTable <- renderTable({ results$cpm_values })  # CPM de??erlerini tablo olarak g??ster
  })
  
  observeEvent(input$findDiffExpression, {
    req(input$diffExpFile)
    results <- findDifferentialExpression(input$diffExpFile$datapath)  # Farkl?? gen ifadelerini tespit et
    output$diffExpResults <- renderTable({
      results$table  # Sonu??lar?? tablo olarak g??ster
    })
  })
  
  observeEvent(input$runSimulation, {
    results <- simulateGeneExpression(input$numGenes, input$numSamples, input$noiseLevel)  # Gen ifadesi sim??lasyonu
    output$simulationPlot <- renderPlot({
      results$plot  # Sim??lasyon sonu??lar??n?? g??rselle??tir
    })
    output$simulatedData <- renderTable({
      results$simulated_data  # Sim??le edilmi?? veriyi tablo olarak g??ster
    })
  })
  GeneticNetworkModuleServer("geneticNetwork1")
  MutationSimulationModuleServer("mutationSimulation")
  ReportingModuleServer("reporting")
  StructureViewerServer("structureViewer")
  CRISPRModuleServer("crispr1")
  AllelFrequencyModuleServer("alleleFreq1")
  DifferentialExpressionServer("DE1")
  HardyWeinbergServer("hwTestModule")
  histoneModificationsServer("histone1")
  RNASeqModuleServer("rnaSeq1")
  RNASeqModuleServer("geneExpressionSimulation")
  callModule(auth_server, "auth")
  epigeneticsServer("epigenetics1")
  rnaAnalysisServer("rnaAnalysisModule")
  cancerCRISPRServer("cancerCRISPRModule")
  timeTreeServer("timeTreeModule")
  coalescentServer("coalescentModule")
  phylogeneticServer("phylogeneticModule")
  dataIntegrationServer("dataIntegrationModule")
  diseasePredictionPPIServer("diseasePredictionPPIModule")
  geneAnnotationServer("geneAnnotationModule")
  
}
# Run the application
shinyApp(ui = ui, server = server)
