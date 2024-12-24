library(shiny)
library(edgeR)
library(ggplot2)

# RNA-Seq Mod??l?? UI
RNASeqModuleUI <- function(id) {
  ns <- NS(id)  # Namespace tan??mlama
  tagList(
    tabsetPanel(
      tabPanel("Upload Data",
               sidebarLayout(
                 sidebarPanel(
                   fileInput(ns("countFile"), "Upload Count Data (CSV)", accept = ".csv"),
                   fileInput(ns("conditionFile"), "Upload Condition Data (CSV)", accept = ".csv"),
                   actionButton(ns("analyze"), "Analyze Data")
                 ),
                 mainPanel(
                   verbatimTextOutput(ns("summary"))
                 )
               )
      ),
      tabPanel("Plots",
               fluidRow(
                 column(6, plotOutput(ns("barPlot"))),
                 column(6, plotOutput(ns("volcanoPlot")))
               )
      ),
      tabPanel("Download Results",
               downloadButton(ns("downloadResults"), "Download Results as CSV")
      )
    )
  )
}

# RNA-Seq Mod??l?? Server
RNASeqModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    reactiveData <- reactiveValues(
      countData = NULL,
      conditionData = NULL,
      results = NULL
    )
    
    # Veri Y??kleme
    observeEvent(input$countFile, {
      req(input$countFile)
      reactiveData$countData <- read.csv(input$countFile$datapath, row.names = 1)
    })
    
    observeEvent(input$conditionFile, {
      req(input$conditionFile)
      reactiveData$conditionData <- read.csv(input$conditionFile$datapath)
    })
    
    # Analizi ??al????t??rma
    observeEvent(input$analyze, {
      req(reactiveData$countData, reactiveData$conditionData)
      reactiveData$results <- analyzeRNASeq(
        count_data = reactiveData$countData,
        condition_data = reactiveData$conditionData,
        fdr_threshold = 0.05,
        logFC_threshold = 0.25
      )
      
    })
    
    # ??zet Sonu??lar
    output$summary <- renderPrint({
      req(reactiveData$results)
      head(reactiveData$results)
    })
    
    # Bar Grafi??i
    output$barPlot <- renderPlot({
      req(reactiveData$results)
      plotRNASeq(reactiveData$results)
    })
    
    # Volkan Grafi??i
    output$volcanoPlot <- renderPlot({
      req(reactiveData$results)
      plotVolcano(reactiveData$results)
    })
    
    # Sonu??lar?? ??ndirme
    output$downloadResults <- downloadHandler(
      filename = function() {
        paste("rna_seq_results", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        req(reactiveData$results)
        write.csv(reactiveData$results, file, row.names = TRUE)
      }
    )
  })
}
