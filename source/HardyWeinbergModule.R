library(genetics)
library(shiny)

HardyWeinbergUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Hardy-Weinberg Balance Test"),
    sidebarLayout(
      sidebarPanel(
        fileInput(ns("genotypeFile"), "Upload Genotype Data (CSV Format)", accept = c(".csv")),
        textInput(ns("genotypeColumn"), "Genotype Column Name", value = "Genotype"),
        numericInput(ns("significanceLevel"), "Significance Level (p)", value = 0.05, min = 0.001, max = 1, step = 0.001),
        actionButton(ns("performHWTest"), "Apply the Hardy-Weinberg Test", class = "btn-primary"),
        p("Verify that the genotype data in the CSV file is under the 'Genotype' heading and that there is genotype data such as 'AA', 'Aa', 'aa'.")
      ),
      mainPanel(
        verbatimTextOutput(ns("hwTestResults")),
        plotOutput(ns("genotypeDistributionPlot")),
        downloadButton(ns("downloadReport"), "Download Results")
      )
    )
  )
}

HardyWeinbergServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    hwTest <- reactive({
      req(input$genotypeFile)
      
      # CSV dosyas??n?? oku
      genotypeData <- tryCatch({
        read.csv(input$genotypeFile$datapath)
      }, error = function(e) {
        stop("An error occurred while reading the CSV file. Please check the file format: ", e$message)
      })
      
      genotypeColumn <- input$genotypeColumn
      if (!genotypeColumn %in% colnames(genotypeData)) {
        stop(paste("The selected column was not found:", genotypeColumn))
      }
      
      # Genotip verilerini kontrol et ve d??n????t??r
      genotypes <- tryCatch({
        factor(genotypeData[[genotypeColumn]])
      }, error = function(e) {
        stop("Genotype data could not be read. Please check the format:", e$message)
      })
      
      # NA kontrol??
      if (any(is.na(genotypes))) {
        stop("Missing (NA) values were found in genotype data. Please check the data.")
      }
      
      # Genotipleri uygun formata d??n????t??r
      genotypes <- tryCatch({
        genetics::genotype(genotypes)
      }, error = function(e) {
        stop("Genotypes could not be converted to the appropriate format. Check the data: ", e$message)
      })
      
      # Hardy-Weinberg Testini yap
      hw_result <- tryCatch({
        genetics::HWTest(genotypes)
      }, error = function(e) {
        stop("An error occurred while applying the Hardy-Weinberg test. Check the data: ", e$message)
      })
      
      return(list(hw_result = hw_result, genotypes = genotypes))
    })
    
    output$hwTestResults <- renderPrint({
      req(hwTest())
      hw_result <- hwTest()$hw_result
      cat("Hardy-Weinberg Balance Test Results:\n")
      print(summary(hw_result))
      if (summary(hw_result)["p.value"] < input$significanceLevel) {
        cat("\nNot: These data are NOT in Hardy-Weinberg equilibrium (p <", input$significanceLevel, ").\n")
      } else {
        cat("\nNot: These data are in Hardy-Weinberg equilibrium (p >=", input$significanceLevel, ").\n")
      }
    })
    
    output$genotypeDistributionPlot <- renderPlot({
      req(input$genotypeFile)
      genotypeData <- read.csv(input$genotypeFile$datapath)
      genotypeColumn <- input$genotypeColumn
      
      if (!genotypeColumn %in% colnames(genotypeData)) {
        stop("The selected column was not found.")
      }
      
      genotypes <- table(genotypeData[[genotypeColumn]])
      barplot(genotypes, 
              main = "Genotype Distribution", 
              col = "skyblue", 
              xlab = "Genotype", 
              ylab = "Frequency", 
              border = "darkblue")
    })
    
    output$downloadReport <- downloadHandler(
      filename = function() {
        paste("HardyWeinberg_Report_", Sys.Date(), ".txt", sep = "")
      },
      content = function(file) {
        hw_result <- hwTest()$hw_result
        genotypes <- hwTest()$genotypes
        
        report <- paste(
          "Hardy-Weinberg Equilibrium Test Results\n",
          "----------------------------------\n",
          paste(capture.output(summary(hw_result)), collapse = "\n"),
          "\n\nGenotype Frequencies:\n",
          paste(names(table(genotypes)), table(genotypes), sep = ": ", collapse = "\n"),
          sep = ""
        )
        
        writeLines(report, file)
      }
    )
  })
}
