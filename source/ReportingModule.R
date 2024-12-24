ReportingModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    textInput(ns("reportTitle"), "Report Title", value = "Analysis Report"),
    textInput(ns("reportAuthor"), "Author", placeholder = "Enter your name..."),
    textAreaInput(ns("reportContent"), "Report Content", rows = 6, placeholder = "Write the report details here..."),
    dateInput(ns("reportDate"), "Report Date", value = Sys.Date()),
    selectInput(ns("reportTheme"), "Select Theme", choices = c("Default", "Professional", "Modern")),
    actionButton(ns("generateReport"), "Preview Report"),
    verbatimTextOutput(ns("previewReport")),
    downloadButton(ns("downloadReportTxt"), "Download as TXT"),
    downloadButton(ns("downloadReportPdf"), "Download as PDF")
  )
}

ReportingModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    library(rmarkdown)
    
    # Reactive report content
    reportContent <- reactive({
      paste(
        "Title: ", input$reportTitle, "\n",
        "Author: ", input$reportAuthor, "\n",
        "Date: ", input$reportDate, "\n",
        "Theme: ", input$reportTheme, "\n\n",
        input$reportContent
      )
    })
    
    # Report preview
    output$previewReport <- renderText({
      reportContent()
    })
    
    # Download as TXT
    output$downloadReportTxt <- downloadHandler(
      filename = function() {
        paste0("Report_", Sys.Date(), ".txt")
      },
      content = function(file) {
        writeLines(reportContent(), file)
      }
    )
    
    # Download as PDF
    output$downloadReportPdf <- downloadHandler(
      filename = function() {
        paste0("Report_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        tempReport <- tempfile(fileext = ".Rmd")
        fileContent <- paste(
          "---",
          "title: \"", input$reportTitle, "\"",
          "author: \"", input$reportAuthor, "\"",
          "date: \"", input$reportDate, "\"",
          "output: pdf_document",
          "---\n\n",
          input$reportContent,
          sep = "\n"
        )
        writeLines(fileContent, tempReport)
        rmarkdown::render(tempReport, output_file = file)
      }
    )
  })
}