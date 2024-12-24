MutationSimulationModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    numericInput(ns("numGenes"), "Number of Genes to Simulate:", value = 100, min = 10),
    numericInput(ns("mutationRate"), "Mutation Rate (%):", value = 5, min = 0, max = 100),
    selectInput(ns("mutationType"), "Select Mutation Type:", 
                choices = c("All", "Insertion", "Deletion", "Substitution"),
                selected = "All"),
    actionButton(ns("runSimulation"), "Start Simulation"),
    downloadButton(ns("downloadData"), "Download Simulation Data"),
    plotOutput(ns("mutationSimulationPlot")),
    tableOutput(ns("mutationSimulationTable"))
  )
}

MutationSimulationModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    simulationData <- reactiveVal(NULL)  # To store simulation data
    
    observeEvent(input$runSimulation, {
      # Simulasyon verileri olusturma
      numGenes <- input$numGenes
      mutationRate <- input$mutationRate / 100
      
      set.seed(123)  # Reproducibility
      genes <- paste0("Gene", seq_len(numGenes))
      mutated <- sample(c(TRUE, FALSE), numGenes, replace = TRUE, prob = c(mutationRate, 1 - mutationRate))
      mutationType <- sample(c("Insertion", "Deletion", "Substitution"), numGenes, replace = TRUE)
      
      data <- data.frame(Gene = genes, Mutated = mutated, MutationType = ifelse(mutated, mutationType, "None"))
      
      # Store the data
      simulationData(data)
      
      # Gorsellestirme
      output$mutationSimulationPlot <- renderPlot({
        barplot(table(data$Mutated), 
                main = "Mutation Simulation", 
                col = c("blue", "red"), 
                names.arg = c("Normal", "Mutated"),
                ylim = c(0, max(table(data$Mutated)) * 1.2))
        legend("topright", legend = c("Normal", "Mutated"), fill = c("blue", "red"))
      })
      
      # Tablo
      output$mutationSimulationTable <- renderTable({
        if (input$mutationType == "All") {
          data
        } else {
          data[data$MutationType == input$mutationType, ]
        }
      })
    })
    
    # Data Download
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("mutation_simulation_data", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(simulationData(), file, row.names = FALSE)
      }
    )
  })
}
