# Server Fonksiyonu
histoneModificationsServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    data <- reactive({
      req(input$histoneFile)
      histone_data <- read.csv(input$histoneFile$datapath)
      
      # Eksik verileri doldurma
      histone_data$Modification[is.na(histone_data$Modification)] <- "Unknown"
      histone_data$Intensity[is.na(histone_data$Intensity)] <- 0
      histone_data$Modification <- factor(histone_data$Modification, 
                                          levels = c("Acetylation", "Methylation", 
                                                     "Phosphorylation", "Ubiquitination", "Unknown"))
      histone_data
    })
    
    filteredData <- reactive({
      req(data())
      if (input$modification_filter == "All") {
        data()
      } else {
        subset(data(), Modification == input$modification_filter)
      }
    })
    
    observeEvent(input$visualizeHistones, {
      output$histonePlot <- renderPlotly({
        req(filteredData())
        
        p <- ggplot(filteredData(), aes(x = Start, y = Modification, color = Modification)) +
          geom_point(size = 3, alpha = 0.8) +
          scale_color_manual(values = c("Acetylation" = "blue", 
                                        "Methylation" = "green", 
                                        "Phosphorylation" = "purple", 
                                        "Ubiquitination" = "orange", 
                                        "Unknown" = "grey")) +
          labs(title = "Histone Modifications Start Points",
               x = "Genomic Position", y = "Modification Type") +
          theme_minimal() +
          theme(
            plot.title = element_text(face = "bold", size = 16),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right"
          )
        
        ggplotly(p) %>% layout(dragmode = "zoom")
      })
    })
  })
}
