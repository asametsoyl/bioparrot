# Server Fonksiyonu
GeneticNetworkModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    observeEvent(input$generateNetwork, {
      req(input$networkFile)
      
      # Veri dosyas??n?? y??kle
      networkData <- tryCatch({
        read.csv(input$networkFile$datapath)
      }, error = function(e) {
        showModal(modalDialog(
          title = "Error",
          "Failed to read the network data. Please ensure the file is in CSV format.",
          easyClose = TRUE
        ))
        return(NULL)
      })
      
      # Verinin uygun olup olmad??????n?? kontrol et
      if (is.null(networkData)) return()
      
      # A?? verisini kontrol et
      if (!all(c("Source", "Target") %in% colnames(networkData))) {
        showModal(modalDialog(
          title = "Error",
          "The CSV file must contain 'Source' and 'Target' columns.",
          easyClose = TRUE
        ))
        return()
      }
      
      # igraph paketini y??kle ve a?? olu??tur
      library(igraph)
      geneticNetwork <- graph_from_data_frame(networkData, directed = TRUE)
      
      # A????n istatistiklerini hesapla
      networkStats <- list(
        num_nodes = vcount(geneticNetwork),
        num_edges = ecount(geneticNetwork),
        density = graph.density(geneticNetwork),  # Do??ru fonksiyon ad?? kullan??ld??
        avg_degree = mean(degree(geneticNetwork))
      )
      
      # A?? istatistiklerini g??ster
      output$networkStats <- renderPrint({
        paste(
          "Number of Nodes: ", networkStats$num_nodes, "\n",
          "Number of Edges: ", networkStats$num_edges, "\n",
          "Network Density: ", round(networkStats$density, 4), "\n",
          "Average Degree: ", round(networkStats$avg_degree, 2)
        )
      })
      
      # A?? g??rselle??tirmesini olu??tur
      output$geneticNetworkPlot <- renderPlot({
        plot(geneticNetwork, vertex.size = 15, vertex.label.cex = 0.8, 
             edge.arrow.size = 0.5, main = "Genetic Network",
             vertex.color = "lightblue", edge.color = "gray")
      })
    })
  })
}

