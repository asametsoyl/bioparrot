# Allel Frekans?? Mod??l?? UI
AllelFrequencyModuleUI <- function(id) {
  ns <- NS(id)  # mod??l id'sini kullanarak g??venli isimler olu??tur
  tagList(
    fileInput(ns("genotypeFile"), "Upload Genotype Data (CSV)", accept = ".csv"),
    textInput(ns("alleleColumn"), "Enter Allele Column", placeholder = "Example: Genotype"),
    actionButton(ns("analyzeButton"), "Calculate Allele Frequency"),
    tableOutput(ns("frequencyTable")),
    plotOutput(ns("frequencyPlot"))
  )
}

# Allel Frekans?? Mod??l?? Server
AllelFrequencyModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # Allel frekans?? hesaplamak i??in reaktif veri
    alleleFreq <- reactive({
      req(input$genotypeFile)
      
      # Dosya verisini y??kle
      genotypeData <- read.csv(input$genotypeFile$datapath)
      
      # Kullan??c??n??n se??ti??i allele kolonunu al
      alleleColumn <- input$alleleColumn
      
      # Kolonun mevcut olup olmad??????n?? kontrol et
      if (!alleleColumn %in% colnames(genotypeData)) {
        showModal(modalDialog(
          title = "Error",
          "The selected column was not found. Please enter a valid column name.",
          easyClose = TRUE
        ))
        return(NULL)
      }
      
      # Allel verisini ay??kla
      alleles <- unlist(strsplit(as.character(genotypeData[[alleleColumn]]), split = "[/]"))
      
      # Allellerin frekans??n?? hesapla
      freqTable <- table(alleles) / length(alleles)
      freqData <- as.data.frame(freqTable)
      colnames(freqData) <- c("Allele", "Frequency")
      
      return(freqData)
    })
    
    # Frekans tablosunu render et
    output$frequencyTable <- renderTable({
      req(alleleFreq())
      alleleFreq()
    })
    
    # Frekans grafi??ini render et
    output$frequencyPlot <- renderPlot({
      req(alleleFreq())
      freqData <- alleleFreq()
      
      if (is.null(freqData)) return(NULL)  # E??er ge??erli veri yoksa hi??bir ??ey render etme
      
      # Barplot olu??turma
      ggplot(freqData, aes(x = Allele, y = Frequency, fill = Allele)) +
        geom_bar(stat = "identity") +
        labs(title = "Allel Frekans??", x = "Allel", y = "Frekans") +
        scale_fill_brewer(palette = "Set3") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))  # X ekseni etiketlerini d??nd??r
    })
    
  })
}
