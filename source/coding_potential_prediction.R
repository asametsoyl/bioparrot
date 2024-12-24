# Gerekli K??t??phaneler
library(shiny)
library(randomForest)
library(Biostrings)
library(DT)
library(ggplot2)

# === FONKS??YONLAR ===

# Kodlama Potansiyeli ??zellik ????kar??m??
extractCodingPotentialFeatures <- function(sequences) {
  features <- data.frame(
    GC_Content = sapply(sequences, function(seq) {
      gc <- sum(letterFrequency(DNAString(seq), letters = c("G", "C")))
      total <- nchar(seq)
      gc / total
    }),
    Length = sapply(sequences, nchar)
  )
  
  # Kodlama olas??l?????? (??rnek bir hesaplama: ATG ile ba??lama)
  features$Start_Codon = sapply(sequences, function(seq) {
    substr(seq, 1, 3) == "ATG"
  })
  
  return(features)
}

# Model E??itimi
trainCodingPotentialModel <- function(features, labels) {
  model <- randomForest(as.factor(labels) ~ ., data = features, ntree = 100)
  return(model)
}

# Kodlama Potansiyeli Tahmini
predictCodingPotential <- function(model, features) {
  predictions <- predict(model, newdata = features)
  return(predictions)
}

# === UI VE SERVER ===

codingPotentialUI <- function(id) {
  ns <- NS(id)
  
  tabsetPanel(
    tabPanel("Train Model",
             fileInput(ns("trainFile"), "Upload Training Data (FASTA format)"),
             textInput(ns("labels"), "Enter Labels (comma-separated)", value = "Coding,Non-coding"),
             actionButton(ns("trainModel"), "Train Model"),
             verbatimTextOutput(ns("modelStatus"))
    ),
    tabPanel("Predict Coding Potential",
             fileInput(ns("testFile"), "Upload Test Data (FASTA format)"),
             actionButton(ns("predict"), "Predict Coding Potential"),
             DTOutput(ns("predictionResults")),
             plotOutput(ns("predictionPlot"))
    )
  )
}

codingPotentialServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # E??itim Modeli
    trained_model <- reactiveVal(NULL)
    
    observeEvent(input$trainModel, {
      req(input$trainFile, input$labels)
      
      # E??itim Verisini Y??kleme
      train_sequences <- tryCatch({
        readDNAStringSet(input$trainFile$datapath)
      }, error = function(e) {
        showNotification(paste("Error loading training data:", e$message), type = "error")
        return(NULL)
      })
      
      # ??zellik ????kartma
      features <- tryCatch({
        extractCodingPotentialFeatures(train_sequences)
      }, error = function(e) {
        showNotification(paste("Error extracting features:", e$message), type = "error")
        return(NULL)
      })
      
      # Etiketleri Al
      labels <- strsplit(input$labels, ",")[[1]]
      if (length(labels) != nrow(features)) {
        showNotification("Number of labels must match the number of sequences.", type = "error")
        return(NULL)
      }
      
      # Modeli E??itme
      model <- tryCatch({
        trainCodingPotentialModel(features, labels)
      }, error = function(e) {
        showNotification(paste("Error training model:", e$message), type = "error")
        return(NULL)
      })
      
      trained_model(model)
      output$modelStatus <- renderText({
        if (is.null(model)) {
          "Failed to train model."
        } else {
          "Model trained successfully."
        }
      })
    })
    
    # Kodlama Potansiyeli Tahmini
    observeEvent(input$predict, {
      req(input$testFile, trained_model())
      
      # Test Verisini Y??kleme
      test_sequences <- tryCatch({
        readDNAStringSet(input$testFile$datapath)
      }, error = function(e) {
        showNotification(paste("Error loading test data:", e$message), type = "error")
        return(NULL)
      })
      
      # ??zellik ????kartma
      test_features <- tryCatch({
        extractCodingPotentialFeatures(test_sequences)
      }, error = function(e) {
        showNotification(paste("Error extracting features:", e$message), type = "error")
        return(NULL)
      })
      
      # Tahmin Yapma
      predictions <- tryCatch({
        predictCodingPotential(trained_model(), test_features)
      }, error = function(e) {
        showNotification(paste("Error during prediction:", e$message), type = "error")
        return(NULL)
      })
      
      # Sonu??lar?? Tablo ve G??rselle??tirme Olarak G??ster
      output$predictionResults <- renderDT({
        data.frame(Sequence = names(test_sequences), Prediction = predictions)
      })
      
      output$predictionPlot <- renderPlot({
        prediction_data <- data.frame(
          Sequence = names(test_sequences),
          Prediction = predictions
        )
        ggplot(prediction_data, aes(x = Prediction, fill = Prediction)) +
          geom_bar() +
          labs(title = "Coding Potential Predictions", x = "Prediction", y = "Count") +
          theme_minimal()
      })
    })
  })
}
