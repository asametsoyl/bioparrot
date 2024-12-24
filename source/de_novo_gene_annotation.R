# Gerekli K??t??phaneler
library(shiny)
library(Biostrings)
library(DT)
library(ggplot2)

# === FONKS??YONLAR ===

# ORF Bulma Fonksiyonu
findORFs <- function(sequence, min_length = 100) {
  start_codon <- "ATG"
  stop_codons <- c("TAA", "TAG", "TGA")
  
  starts <- unlist(gregexpr(start_codon, sequence))
  stops <- unlist(lapply(stop_codons, function(codon) gregexpr(codon, sequence)))
  stops <- sort(unlist(stops))
  
  orfs <- data.frame(Start = integer(), Stop = integer(), Length = integer())
  for (start in starts) {
    stop <- stops[stops > start & (stops - start) %% 3 == 0]
    if (length(stop) > 0) {
      orf_length <- stop[1] - start + 3
      if (orf_length >= min_length) {
        orfs <- rbind(orfs, data.frame(Start = start, Stop = stop[1], Length = orf_length))
      }
    }
  }
  
  return(orfs)
}

# Promoter B??lgesi Tahmini
findPromoterRegions <- function(sequence, promoter_motif = "TATA") {
  positions <- gregexpr(promoter_motif, sequence)[[1]]
  if (positions[1] == -1) return(data.frame(Position = integer(0)))
  return(data.frame(Position = positions))
}

# GC ????eri??i Hesaplama
calculateGCContent <- function(sequence) {
  gc_count <- sum(str_count(sequence, c("G", "C")))
  total_count <- nchar(sequence)
  return(gc_count / total_count)
}

# === UI VE SERVER ===

geneAnnotationUI <- function(id) {
  ns <- NS(id)
  
  tabsetPanel(
    tabPanel("ORF Analysis",
             fileInput(ns("orfFile"), "Upload DNA Sequence (FASTA format)"),
             numericInput(ns("minORFLength"), "Minimum ORF Length", value = 100, min = 1),
             actionButton(ns("findORFs"), "Find ORFs"),
             DTOutput(ns("orfResults")),
             plotOutput(ns("orfPlot"))
    ),
    tabPanel("Promoter Regions",
             fileInput(ns("promoterFile"), "Upload DNA Sequence (FASTA format)"),
             textInput(ns("promoterMotif"), "Promoter Motif (e.g., TATA)", value = "TATA"),
             actionButton(ns("findPromoters"), "Find Promoter Regions"),
             DTOutput(ns("promoterResults")),
             plotOutput(ns("promoterPlot"))
    ),
    tabPanel("GC Content",
             fileInput(ns("gcFile"), "Upload DNA Sequence (FASTA format)"),
             actionButton(ns("calculateGC"), "Calculate GC Content"),
             verbatimTextOutput(ns("gcContent"))
    )
  )
}

geneAnnotationServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # ORF Analizi
    observeEvent(input$findORFs, {
      req(input$orfFile)
      
      sequence_data <- tryCatch({
        readDNAStringSet(input$orfFile$datapath)
      }, error = function(e) {
        showNotification(paste("Error loading sequence data:", e$message), type = "error")
        return(NULL)
      })
      
      orf_results <- tryCatch({
        sequence <- as.character(sequence_data[[1]])
        findORFs(sequence, input$minORFLength)
      }, error = function(e) {
        showNotification(paste("Error finding ORFs:", e$message), type = "error")
        return(NULL)
      })
      
      output$orfResults <- renderDT({
        if (is.null(orf_results) || nrow(orf_results) == 0) {
          data.frame(Message = "No ORFs found.")
        } else {
          orf_results
        }
      })
      
      output$orfPlot <- renderPlot({
        if (is.null(orf_results) || nrow(orf_results) == 0) {
          ggplot() + labs(title = "No ORFs Found") + theme_void()
        } else {
          ggplot(orf_results, aes(x = Start, xend = Stop, y = 1, yend = 1)) +
            geom_segment(size = 2, color = "blue") +
            geom_text(aes(label = paste0("Length: ", Length), x = (Start + Stop) / 2, y = 1.1)) +
            labs(title = "ORFs in DNA Sequence", x = "Position", y = "") +
            theme_minimal() +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
        }
      })
    })
    
    # Promoter B??lgesi Tahmini
    observeEvent(input$findPromoters, {
      req(input$promoterFile)
      
      sequence_data <- tryCatch({
        readDNAStringSet(input$promoterFile$datapath)
      }, error = function(e) {
        showNotification(paste("Error loading sequence data:", e$message), type = "error")
        return(NULL)
      })
      
      promoter_results <- tryCatch({
        sequence <- as.character(sequence_data[[1]])
        findPromoterRegions(sequence, input$promoterMotif)
      }, error = function(e) {
        showNotification(paste("Error finding promoters:", e$message), type = "error")
        return(NULL)
      })
      
      output$promoterResults <- renderDT({
        promoter_results
      })
      
      output$promoterPlot <- renderPlot({
        if (is.null(promoter_results) || nrow(promoter_results) == 0) {
          ggplot() + labs(title = "No Promoter Regions Found") + theme_void()
        } else {
          ggplot(promoter_results, aes(x = Position, y = 1)) +
            geom_point(size = 4, color = "red") +
            labs(title = "Promoter Regions in DNA Sequence", x = "Position", y = "") +
            theme_minimal() +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
        }
      })
    })
    
    # GC ????eri??i Hesaplama
    observeEvent(input$calculateGC, {
      req(input$gcFile)
      
      sequence_data <- tryCatch({
        readDNAStringSet(input$gcFile$datapath)
      }, error = function(e) {
        showNotification(paste("Error loading sequence data:", e$message), type = "error")
        return(NULL)
      })
      
      gc_content <- tryCatch({
        sequence <- as.character(sequence_data[[1]])
        calculateGCContent(sequence)
      }, error = function(e) {
        showNotification(paste("Error calculating GC content:", e$message), type = "error")
        return(NULL)
      })
      
      output$gcContent <- renderText({
        paste("GC Content:", round(gc_content * 100, 2), "%")
      })
    })
  })
}
