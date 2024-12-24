# CRISPR Mod??l??n??n UI K??sm??
CRISPRModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    textInput(ns("manualInput"), "Manual DNA Sequence Input", placeholder = "Enter DNA sequence here"),
    fileInput(ns("fileInput"), "Upload DNA Sequence (FASTA format)", accept = c(".fasta")),  # Dosya y??kleme
    textInput(ns("pamInput"), "PAM Sequence", placeholder = "Default: NGG"),
    actionButton(ns("analyzeButton"), "Find CRISPR Targets"),
    DTOutput(ns("targetTable")),
    plotOutput(ns("targetPlot"))
  )
}
# CRISPR Mod??l??n??n Server K??sm??
CRISPRModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # CRISPR hedef b??lgelerini bulma fonksiyonu
    findCRISPRTargets <- function(dna_sequence, pam = "NGG") {
      pam_regex <- gsub("N", "[ATGC]", pam)  # PAM dizisini regex format??na ??evir
      pam_positions <- gregexpr(pam_regex, dna_sequence, perl = TRUE)[[1]]
      
      if (pam_positions[1] == -1) {
        return(data.frame(Target = character(), Start = integer(), End = integer()))
      }
      
      # Hedef b??lgeleri belirleme
      targets <- data.frame(Target = character(), Start = integer(), End = integer())
      for (pos in pam_positions) {
        start_pos <- max(1, pos - 20)  # PAM ??ncesindeki 20 bazl??k hedef b??lge
        end_pos <- pos + attr(pam_positions, "match.length") - 1
        target_sequence <- substr(dna_sequence, start_pos, end_pos)
        
        targets <- rbind(targets, data.frame(
          Target = target_sequence,
          Start = start_pos,
          End = end_pos
        ))
      }
      return(targets)
    }
    
    # Kullan??c?? CRISPR analizini ba??latt??????nda
    observeEvent(input$analyzeButton, {
      req(input$manualInput)  # Manuel DNA dizisinin girilmi?? olmas?? gerekti??ini kontrol et
      req(input$fileInput)  # Veya dosya y??klenmi?? olmas?? gerekti??ini kontrol et
      
      # Kullan??c??dan DNA dizisini al
      dna_sequence <- toupper(input$manualInput)
      
      # E??er dosya y??klenmi??se, y??klenen dosyay?? i??le
      if (!is.null(input$fileInput)) {
        file_path <- input$fileInput$datapath
        dna_sequence <- paste(readLines(file_path), collapse = "")
      }
      
      pam_sequence <- ifelse(nzchar(input$pamInput), toupper(input$pamInput), "NGG")
      
      # CRISPR hedef b??lgelerini bul
      targets <- findCRISPRTargets(dna_sequence, pam_sequence)
      
      if (nrow(targets) == 0) {
        showModal(modalDialog(
          title = "No Targets Found",
          "No CRISPR targets found for the given DNA sequence and PAM sequence.",
          easyClose = TRUE
        ))
        return()
      }
      
      # Hedef b??lgeleri tablo olarak g??ster
      output$targetTable <- renderDT({
        datatable(targets, options = list(pageLength = 5), rownames = FALSE)
      })
      
      # Hedef b??lgeleri g??rselle??tir
      output$targetPlot <- renderPlot({
        ggplot(targets, aes(x = Start, xend = End, y = 1, yend = 1)) +
          geom_segment(size = 3, color = "blue", lineend = "round") +
          geom_point(aes(x = Start, y = 1), color = "green", size = 3) +
          geom_point(aes(x = End, y = 1), color = "red", size = 3) +
          labs(title = "CRISPR Target Regions", x = "Position", y = "") +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),  # X eksenindeki metinleri d??nd??r
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "none"  # Legend'?? kald??r
          ) +
          geom_text(aes(x = (Start + End) / 2, y = 1.2, label = Target), angle = 45, hjust = 1)
      })
    })
  })
}
