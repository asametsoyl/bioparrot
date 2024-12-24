library(shiny)

# === FONKS??YONLAR ===

# Mod??l Kullan??m Verilerini Depolama
trackModuleUsage <- function(module_name) {
  # Kullan??m verilerini bir dosyada sakla
  usage_file <- "module_usage.csv"
  
  if (!file.exists(usage_file)) {
    usage_data <- data.frame(Module = character(), Count = numeric())
  } else {
    usage_data <- read.csv(usage_file, stringsAsFactors = FALSE)
  }
  
  if (module_name %in% usage_data$Module) {
    usage_data$Count[usage_data$Module == module_name] <- usage_data$Count[usage_data$Module == module_name] + 1
  } else {
    usage_data <- rbind(usage_data, data.frame(Module = module_name, Count = 1))
  }
  
  write.csv(usage_data, usage_file, row.names = FALSE)
}

# En S??k Kullan??lan Mod??lleri Getir
getFrequentModules <- function() {
  usage_file <- "module_usage.csv"
  
  if (!file.exists(usage_file)) {
    return(data.frame(Module = character(), Count = numeric()))
  }
  
  usage_data <- read.csv(usage_file, stringsAsFactors = FALSE)
  usage_data <- usage_data[order(-usage_data$Count), ]
  head(usage_data, 5)  # ??lk 5 mod??l
}

# === UI VE SERVER ===

frequentModulesUI <- function(id) {
  ns <- NS(id)
  
  tabPanel(
    "Frequent Modules",
    verbatimTextOutput(ns("frequentModules")),
    actionButton(ns("refreshModules"), "Refresh")
  )
}

frequentModulesServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Kullan??m Verilerini Yenile
    observeEvent(input$refreshModules, {
      frequent_modules <- getFrequentModules()
      output$frequentModules <- renderText({
        if (nrow(frequent_modules) == 0) {
          "No modules have been tracked yet."
        } else {
          paste("Top Frequent Modules:\n", paste(frequent_modules$Module, ":", frequent_modules$Count, "uses"))
        }
      })
    })
  })
}
