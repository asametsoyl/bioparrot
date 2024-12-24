library(shiny)

# === CSS ===
darkModeCSS <- "
body {
  background-color: #121212 !important;
  color: #ffffff !important;
}
.navbar, .tab-content, .btn, .sidebar {
  background-color: #1f1f1f !important;
  border-color: #333333 !important;
  color: #ffffff !important;
}
.btn:hover {
  background-color: #333333 !important;
}
"

lightModeCSS <- "
body {
  background-color: #ffffff !important;
  color: #000000 !important;
}
.navbar, .tab-content, .btn, .sidebar {
  background-color: #f8f9fa !important;
  border-color: #dddddd !important;
  color: #000000 !important;
}
.btn:hover {
  background-color: #e9ecef !important;
}
"

# === UI VE SERVER ===

darkModeUI <- function(id) {
  ns <- NS(id)
  
  checkboxInput(ns("darkModeToggle"), "Enable Dark Mode", value = FALSE)
}

darkModeServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    observe({
      if (input$darkModeToggle) {
        tags$style(HTML(darkModeCSS)) %>% insertUI(selector = "head")
      } else {
        tags$style(HTML(lightModeCSS)) %>% insertUI(selector = "head")
      }
    })
  })
}

