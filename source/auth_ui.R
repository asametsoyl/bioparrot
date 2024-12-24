auth_ui <- function(id) {
  ns <- NS(id)  # Mod??l i??in namespace
  
  fluidPage(
    titlePanel("Kullan??c?? Giri?? ve Kay??t Sistemi"),
    
    # Giri?? Paneli
    uiOutput(ns("auth_ui"))
  )
}
