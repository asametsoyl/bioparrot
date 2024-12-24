auth_server <- function(input, output, session) {
  
  # Reactive de??er olarak giri?? durumunu tutuyoruz
  is_logged_in <- reactiveVal(FALSE)
  
  # Kullan??c?? giri?? ve kay??t formunun dinamik olarak g??sterilmesi
  output$auth_ui <- renderUI({
    if (is_logged_in()) {
      return(NULL)  # Giri?? yap??lm????sa, giri?? ekran??n?? gizle
    }
    
    fluidRow(
      column(4, offset = 4,
             wellPanel(
               textInput("username", "Kullan??c?? Ad??"),
               passwordInput("password", "??ifre"),
               actionButton("login_btn", "Giri?? Yap"),
               actionButton("register_btn", "??ye Ol")
             )
      )
    )
  })
  
  # Kullan??c?? veritaban?? i??lemleri i??in fonksiyonlar
  register_user <- function(username, password) {
    con <- dbConnect(RSQLite::SQLite(), "users.db")
    existing_user <- dbGetQuery(con, paste0("SELECT * FROM users WHERE username = '", username, "'"))
    
    if (nrow(existing_user) > 0) {
      dbDisconnect(con)
      return("Kullan??c?? ad?? zaten mevcut.")
    }
    
    dbExecute(con, "INSERT INTO users (username, password) VALUES (?, ?)", params = list(username, password))
    dbDisconnect(con)
    return("Kay??t ba??ar??l??!")
  }
  
  login_user <- function(username, password) {
    con <- dbConnect(RSQLite::SQLite(), "users.db")
    query <- dbGetQuery(con, paste0("SELECT * FROM users WHERE username = '", username, "' AND password = '", password, "'"))
    dbDisconnect(con)
    
    if (nrow(query) > 0) {
      return(TRUE)  # Ba??ar??l?? giri??
    } else {
      return(FALSE)  # Ba??ar??s??z giri??
    }
  }
  
  # Giri?? i??lemi
  observeEvent(input$login_btn, {
    username <- input$username
    password <- input$password
    
    if (login_user(username, password)) {
      # Ba??ar??l?? giri??
      is_logged_in(TRUE)  # Giri?? ba??ar??l??ysa giri?? durumunu g??ncelle
    } else {
      # Giri?? ba??ar??s??z
      showModal(modalDialog(
        title = "Hata",
        "Yanl???? kullan??c?? ad?? veya ??ifre.",
        easyClose = TRUE
      ))
    }
  })
  
  # Kay??t i??lemi
  observeEvent(input$register_btn, {
    username <- input$username
    password <- input$password
    
    result <- register_user(username, password)
    showModal(modalDialog(
      title = "Sonu??",
      result,
      easyClose = TRUE
    ))
  })
  
  # ????k???? i??lemi
  observeEvent(input$logout_btn, {
    is_logged_in(FALSE)  # ????k???? yap??ld??????nda giri?? durumu s??f??rlan??r
  })
  
  # Giri?? durumu geri d??nd??r??l??r
  return(is_logged_in)
}
