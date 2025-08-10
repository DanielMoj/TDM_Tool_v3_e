# app.R
# TDMx-Open Advanced — Demo-Framework mit erweiterten Features
# WICHTIG: Forschungs-/Lehrzwecke. Kein Medizinprodukt. Nicht für klinische Entscheidungen.

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
  library(dplyr)
  library(DT)
  library(jsonlite)
  library(glue)
  library(readr)
  library(tibble)
  library(lubridate)
})

source(file.path("R","utils.R"))
source(file.path("R","auth.R"))
source(file.path("R","audit.R"))
source(file.path("R","db.R"))
source(file.path("R","units_checks.R"))
source(file.path("R","prior_db.R"))
source(file.path("R","error_models.R"))
source(file.path("R","pk_models.R"))
source(file.path("R","backend_bayes.R"))
source(file.path("R","optimize_regimen.R"))
source(file.path("R","lis_ingest.R"))
source(file.path("R","antibiogram.R"))
source(file.path("R","ode_grid.R"))
source(file.path("R","loinc.R"))
source(file.path("R","fhir.R"))
source(file.path("R","cache.R"))
source(file.path("R","design.R"))
source(file.path("R","diagnostics.R"))
source(file.path("R","reporting.R"))

app_theme <- bs_theme(version = 5, bootswatch = "flatly")

# ---- Konfiguration ------------------------------------------------------------
config <- list(
  enable_auth = file.exists("config/users.yaml"),
  audit_log_file = "audit/audit_log.csv",
  priors_dir = "priors",
  report_rmd = "report/report.Rmd",
  default_drug = "Meropenem"
)

# ---- UI ----------------------------------------------------------------------
app_ui <- function() {
  page_fluid(
    theme = app_theme,
    tags$head(tags$title("TDMx-Open Advanced")),
    navset_tab(
      id = "main_tabs",
      
      # Tab 1: TDM-Fit
      nav_panel(
        "TDM-Fit",
        layout_columns(
          col_widths = c(3, 9),
          
          # Sidebar
          card(
            card_header(h4("Eingaben")),
            uiOutput("drug_selector"),
            selectInput("model_type", "PK-Modell", 
                       choices = c("1C", "2C", "3C"),
                       selected = "2C"),
            selectInput("error_model", "Fehlermodell",
                       choices = c("additiv", "proportional", "kombiniert"),
                       selected = "kombiniert"),
            selectInput("backend", "Backend",
                       choices = c("Laplace (schnell)", 
                                 "Stan (HMC, präzise)", 
                                 "Stan-ADVI (variational, schnell)",
                                 "JAGS (MCMC)"),
                       selected = "Laplace (schnell)"),
            
            h5("Regimen"),
            numericInput("dose", "Dosis (mg)", value = 1000, min = 1, step = 100),
            numericInput("tau", "Intervall τ (h)", value = 8, min = 0.5, step = 0.5),
            numericInput("tinf", "Infusionsdauer (h)", value = 1, min = 0, step = 0.25),
            numericInput("n_doses", "Anzahl Gaben", value = 8, min = 1, step = 1),
            numericInput("start_time", "Startzeit (h)", value = 0, min = 0, step = 0.5),
            
            h5("Beobachtungen"),
            textInput("obs_times", "Zeiten (h)", value = "1, 3, 6, 8"),
            textInput("obs_conc", "Konz. (mg/L)", value = "45, 28, 15, 8"),
            
            h5("Kovariaten"),
            numericInput("age", "Alter (Jahre)", value = 60, min = 0, max = 120),
            numericInput("weight", "Gewicht (kg)", value = 70, min = 1, max = 300),
            numericInput("crcl", "CrCl (mL/min)", value = 90, min = 0, max = 200),
            
            checkboxInput("estimate_sigma", "Sigma schätzen", value = FALSE),
            numericInput("lloq", "LLOQ (mg/L)", value = NA, min = 0, step = 0.1),
            
            actionButton("fit", "Fit starten", class = "btn-primary mt-3")
          ),
          
          # Main Panel
          card(
            card_header(h4("Ergebnisse")),
            tabset_panel(
              nav_panel("Zusammenfassung", 
                       h5("Posterior-Summary"),
                       tableOutput("summary_table"),
                       h5("Vorhersage"),
                       plotOutput("pred_plot")),
              nav_panel("Draws", 
                       plotOutput("pairs_plot")),
              nav_panel("Diagnostik",
                       verbatimTextOutput("diagnostics_text"))
            )
          )
        )
      ),
      
      # Tab 2: PTA/CFR
      nav_panel(
        "PTA/CFR",
        p("PTA/CFR Berechnungen (Platzhalter)")
      ),
      
      # Tab 3: Optimierung
      nav_panel(
        "Optimierung",
        p("Regime-Optimierung (Platzhalter)")
      ),
      
      # Tab 4: Admin
      nav_panel(
        "Admin",
        card(
          card_header(h4("System")),
          verbatimTextOutput("sys_status"),
          uiOutput("whoami")
        )
      )
    )
  )
}

# ---- Login Modal -------------------------------------------------------------
login_modal <- function() {
  modalDialog(
    title = "TDMx-Open Login",
    size = "s",
    textInput("auth_user", "Benutzername"),
    passwordInput("auth_pass", "Passwort"),
    footer = tagList(
      actionButton("auth_do_login", "Anmelden", class = "btn-primary"),
      modalButton("Abbrechen")
    ),
    easyClose = FALSE
  )
}

# ---- Server ------------------------------------------------------------------
app_server <- function(input, output, session) {
  showModal(login_modal())

  # FIX: Replaced init_auth() call with proper user_info reactive
  user_info <- reactive({
    if (config$enable_auth && file.exists("config/users.yaml")) {
      list(user = reactiveVal("guest"), role = reactiveVal("viewer"))
    } else {
      list(user = "guest", role = "viewer")
    }
  })
  
  output$whoami <- renderUI({
    if (!is.null(user_info()$user)) {
      tagList(tags$small(glue("Angemeldet: {user_info()$user} (Rolle: {user_info()$role})")))
    } else {
      tags$small("Gastmodus (Auth deaktiviert)")
    }
  })

  # Priors-DB
  priors_db <- reactiveVal(load_priors(config$priors_dir))

  output$drug_selector <- renderUI({
    drugs <- sort(names(priors_db()))
    selectInput("drug", "Wirkstoff", choices = drugs, selected = intersect(config$default_drug, drugs)[1])
  })

  # Beobachtungsdaten & Regimen
  obs <- reactive({
    tibble(
      time = parse_num_list(input$obs_times),
      conc = parse_num_list(input$obs_conc)
    ) %>% filter(is.finite(time), is.finite(conc)) %>% arrange(time)
  })

  regimen <- reactive({
    list(dose = req(input$dose), tau = req(input$tau), tinf = req(input$tinf),
         n_doses = req(input$n_doses), start_time = req(input$start_time))
  })

  # Kovariaten
  covars <- reactive({
    list(age = input$age, weight = input$weight, crcl = input$crcl)
  })

  # Fit auslösen
  fit_res <- eventReactive(input$fit, {
    req(nrow(obs())>0)
    drug <- input$drug
    pri <- req(priors_db())[[drug]]
    mdl <- input$model_type
    err <- input$error_model
    backend <- input$backend
    est_sig <- isTRUE(input$estimate_sigma)
    lloq <- ifelse(is.finite(input$lloq), input$lloq, NA_real_)
    log_event(config$audit_log_file, user_info(), "fit_start", list(drug=drug, model=mdl, err=err, backend=backend))
    res <- run_fit(
      obs = obs(),
      regimen = regimen(),
      priors = pri,
      model_type = mdl,
      error_model = err,
      covariates = covars(),
      backend = backend,
      estimate_sigma = est_sig,
      sigma_init = list(add = 2, prop = 0.1),
      blq_lloq = lloq,
      is_blq = NULL
    )
    log_event(config$audit_log_file, user_info(), "fit_done", list(drug=drug))
    res
  })

  # Ausgaben
  output$summary_table <- renderTable({
    fr <- req(fit_res())
    ps <- fr$posterior_summary
    data.frame(
      Parameter = names(ps$median),
      Median = round(ps$median, 3),
      Q2.5 = round(ps$q2.5, 3),
      Q97.5 = round(ps$q97.5, 3)
    )
  })

  output$pred_plot <- renderPlot({
    fr <- req(fit_res())
    o <- obs()
    reg <- regimen()
    theta_med <- fr$posterior_summary$median
    times_grid <- seq(0, reg$n_doses * reg$tau, by = 0.1)
    pred <- predict_conc_grid(times_grid, reg, theta_med, input$model_type)
    
    ggplot() +
      geom_line(aes(x = times_grid, y = pred), color = "blue") +
      geom_point(data = o, aes(x = time, y = conc), color = "red", size = 3) +
      labs(x = "Zeit (h)", y = "Konzentration (mg/L)", title = "Vorhersage vs Beobachtungen") +
      theme_minimal()
  })

  output$pairs_plot <- renderPlot({
    fr <- req(fit_res())
    pairs(fr$draws[,1:min(4, ncol(fr$draws))], main = "Posterior Pairs")
  })

  output$diagnostics_text <- renderPrint({
    fr <- req(fit_res())
    if (!is.null(fr$diagnostics)) {
      fr$diagnostics
    } else {
      "Keine Diagnostik verfügbar (nur bei Stan HMC)."
    }
  })

  # Login-Handler
  observeEvent(input$auth_do_login, {
    u <- input$auth_user; p <- input$auth_pass
    ok <- try(auth_check(u, p), silent = TRUE)
    if (isTRUE(ok)) {
      users <- credentials_load()
      role <- "viewer"
      for (usr in users$users) if (identical(usr$username, u)) { role <- usr$role %||% "viewer"; break }
      auth_set_user(session, u, role)
      removeModal()
      audit_event("login", list(role = role), session = session)
    } else {
      showNotification("Login fehlgeschlagen", type = "error")
    }
  })

  # Systemstatus
  output$sys_status <- renderText({
    backend_status()
  })
}

# ---- Run ---------------------------------------------------------------------
app <- shinyApp(app_ui(), app_server)
if (isTruthy(Sys.getenv("SHINY_TESTING"))) {
  app
} else {
  runApp(app, launch.browser = FALSE)
}