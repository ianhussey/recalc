# app.R
# Shiny app for recalc::independent_t_test_summary()
# - Shows res$reproduced, plot_multiverse_p(res), plot_multiverse_d(res)
# - Provides a downloadable HTML report via report.Rmd

suppressPackageStartupMessages({
  library(shiny)
  library(recalc)
  library(rmarkdown)
})

ui <- fluidPage(
  titlePanel("recalc: Bounds of an independent t-test's p-value and Cohen's d from M/SD/N"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Summary statistics (required)"),
      numericInput("m1", "Mean group 1", value = NA, step = 0.01),
      numericInput("m2", "Mean group 2", value = NA, step = 0.01),
      numericInput("sd1", "SD group 1", value = NA, step = 0.01),
      numericInput("sd2", "SD group 2", value = NA, step = 0.01),
      numericInput("n1", "N group 1", value = NA, min = 2, step = 1),
      numericInput("n2", "N group 2", value = NA, min = 2, step = 1),
      numericInput("m_digits", "Number of digits means reported to", value = NA, min = 0, step = 1),
      numericInput("sd_digits", "Number of digits SDs reported to", value = NA, min = 0, step = 1),
      
      hr(),
      h4("Reported p-value (optional)"),
      selectInput(
        "p_operator",
        "p operator",
        choices = c("less_than", "more_than", "equals"),
        selected = "less_than"
      ),
      numericInput("p", "Reported p", value = NA, min = 0, max = 1, step = 0.01),
      numericInput("p_digits", "Digits used for p (p_digits)", value = 3, min = 0, step = 1),
      selectInput(
        "p_methods",
        "P-value method",
        choices = c("Both Student's t and Welch's t" = "NULL", 
                    "Student's t" = "student_t", 
                    "Welch's t" = "welch_t"),
        selected = "NULL"
      ),
      numericInput("alpha", "alpha", value = 0.05, min = 0, max = 1, step = 0.01),
      
      hr(),
      h4("Reported Cohen's d and 95% Confidence Intervals (optional)"),
      numericInput("d", "Reported SMD (Cohen's d or Hedges' g)", value = NA, step = 0.01),
      numericInput("d_ci_lower", "Reported d CI lower", value = NA, step = 0.01),
      numericInput("d_ci_upper", "Reported d CI upper", value = NA, step = 0.01),
      
      selectInput(
        "direction",
        "Direction of effect",
        choices = c("m1_minus_m2", "m2_minus_m1", "both"),
        selected = "m1_minus_m2"
      ),
      
      hr(),
      h4("Rounding and method options"),
      
      selectInput(
        "input_rounding",
        "Method of rounding the reported summary statistics",
        choices = c("Unknown" = "NULL",
                    "Rounded (up, down, or bankers)" = "rounded",
                    "Truncated" = "truncated"),
        selected = "NULL"
      ),
      
      selectInput(
        "output_rounding",
        "Output rounding method",
        choices = c("All methods" = "NULL",
                    "Round half up" = "half_up",
                    "Round half down" = "half_down",
                    "Bankers rounding" = "bankers",
                    "Truncation" = "trunc"),
        selected = "NULL"
      ),
      
      numericInput("d_digits", "Number of digits Cohens' d reported to", value = 2, min = 0, step = 1),
      
      selectInput(
        "hedges_correction",
        "Apply Hedges' correction",
        choices = c("Both methods" = "NULL",
                    "Apply"     = "TRUE",
                    "Do not apply" = "FALSE"),
        selected = "NULL"
      ),
      
      selectInput(
        "ci_methods",
        "CI method",
        choices = c("All methods" = "NULL",
                    "Wald" = "wald_t",
                    "Welch" = "welch_t",
                    "Non-central t" = "nct"),
        selected = "NULL"
      ),
      
      checkboxInput(
        "include_se_sd_confusion",
        "Treat SDs as if they may have been SEs confused for SDs",
        value = FALSE
      ),
      
      hr(),
      actionButton("run", "Run analysis"),
      hr(),
      downloadButton("downloadReport", "Download HTML report")
    ),
    
    mainPanel(
      tableOutput("tbl_reproduced"),
      br(),
      plotOutput("plot_p", height = "400px"),
      br(),
      plotOutput("plot_d", height = "400px")
    )
  )
)

server <- function(input, output, session) {
  
  # Helper: map selectInput strings for optional arguments to real values
  get_input_rounding <- reactive({
    if (identical(input$input_rounding, "NULL")) NULL else input$input_rounding
  })
  
  get_output_rounding <- reactive({
    if (identical(input$output_rounding, "NULL")) NULL else input$output_rounding
  })
  
  get_p_methods <- reactive({
    if (identical(input$p_methods, "NULL")) NULL else input$p_methods
  })
  
  get_hedges_correction <- reactive({
    switch(
      input$hedges_correction,
      "NULL"  = NULL,
      "TRUE"  = TRUE,
      "FALSE" = FALSE
    )
  })
  
  get_ci_methods <- reactive({
    if (identical(input$ci_methods, "NULL")) NULL else input$ci_methods
  })
  
  # Only runs when the button is clicked
  recalc_res <- eventReactive(
    input$run,
    {
      # Validate that all necessary numeric inputs are filled in
      req(
        !is.na(input$m1),
        !is.na(input$m2),
        !is.na(input$sd1),
        !is.na(input$sd2),
        !is.na(input$n1),
        !is.na(input$n2)
      )
      
      recalc::independent_t_test_summary(
        m1 = input$m1,
        m2 = input$m2,
        sd1 = input$sd1,
        sd2 = input$sd2,
        n1 = input$n1,
        n2 = input$n2,
        m_digits = input$m_digits,
        sd_digits = input$sd_digits,
        
        p_operator = input$p_operator,
        p = input$p,
        p_digits = input$p_digits,
        p_methods = get_p_methods(),
        alpha = input$alpha,
        
        d = input$d,
        d_ci_lower = input$d_ci_lower,
        d_ci_upper = input$d_ci_upper,
        direction = input$direction,
        
        input_rounding = get_input_rounding(),
        output_rounding = get_output_rounding(),
        d_digits = input$d_digits,
        hedges_correction = get_hedges_correction(),
        ci_methods = get_ci_methods(),
        include_se_sd_confusion = isTRUE(input$include_se_sd_confusion)
      )
    },
    ignoreNULL = TRUE   # before first click, recalc_res() returns NULL
  )
  
  # Outputs depend ONLY on recalc_res()
  output$tbl_reproduced <- renderTable({
    res <- recalc_res()
    req(res)               # ensures the button has been clicked successfully
    res$reproduced |>
      select(p_operator,
             p_reported = p,
             p_min_rounded = min_p_rounded,
             p_max_rounded = max_p_rounded,
             p_in_bounds_given_operator = p_inbounds,
             d_reported = d,
             d_min_rounded = min_d_rounded,
             d_max_rounded = max_d_rounded,
             d_in_bounds = d_inbounds)
  })
  
  output$plot_p <- renderPlot({
    res <- recalc_res()
    req(res)
    plot_multiverse_p(res)
  })
  
  output$plot_d <- renderPlot({
    res <- recalc_res()
    req(res)
    plot_multiverse_d(res)
  })
  
  output$downloadReport <- downloadHandler(
    filename = function() {
      paste0("recalc_report_", Sys.Date(), ".html")
    },
    content = function(file) {
      # Ensure report.Rmd is in the working directory used by render()
      src <- normalizePath("report.Rmd")
      owd <- setwd(tempdir())
      on.exit(setwd(owd), add = TRUE)
      file.copy(src, "report.Rmd", overwrite = TRUE)
      
      params <- list(
        m1 = input$m1,
        m2 = input$m2,
        sd1 = input$sd1,
        sd2 = input$sd2,
        n1 = input$n1,
        n2 = input$n2,
        m_digits = input$m_digits,
        sd_digits = input$sd_digits,
        
        p_operator = input$p_operator,
        p = input$p,
        p_digits = input$p_digits,
        p_methods = input$p_methods,
        alpha = input$alpha,
        
        d = input$d,
        d_ci_lower = input$d_ci_lower,
        d_ci_upper = input$d_ci_upper,
        direction = input$direction,
        
        input_rounding = get_input_rounding(),
        output_rounding = get_output_rounding(),
        d_digits = input$d_digits,
        hedges_correction = get_hedges_correction(),
        ci_methods = get_ci_methods(),
        include_se_sd_confusion = isTRUE(input$include_se_sd_confusion)
      )
      
      out <- rmarkdown::render(
        "report.Rmd",
        params = params,
        envir = new.env(parent = globalenv())
      )
      
      file.copy(out, file, overwrite = TRUE)
    }
  )
}

shinyApp(ui = ui, server = server)