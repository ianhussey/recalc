# app.R
# Shiny app for descriptive-statistic consistency between a whole sample and
# its subgroups, backed by:
#   recalc::recalc_total_from_subgroups()    -- E1
#   recalc::recalc_missing_subgroup()        -- E2
# Replaces the previous pooled-SD implementation, which was only correct when
# subgroup means coincide. The recalc package uses the total-variance
# decomposition (within-SS + between-SS) instead.

suppressPackageStartupMessages({
  library(shiny)
  library(DT)
  library(openxlsx)
  library(dplyr)
  library(recalc)
})

# ---------------------------------------------------------------------------
# UI helpers
# ---------------------------------------------------------------------------

subgroup_inputs <- function(prefix, k, compare_sds,
                            default_n = 16, default_m = 4.7, default_sd = 2.4) {
  lapply(seq_len(k), function(i) {
    tagList(
      h5(paste("Subgroup", i)),
      numericInput(paste0(prefix, "_n_",    i), "N",    value = default_n,  step = 1),
      numericInput(paste0(prefix, "_mean_", i), "Mean", value = default_m,  step = 0.01),
      conditionalPanel(
        condition = sprintf("input.%s_compare_sds == true", prefix),
        numericInput(paste0(prefix, "_sd_", i), "SD",   value = default_sd, step = 0.01)
      )
    )
  })
}

# Collect subgroup vectors from numeric inputs. Returns NULL for any vector
# whose entries are still NULL (UI not yet rendered).
collect_subgroups <- function(prefix, k, input, compare_sds) {
  ns    <- vapply(seq_len(k),
                  function(i) input[[paste0(prefix, "_n_",    i)]] %||% NA_real_,
                  numeric(1))
  means <- vapply(seq_len(k),
                  function(i) input[[paste0(prefix, "_mean_", i)]] %||% NA_real_,
                  numeric(1))
  sds <- if (compare_sds) {
    vapply(seq_len(k),
           function(i) input[[paste0(prefix, "_sd_", i)]] %||% NA_real_,
           numeric(1))
  } else NULL
  list(ns = ns, means = means, sds = sds)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Translate the recalc result tibble into a user-friendly data frame for
# display and download:
#   - rename `check` to a plain-English description per tab
#   - drop reported_lower / reported_upper (rounding interval of the reported
#     value; the precision is already shown in the digits inputs)
#   - drop `reported` on tab 2 (E2 has no reported missing-group value)
#   - stringify recalculated_lower / recalculated_upper so NaN renders as
#     "impossible" rather than as an empty DT cell, with numeric values
#     formatted to `display_digits`
prettify_result <- function(res, tab, display_digits) {
  labels_t1 <- c(
    "E1: N = sum n_g"                                                   = "Whole sample N reproduces from subsets",
    "E1: M = (sum n_g M_g) / N"                                         = "Whole sample mean reproduces from subsets",
    "E1: SD = sqrt((sum (n_g-1) s_g^2 + sum n_g (M_g - M)^2) / (N - 1))" = "Whole sample SD reproduces from subsets"
  )
  labels_t2 <- c(
    "E2: n_miss = N - sum n_g"                                              = "Implied missing subgroup N",
    "E2: M_miss = (N M - sum n_g M_g) / n_miss"                             = "Implied missing subgroup mean",
    "E2: SD_miss = sqrt(((N-1) SD^2 - within_g - between_g) / (n_miss - 1))" = "Implied missing subgroup SD"
  )
  labels <- if (tab == "t1") labels_t1 else labels_t2

  fmt <- function(x) {
    ifelse(is.nan(x), "impossible",
           ifelse(is.na(x), "",
                  formatC(x, digits = display_digits, format = "f")))
  }

  out <- data.frame(
    Check               = unname(labels[res$check]),
    `Recalculated lower`  = fmt(res$recalculated_lower),
    `Recalculated upper`  = fmt(res$recalculated_upper),
    Consistent          = res$consistent,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if (tab == "t1") {
    out <- cbind(
      out[, "Check", drop = FALSE],
      Reported = ifelse(is.na(res$reported), "",
                        formatC(res$reported, digits = display_digits,
                                format = "f")),
      out[, c("Recalculated lower", "Recalculated upper", "Consistent"),
          drop = FALSE]
    )
  }
  out
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

ui <- fluidPage(
  titlePanel(
    "ANCHOR: numerical consistency between whole sample and subgroups"
  ),

  tabsetPanel(
    id = "main_tab",

    # ----- Tab 1 ----------------------------------------------------------
    tabPanel(
      "Whole vs subgroups",
      sidebarLayout(
        sidebarPanel(
          h4("Upload / download"),
          fileInput("t1_upload",
                    "Upload .xlsx (download a results file to see the format)",
                    accept = ".xlsx"),
          downloadButton("t1_download", "Download results"),

          hr(),
          h4("Reporting precision"),
          numericInput("t1_n_digits",    "Digits for N",    value = 0, min = 0, step = 1),
          numericInput("t1_mean_digits", "Digits for Mean", value = 2, min = 0, step = 1),
          numericInput("t1_sd_digits",   "Digits for SD",   value = 2, min = 0, step = 1),

          hr(),
          h4("Subgroups"),
          numericInput("t1_k", "Number of subgroups", value = 2, min = 2, step = 1),
          checkboxInput("t1_compare_sds", "Compare SDs?", value = TRUE),

          hr(),
          h4("Overall sample"),
          numericInput("t1_overall_n",    "Overall N",    value = 32,  step = 1),
          numericInput("t1_overall_mean", "Overall Mean", value = 4.7, step = 0.01),
          conditionalPanel(
            condition = "input.t1_compare_sds == true",
            numericInput("t1_overall_sd", "Overall SD", value = 2.4, step = 0.01)
          ),

          hr(),
          uiOutput("t1_subgroup_ui")
        ),
        mainPanel(
          h4("Result"),
          p("Each row checks a reported overall statistic against the value ",
            "recalculated from the subgroups. Rounding intervals are propagated ",
            "through the identity; ", code("consistent == TRUE"),
            " if the reported and recalculated intervals overlap."),
          DTOutput("t1_table"),
          br(),
          tags$details(
            tags$summary("Identities used"),
            withMathJax(HTML(
              "\\(N = \\sum_g n_g\\)<br>",
              "\\(M = \\sum_g n_g M_g \\,/\\, N\\)<br>",
              "\\(\\mathrm{SD}^2 = \\big[\\sum_g (n_g - 1) s_g^2 + ",
              "\\sum_g n_g (M_g - M)^2\\big] \\,/\\, (N - 1)\\)<br>",
              "<br>",
              "The SD identity is the total-variance decomposition ",
              "(within + between). The textbook \"pooled SD\" formula ",
              "\\(\\sqrt{\\sum_g (n_g - 1) s_g^2 / \\sum_g (n_g - 1)}\\) ",
              "equals total SD only when subgroup means coincide; using it ",
              "as a check would false-flag any paper whose subgroup means ",
              "differ."
            ))
          )
        )
      )
    ),

    # ----- Tab 2 ----------------------------------------------------------
    tabPanel(
      "Implied missing subgroup",
      sidebarLayout(
        sidebarPanel(
          h4("Upload / download"),
          fileInput("t2_upload",
                    "Upload .xlsx (download a results file to see the format)",
                    accept = ".xlsx"),
          downloadButton("t2_download", "Download results"),

          hr(),
          h4("Reporting precision"),
          numericInput("t2_n_digits",    "Digits for N",    value = 0, min = 0, step = 1),
          numericInput("t2_mean_digits", "Digits for Mean", value = 2, min = 0, step = 1),
          numericInput("t2_sd_digits",   "Digits for SD",   value = 2, min = 0, step = 1),

          hr(),
          h4("Reported subgroups"),
          numericInput("t2_k", "Number of reported subgroups",
                       value = 1, min = 1, step = 1),
          checkboxInput("t2_compare_sds", "Include SDs?", value = TRUE),

          hr(),
          h4("Overall sample"),
          numericInput("t2_overall_n",    "Overall N",    value = 32,  step = 1),
          numericInput("t2_overall_mean", "Overall Mean", value = 4.7, step = 0.01),
          conditionalPanel(
            condition = "input.t2_compare_sds == true",
            numericInput("t2_overall_sd", "Overall SD", value = 2.4, step = 0.01)
          ),

          hr(),
          h4("Scale bounds (optional)"),
          p(tags$small("If the measurement scale has known logical endpoints ",
                       "(e.g. 1 and 7 for a 7-point Likert), enter both to ",
                       "enable bounds tests on the implied missing-group ",
                       "mean and SD (Bhatia-Davis bound).")),
          numericInput("t2_scale_min", "Scale minimum", value = NA, step = 1),
          numericInput("t2_scale_max", "Scale maximum", value = NA, step = 1),

          hr(),
          uiOutput("t2_subgroup_ui")
        ),
        mainPanel(
          h4("Result"),
          p("Each row gives the implied N, mean, or SD of an unreported ",
            "subgroup, propagated through input rounding. ",
            code("consistent"),
            " reports physical-feasibility checks (rather than an ",
            "overlap check, since there is no reported missing-group value)."),
          DTOutput("t2_table"),
          br(),
          tags$details(
            tags$summary("Identities and feasibility checks used"),
            withMathJax(HTML(
              "\\(n_\\star = N - \\sum_g n_g\\)<br>",
              "\\(M_\\star = (N M - \\sum_g n_g M_g) \\,/\\, n_\\star\\)<br>",
              "\\(s_\\star^2 = \\big[(N - 1) \\mathrm{SD}^2 - ",
              "\\sum_g (n_g - 1) s_g^2 - n_\\star (M_\\star - M)^2 - ",
              "\\sum_g n_g (M_g - M)^2\\big] \\,/\\, (n_\\star - 1)\\)<br>",
              "<br>",
              "<b>Feasibility checks.</b><br>",
              "<i>n</i> row: passes iff the recalculated interval reaches ",
              "\\(\\geq 1\\) (some rounding makes the missing group ",
              "non-empty).<br>",
              "<i>M</i> row: with scale endpoints, passes iff the recalculated ",
              "interval intersects \\([a, b]\\); otherwise NA.<br>",
              "<i>SD</i> row: fails if the propagated interval contains ",
              "NaN (variance negative at some corner). With scale endpoints, ",
              "additionally fails if the implied SD exceeds the ",
              "Bhatia-Davis bound ",
              "\\(\\sqrt{n_\\star/(n_\\star-1)\\,(M_\\star-a)(b-M_\\star)}\\) ",
              "(the largest SD any \\(n_\\star\\)-sample on \\([a,b]\\) with ",
              "mean \\(M_\\star\\) can attain; Bhatia & Davis, 2000)."
            ))
          )
        )
      )
    )
  ),

  tags$hr(),
  tags$div(
    "Developed by Ian Hussey – ",
    tags$a(href  = "https://github.com/ianhussey/recalc",
           "{recalc} GitHub repository", target = "_blank"),
    style = "text-align: center; padding: 10px 0 5px 0;
             font-size: 0.9em; color: #555;"
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  # --- Subgroup-input UIs ---------------------------------------------------
  output$t1_subgroup_ui <- renderUI({
    do.call(tagList, subgroup_inputs("t1", input$t1_k, input$t1_compare_sds))
  })
  output$t2_subgroup_ui <- renderUI({
    do.call(tagList, subgroup_inputs("t2", input$t2_k, input$t2_compare_sds,
                                     default_n = 16, default_m = 4.5,
                                     default_sd = 2.4))
  })

  # --- Tab 1: recalc_total_from_subgroups() --------------------------------
  t1_result <- reactive({
    req(input$t1_k, input$t1_overall_n, input$t1_overall_mean,
        input$t1_n_digits, input$t1_mean_digits)
    sg <- collect_subgroups("t1", input$t1_k, input, input$t1_compare_sds)
    if (any(is.na(sg$ns)) || any(is.na(sg$means))) return(NULL)
    if (input$t1_compare_sds && (is.null(sg$sds) || any(is.na(sg$sds)))) return(NULL)

    args <- list(
      subgroup_ns    = sg$ns,
      subgroup_means = sg$means,
      subgroup_sds   = if (input$t1_compare_sds) sg$sds else NULL,
      overall_n      = input$t1_overall_n,
      overall_mean   = input$t1_overall_mean,
      overall_sd     = if (input$t1_compare_sds) input$t1_overall_sd else NULL,
      n_digits       = input$t1_n_digits,
      mean_digits    = input$t1_mean_digits,
      sd_digits      = if (input$t1_compare_sds) input$t1_sd_digits else NULL
    )
    do.call(recalc::recalc_total_from_subgroups, args)
  })

  t1_display_digits <- reactive({
    max(input$t1_n_digits, input$t1_mean_digits,
        if (input$t1_compare_sds) input$t1_sd_digits else 0)
  })
  t1_pretty <- reactive({
    res <- t1_result(); req(res)
    prettify_result(res, "t1", t1_display_digits())
  })
  output$t1_table <- renderDT({
    datatable(t1_pretty(), options = list(dom = "t"), rownames = FALSE)
  })

  # --- Tab 2: recalc_missing_subgroup() ------------------------------------
  t2_result <- reactive({
    req(input$t2_k, input$t2_overall_n, input$t2_overall_mean,
        input$t2_n_digits, input$t2_mean_digits)
    sg <- collect_subgroups("t2", input$t2_k, input, input$t2_compare_sds)
    if (any(is.na(sg$ns)) || any(is.na(sg$means))) return(NULL)
    if (input$t2_compare_sds && (is.null(sg$sds) || any(is.na(sg$sds)))) return(NULL)

    has_scale <- isTRUE(is.finite(input$t2_scale_min)) &&
                 isTRUE(is.finite(input$t2_scale_max))

    args <- list(
      reported_ns    = sg$ns,
      reported_means = sg$means,
      reported_sds   = if (input$t2_compare_sds) sg$sds else NULL,
      overall_n      = input$t2_overall_n,
      overall_mean   = input$t2_overall_mean,
      overall_sd     = if (input$t2_compare_sds) input$t2_overall_sd else NULL,
      scale_min      = if (has_scale) input$t2_scale_min else NULL,
      scale_max      = if (has_scale) input$t2_scale_max else NULL,
      n_digits       = input$t2_n_digits,
      mean_digits    = input$t2_mean_digits,
      sd_digits      = if (input$t2_compare_sds) input$t2_sd_digits else NULL
    )
    tryCatch(
      suppressWarnings(do.call(recalc::recalc_missing_subgroup, args)),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error")
        NULL
      }
    )
  })

  t2_display_digits <- reactive({
    max(input$t2_n_digits, input$t2_mean_digits,
        if (input$t2_compare_sds) input$t2_sd_digits else 0)
  })
  t2_pretty <- reactive({
    res <- t2_result(); req(res)
    prettify_result(res, "t2", t2_display_digits())
  })
  output$t2_table <- renderDT({
    datatable(t2_pretty(), options = list(dom = "t"), rownames = FALSE)
  })

  # --- Excel I/O (shared format) -------------------------------------------
  # Sheets: "Overall" (Statistic/Value), "Subgroups" (one row per group),
  # optionally "Settings" (Statistic/Value with scale_min/scale_max for tab 2),
  # and (on download) "Results" (the output tibble).

  build_overall_df <- function(tab) {
    if (tab == "t1") {
      df <- data.frame(
        Statistic = c("Overall N", "Overall Mean"),
        Value     = c(input$t1_overall_n, input$t1_overall_mean)
      )
      if (input$t1_compare_sds)
        df <- rbind(df, data.frame(Statistic = "Overall SD",
                                   Value = input$t1_overall_sd))
    } else {
      df <- data.frame(
        Statistic = c("Overall N", "Overall Mean"),
        Value     = c(input$t2_overall_n, input$t2_overall_mean)
      )
      if (input$t2_compare_sds)
        df <- rbind(df, data.frame(Statistic = "Overall SD",
                                   Value = input$t2_overall_sd))
    }
    df
  }

  build_subgroup_df <- function(tab) {
    prefix <- tab
    k <- if (tab == "t1") input$t1_k else input$t2_k
    sg <- collect_subgroups(prefix, k, input,
                            compare_sds = if (tab == "t1") input$t1_compare_sds
                                          else input$t2_compare_sds)
    df <- data.frame(Subgroup = paste0("Subgroup ", seq_len(k)),
                     N    = sg$ns,
                     Mean = sg$means)
    if (!is.null(sg$sds)) df$SD <- sg$sds
    df
  }

  build_settings_df_t2 <- function() {
    data.frame(
      Statistic = c("Digits N", "Digits Mean", "Digits SD",
                    "Scale min", "Scale max"),
      Value = c(input$t2_n_digits, input$t2_mean_digits, input$t2_sd_digits,
                input$t2_scale_min, input$t2_scale_max)
    )
  }

  download_handler <- function(tab) {
    downloadHandler(
      filename = function()
        sprintf("ANCHOR_%s_%s.xlsx",
                if (tab == "t1") "whole_vs_subgroups" else "missing_subgroup",
                Sys.Date()),
      content = function(file) {
        wb <- createWorkbook()
        addWorksheet(wb, "Overall");   writeData(wb, "Overall",   build_overall_df(tab))
        addWorksheet(wb, "Subgroups"); writeData(wb, "Subgroups", build_subgroup_df(tab))
        if (tab == "t2") {
          addWorksheet(wb, "Settings"); writeData(wb, "Settings", build_settings_df_t2())
        }
        res_pretty <- if (tab == "t1") t1_pretty() else t2_pretty()
        if (!is.null(res_pretty)) {
          addWorksheet(wb, "Results")
          writeData(wb, "Results", res_pretty)
        }
        saveWorkbook(wb, file, overwrite = TRUE)
      }
    )
  }
  output$t1_download <- download_handler("t1")
  output$t2_download <- download_handler("t2")

  # Upload handler (shared format reader).
  read_workbook <- function(path, tab) {
    overall   <- read.xlsx(path, sheet = "Overall")
    subgroups <- read.xlsx(path, sheet = "Subgroups")
    get_overall <- function(stat)
      as.numeric(overall$Value[overall$Statistic == stat])
    sd_present <- "Overall SD" %in% overall$Statistic &&
                  "SD" %in% colnames(subgroups)
    k <- nrow(subgroups)

    updateCheckboxInput(session, paste0(tab, "_compare_sds"), value = sd_present)
    updateNumericInput(session,  paste0(tab, "_k"),           value = k)
    updateNumericInput(session,  paste0(tab, "_overall_n"),
                       value = get_overall("Overall N"))
    updateNumericInput(session,  paste0(tab, "_overall_mean"),
                       value = get_overall("Overall Mean"))
    if (sd_present)
      updateNumericInput(session, paste0(tab, "_overall_sd"),
                         value = get_overall("Overall SD"))

    if (tab == "t2" && "Settings" %in% getSheetNames(path)) {
      settings <- read.xlsx(path, sheet = "Settings")
      get_setting <- function(s)
        as.numeric(settings$Value[settings$Statistic == s])
      updateNumericInput(session, "t2_scale_min", value = get_setting("Scale min"))
      updateNumericInput(session, "t2_scale_max", value = get_setting("Scale max"))
    }

    session$onFlushed(function() {
      for (i in seq_len(k)) {
        updateNumericInput(session, paste0(tab, "_n_",    i),
                           value = as.numeric(subgroups$N[i]))
        updateNumericInput(session, paste0(tab, "_mean_", i),
                           value = as.numeric(subgroups$Mean[i]))
        if (sd_present)
          updateNumericInput(session, paste0(tab, "_sd_", i),
                             value = as.numeric(subgroups$SD[i]))
      }
    }, once = TRUE)
  }
  observeEvent(input$t1_upload, { req(input$t1_upload)
    read_workbook(input$t1_upload$datapath, "t1") })
  observeEvent(input$t2_upload, { req(input$t2_upload)
    read_workbook(input$t2_upload$datapath, "t2") })
}

shinyApp(ui = ui, server = server)
