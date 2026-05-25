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
  # devtools::install_github("ianhussey/recalc")
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

  fmt_consistent <- function(x) {
    ifelse(is.na(x), "not checked",
           ifelse(x, "TRUE", "FALSE"))
  }

  out <- data.frame(
    Check               = unname(labels[res$check]),
    `Recalculated lower`  = fmt(res$recalculated_lower),
    `Recalculated upper`  = fmt(res$recalculated_upper),
    Consistent          = fmt_consistent(res$consistent),
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
  # Load MathJax for the identities sections (renderUI(withMathJax(...))
  # re-typesets on each render).
  withMathJax(),

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
          selectInput(
            "t1_rounding",
            "Rounding mode for reported values",
            choices = c("Either half-up or bankers (default)" = "either",
                        "Round half up" = "half_up",
                        "Bankers (round half to even)" = "bankers",
                        "Truncated" = "truncate"),
            selected = "either"
          ),

          hr(),
          h4("Subgroups"),
          numericInput("t1_k", "Number of subgroups", value = 2, min = 2, step = 1),
          checkboxInput("t1_compare_sds", "Include SDs?", value = TRUE),

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
          p(tags$small(
            "See the ", tags$em("Guidance"),
            " tab for the formulas, notation, and feasibility-check ",
            "definitions used by both tabs."
          ))
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
          selectInput(
            "t2_rounding",
            "Rounding mode for reported values",
            choices = c("Either half-up or bankers (default)" = "either",
                        "Round half up" = "half_up",
                        "Bankers (round half to even)" = "bankers",
                        "Truncated" = "truncate"),
            selected = "either"
          ),

          hr(),
          h4("Subgroups"),
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
          p(tags$small(
            "See the ", tags$em("Guidance"),
            " tab for the formulas, notation, and feasibility-check ",
            "definitions used by both tabs."
          ))
        )
      )
    ),

    # ----- Guidance tab ----------------------------------------------------
    tabPanel(
      "Guidance",
      div(
        style = "max-width: 900px; margin: 0 auto; padding: 1.5em 1em;",
        h4("Guidance on app usage"),
        p(
          "The app has two use cases, partially nested. Most users will ",
          "want to start with the first, and only proceed to the second if ",
          "a subgroup appears under-reported or missing."
        ),
        tags$ol(
          tags$li(
            tags$b("Whole sample vs subgroups (Tab 1)."),
            " Enter the reported overall ", tags$em("N"), ", ",
            tags$em("mean"), ", and (optionally) ", tags$em("SD"),
            " of a sample, plus the per-subgroup ", tags$em("n"), ", ",
            tags$em("mean"), ", and ", tags$em("SD"), " for each reported ",
            "subgroup. The app checks whether the reported overall ",
            "statistics could have been computed from the reported ",
            "subgroups, propagating the rounding interval of every input ",
            "through the aggregation identities (see Tab 1's ",
            tags$em("Identities used"),
            " for the formulas). A flagged row indicates that no ",
            "rounding of the inputs is consistent with the reported ",
            "overall value — i.e., the table is internally incoherent."
          ),
          tags$li(
            tags$b("Implied missing subgroup (Tab 2)."),
            " If Tab 1 surfaces an inconsistency, or if a paper otherwise ",
            "suggests that one or more subgroups are missing or under-",
            "reported (for example, a third gender category like ",
            "gender-diverse participants who are mentioned in the text but ",
            "do not appear in the breakdown), enter the overall and the ",
            "reported subgroups in Tab 2. The app then derives the ",
            tags$em("implied"), " ", tags$em("n"), ", ", tags$em("mean"),
            ", and ", tags$em("SD"), " that the missing subgroup would ",
            "have to have for the reported overall to be self-consistent. ",
            "Supplying the scale's logical minimum and maximum ",
            "(e.g., 1 and 7 for a 7-point Likert response) activates ",
            "additional plausibility checks on the implied mean and SD ",
            "via the Bhatia-Davis bound."
          )
        ),
        p(
          tags$b("Rounding mode."),
          " By default the app assumes reported values were rounded under ",
          "either round-half-up or bankers (round-half-to-even). For a ",
          "paper that specifies a different convention, use the ",
          tags$em("Rounding mode"), " dropdown on the relevant tab. The ",
          tags$code("Truncated"), " option truncated values toward zero. ",
          "This produces wider compatability intervals, but is uncommon ",
          "in the literature."
        ),
        p(
          tags$b("Upload / download."),
          " Each tab can download a results spreadsheet that records the ",
          "inputs, the settings (digits and rounding mode), and the ",
          "consistency table. Uploading that spreadsheet on the same tab ",
          "restores the full app state, so a results file functions as a ",
          "shareable template for re-running the check."
        ),
        hr(),
        h4("Math"),
        uiOutput("math_section")
      )
    ),

    # ----- About tab -------------------------------------------------------
    tabPanel(
      "About",
      div(
        style = "max-width: 900px; margin: 0 auto; padding: 1.5em 1em;",
        h4("About"),
        p(
          "The checks in this app correspond to ",
          tags$a(
            href = "https://inspect.sr/chapters/check_4_10.html",
            target = "_blank",
            "INSPECT-SR check 4.10"
          ),
          " (consistency between whole-sample and subgroup summary ",
          "statistics). INSPECT-SR (",
          tags$a(
            href = "https://www.medrxiv.org/content/10.1101/2025.09.03.25334905v3",
            target = "_blank",
            "Wilkinson et al., 2025"
          ),
          ") is a framework for assessing the trustworthiness of ",
          "randomised controlled trials in systematic reviews; the same ",
          "checks apply equally outside the RCT setting."
        ),
        p(
          "The app is a thin Shiny front-end to the ",
          tags$code("recalc_total_from_subgroups()"), " and ",
          tags$code("recalc_missing_subgroup()"),
          " functions of the ", tags$code("{recalc}"),
          " R package. Both functions implement the analysis-of-variance ",
          "sum-of-squares decomposition exactly (with rounding-aware ",
          "interval propagation through the inputs); Tab 2's ",
          "Bhatia-Davis bound on the implied missing-subgroup SD is from ",
          "Bhatia & Davis (2000)."
        ),
        p(
          "Source code is available ",
          tags$a(
            href = "https://github.com/ianhussey/recalc",
            target = "_blank",
            "on GitHub"
          ),
          ". Bug reports, suggestions, and pull requests are welcome."
        ),
        h5("References"),
        tags$ul(
          tags$li(
            "Bhatia, R., & Davis, C. (2000). A better bound on the variance. ",
            tags$em("American Mathematical Monthly"), ", 107(4), 353–357."
          ),
          tags$li(
            "Wilkinson, J., et al. (2025). INSPECT-SR development protocol. ",
            tags$em("BMJ Open"), ", 14, e084164. ",
            tags$a(
              href = "https://inspect.sr",
              target = "_blank",
              "inspect.sr"
            )
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
      sd_digits      = if (input$t1_compare_sds) input$t1_sd_digits else NULL,
      rounding       = input$t1_rounding
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

  # Consolidated math section (rendered in the Guidance tab) covering the
  # notation, the Tab 1 aggregation identities, the pooled-SD counterfactual,
  # and the Tab 2 missing-subgroup identities + feasibility checks.
  output$math_section <- renderUI({
    withMathJax(HTML(
      # ----- Notation -----
      "<h5>Notation</h5>",
      "<ul>",
      "<li>\\(N\\), \\(M\\): whole-sample size and mean. ",
      "\\(s^2\\): whole-sample variance.</li>",
      "<li>\\(n_g\\), \\(M_g\\), \\(s_g^2\\): subgroup \\(g\\)'s size, mean, ",
      "and variance, for \\(g = 1, \\ldots, k\\).</li>",
      "<li>\\(n_\\star\\), \\(M_\\star\\), \\(s_\\star^2\\): implied size, ",
      "mean, and variance of the unreported subgroup (Tab 2 only).</li>",
      "<li>\\([a, b]\\): logical scale endpoints (Tab 2, optional).</li>",
      "</ul>",
      "<p>The corresponding SD is the square root of each variance. ",
      "Formulas are given for variances throughout to keep MathJax ",
      "rendering robust across browsers.</p>",

      # ----- Tab 1 identities -----
      "<h5>Tab 1: whole sample reproduces from subgroups</h5>",
      "$$N = \\sum_g n_g, \\qquad M = \\frac{\\sum_g n_g M_g}{N}$$",
      "$$s^2 = \\frac{\\sum_g (n_g - 1)\\, s_g^2 + ",
      "\\sum_g n_g\\, (M_g - M)^2}{N - 1}$$",
      "<p>The SD identity is the total-variance decomposition ",
      "(within-group SS + between-group SS). Each row's ",
      "<code>consistent</code> flag is <code>TRUE</code> when the rounding ",
      "interval of the reported overall value overlaps the rounding-aware ",
      "interval implied by the subgroups.</p>",

      # ----- Pooled-SD counterfactual -----
      "<h5>Why not the textbook \"pooled SD\" formula?</h5>",
      "$$s_\\text{pooled}^2 = \\frac{\\sum_g (n_g - 1)\\, s_g^2}",
      "{\\sum_g (n_g - 1)}$$",
      "<p>This is the equal-variance \\(t\\)-test's pooled variance ",
      "estimator and the denominator of Cohen's \\(d\\). It equals the ",
      "whole-sample variance only when subgroup means coincide ",
      "(\\(SS_\\text{between} = 0\\)); otherwise it understates by ",
      "\\(SS_\\text{between} / (N - 1)\\). Using it as a consistency ",
      "check would false-flag any table whose subgroups have different ",
      "means.</p>",

      # ----- Tab 2 identities -----
      "<h5>Tab 2: implied missing subgroup</h5>",
      "$$n_\\star = N - \\sum_g n_g, \\qquad ",
      "M_\\star = \\frac{N M - \\sum_g n_g M_g}{n_\\star}$$",
      "$$s_\\star^2 = \\frac{(N - 1)\\, s^2 - ",
      "\\sum_g (n_g - 1)\\, s_g^2 - n_\\star\\, (M_\\star - M)^2 - ",
      "\\sum_g n_g\\, (M_g - M)^2}{n_\\star - 1}$$",
      "<p>Tab 2's <code>consistent</code> flags are feasibility checks ",
      "(rather than overlap checks, since the missing subgroup is not ",
      "reported):</p>",
      "<ul>",
      "<li><i>n</i> row: passes iff the recalculated interval reaches ",
      "\\(\\geq 1\\).</li>",
      "<li><i>M</i> row: with scale endpoints, passes iff the recalculated ",
      "interval intersects \\([a, b]\\); otherwise shown as \"not ",
      "checked\".</li>",
      "<li><i>SD</i> row: fails if the propagated interval contains NaN ",
      "(variance strictly negative at some corner). With scale endpoints, ",
      "additionally fails if the implied SD exceeds the Bhatia-Davis ",
      "upper bound:</li>",
      "</ul>",
      "$$s_\\star^2 \\;\\le\\; \\frac{n_\\star}{n_\\star - 1}\\,",
      "(M_\\star - a)(b - M_\\star)$$",
      "<p>i.e. the largest variance any \\(n_\\star\\)-sample on ",
      "\\([a, b]\\) with mean \\(M_\\star\\) can attain ",
      "(Bhatia &amp; Davis, 2000).</p>"
    ))
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
      sd_digits      = if (input$t2_compare_sds) input$t2_sd_digits else NULL,
      rounding       = input$t2_rounding
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

  # (Tab 1 and Tab 2 identities are rendered in the Guidance tab via
  # output$math_section above.)

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

  # Build the Settings sheet. Shared across tabs; Tab 2 additionally writes
  # the scale-bound entries. Character type so "either"/"half_up"/etc. survive
  # alongside the numeric digit values when openxlsx writes the column.
  build_settings_df <- function(tab) {
    stats_common <- c("Digits N", "Digits Mean", "Digits SD", "Rounding")
    vals_common <- if (tab == "t1") {
      c(input$t1_n_digits, input$t1_mean_digits, input$t1_sd_digits,
        input$t1_rounding)
    } else {
      c(input$t2_n_digits, input$t2_mean_digits, input$t2_sd_digits,
        input$t2_rounding)
    }
    df <- data.frame(Statistic = stats_common,
                     Value     = as.character(vals_common),
                     stringsAsFactors = FALSE)
    if (tab == "t2") {
      df <- rbind(df, data.frame(
        Statistic = c("Scale min", "Scale max"),
        Value     = as.character(c(input$t2_scale_min, input$t2_scale_max)),
        stringsAsFactors = FALSE
      ))
    }
    df
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
        addWorksheet(wb, "Settings"); writeData(wb, "Settings", build_settings_df(tab))
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

    # Restore Settings (digits + rounding, and on Tab 2 also scale bounds).
    # Old downloads that predate the Settings sheet are still accepted: any
    # missing entries are simply left at their current input values.
    if ("Settings" %in% getSheetNames(path)) {
      settings <- read.xlsx(path, sheet = "Settings")
      get_setting <- function(s, parse = c("numeric", "character")) {
        parse <- match.arg(parse)
        v <- settings$Value[settings$Statistic == s]
        if (length(v) == 0L) return(NULL)
        if (parse == "numeric") as.numeric(v) else as.character(v)
      }
      n_d  <- get_setting("Digits N");    if (!is.null(n_d))  updateNumericInput(session, paste0(tab, "_n_digits"),    value = n_d)
      m_d  <- get_setting("Digits Mean"); if (!is.null(m_d))  updateNumericInput(session, paste0(tab, "_mean_digits"), value = m_d)
      s_d  <- get_setting("Digits SD");   if (!is.null(s_d))  updateNumericInput(session, paste0(tab, "_sd_digits"),   value = s_d)
      rnd  <- get_setting("Rounding", "character")
      if (!is.null(rnd) && !is.na(rnd) && nzchar(rnd))
        updateSelectInput(session, paste0(tab, "_rounding"), selected = rnd)
      if (tab == "t2") {
        smin <- get_setting("Scale min"); if (!is.null(smin)) updateNumericInput(session, "t2_scale_min", value = smin)
        smax <- get_setting("Scale max"); if (!is.null(smax)) updateNumericInput(session, "t2_scale_max", value = smax)
      }
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
