library(shiny)
library(rhandsontable)

# main function -----------------------------------------------------------

base_p_cont <- function(df) {
  a <- df

  a <- a %>% dplyr::rename(study = Study, var = Variable)

  #helper to select columns - nm = m& s, nm.x = dec place, nm_x = m & s
  nm <- colnames(a)[
    substr(colnames(a), 1, 1) == "m" |
      (substr(colnames(a), 1, 1) == "s" & colnames(a) != "study")
  ]
  nm_m <- nm[substr(nm, 1, 1) == "m"]
  nm_s <- nm[substr(nm, 1, 1) == "s"]
  nm_n <- colnames(a)[substr(colnames(a), 1, 1) == "n"]
  nm.m <- paste0("dp", nm_m)
  nm.s <- paste0("dp", nm_s)
  nm.sf <- paste0(nm, "_sf")

  #make sure all numeric variables are correct types, but leave p as is.
  a <- a %>% dplyr::mutate(across(all_of(c(nm_m, nm_s, nm_n)), as.numeric))
  a$p <- ifelse(is.na(a$p), "", a$p)
  a$p_num <- as.numeric(a$p)
  a$p_dec.place <- dp(a$p_num)

  #remove impossible rows
  a$missing <- rowSums(is.na(a[nm]))
  missing <- a %>% dplyr::filter(a$missing > 0)
  a <- a %>% dplyr::filter(a$missing == 0) %>% dplyr::select(-missing)
  am <- a %>% dplyr::select(-p_num, -p_dec.place)

  #check pvals are possible
  if (sum(!is.na(a$p)) > 0) {
    a <- chk_pvals(dat = a, am, nm, nm_m, nm_s, nm_n, nm.m, nm.s, nm.sf)
  }

  max_digits <- ifelse(
    sum(!is.na(a$p_dec.place)) > 0,
    max(a$p_dec.place, na.rm = T),
    3
  )

  #add back in missing
  a <- dplyr::bind_rows(
    a,
    missing %>%
      dplyr::select(study, var, p, p_num) %>%
      dplyr::mutate(p_flag = "not done because missing data")
  )
  #change names
  a <- a %>%
    dplyr::rename(
      reported_pvalue = p,
      variable = var,
      max_possible_p = p_max,
      min_possible_p = p_min,
      pvalue_flag = p_flag
    ) %>%
    dplyr::select(-all_of(c(nm_m, nm_s, nm_n)), -p_num, -p_dec.place)

  return(list(a, max_digits))
}

# helper functions --------------------------------------------------------

#generic function to make into long format
lng_f <- function(dat, remvars) {
  b <- dat %>%
    tidyr::gather(v, value, -c(all_of(remvars))) %>%
    tidyr::separate_wider_regex(v, c(var1 = ".", group = "\\d+")) %>%
    tidyr::spread(var1, value) %>% #delete missing data
    dplyr::filter(!is.na(m) & !is.na(s) & !is.na(n))
  return(b)
}

#calculate the pvalues
pval <- function(df, simulate = "n") {
  if (simulate == "n") {
    df <- df %>% dplyr::group_by(study, var)
  } else {
    df <- df %>% dplyr::group_by(sim, study, var)
  }

  df <- df %>%
    dplyr::summarise(
      # Get the # of means and overall mean
      k = dplyr::n(),
      tot_n = sum(n),
      ov_m = sum(n * m) / tot_n,
      # Calculate mean sq between and within
      msb = sum(n * (m - ov_m)^2),
      msw = sum((n - 1) * s^2),
      .groups = "drop"
    ) %>%
    # Filter out cases with insufficient data
    dplyr::filter(k >= 2, tot_n > k) %>%
    dplyr::mutate(
      # Degrees of freedom
      dfb = k - 1,
      dfw = tot_n - k,
      # F-statistic and p-value
      f_value = ifelse(msw > 0, msb / dfb / (msw / dfw), NA),
      p.value = ifelse(
        !is.na(f_value),
        pf(f_value, dfb, dfw, lower.tail = FALSE),
        NA
      )
    )

  if (simulate == "n") {
    df <- df %>% dplyr::select(study, var, p.value)
  } else {
    df <- df %>% dplyr::select(study, var, sim, p.value)
  }
}

#function to count decimal places
dec_places <- function(a, nm, nm.m, nm.s) {
  a1 <- cbind(a, setNames(lapply(a[nm], dp), paste0("dp", nm)))

  #r and excel don't include final 0s (as numbers) so take maximum dp per variable (assuming that if m1=x and m2 = x.1, m1 is actually x.0)
  m <- a1[nm.m]
  m[] <- t(apply(m, 1, function(x) ifelse(is.na(x), NA, max(x, na.rm = TRUE))))
  s <- a1[nm.s]
  s[] <- t(apply(s, 1, function(x) ifelse(is.na(x), NA, max(x, na.rm = TRUE))))

  # Replace the columns back in original dataframe
  a1[nm.m] <- m
  a1[nm.s] <- s

  return(a1)
}

#basic dec places fn
dp <- function(x) {
  ifelse(
    abs(x - round(x)) > .Machine$double.eps^0.5,
    nchar(sub('^\\d+\\.', '', sub('0+$', '', as.character(abs(x))))),
    0
  )
}

#sig figures function
sigfigs <- function(j, dat) {
  x <- dat[, j]

  sapply(x, function(val) {
    if (is.na(val) || !is.finite(val)) {
      return(NA)
    }
    if (val == 0) {
      return(1L)
    }

    # sprintf
    str_val <- sprintf("%g", abs(val))
    cleaned <- gsub("0+$", "", gsub("^0+", "", gsub("\\.", "", str_val)))

    if (cleaned == "") {
      return(1L)
    }
    nchar(cleaned)
  })
}

#function to work out max and min p-value
rng_p_fn <- function(dat, nm, nm_m, nm_s, nm_n, nm.m, nm.s, nm.sf) {
  am <- dat
  #get dp and sig figs
  rng <- dec_places(am, nm, nm.m, nm.s)
  rng[, nm.sf] <- lapply(nm, function(x) sigfigs(x, rng))

  #now get max values to round
  nm.dtm <- paste0(nm_m, "dt")
  nm.dts <- paste0(nm_s, "dt")
  rng[, c(nm.dtm, nm.dts)] <- rng[, c(nm_m, nm_s)]

  #now do each one separately
  rng <- rng %>% #get max values for each
    dplyr::mutate(
      dpm = do.call(pmax, c(dplyr::select(., all_of(nm.m)), na.rm = T)),
      dps = do.call(pmax, c(dplyr::select(., all_of(nm.s)), na.rm = T)),
      sfm = do.call(pmax, c(dplyr::select(., paste0(nm_m, "_sf")), na.rm = T)),
      sfs = do.call(pmax, c(dplyr::select(., paste0(nm_s, "_sf")), na.rm = T))
    ) %>%
    dplyr::select(-c(all_of(nm.m), all_of(nm.s), all_of(nm.sf)))

  #round and select max value for rounding
  rng <- rng %>%
    dplyr::mutate(
      across(all_of(nm.dtm), ~ max_round(.x, dpm, sfm)),
      across(all_of(nm.dts), ~ max_round(.x, dps, sfs))
    ) %>%
    dplyr::mutate(
      dfm = do.call(pmax, c(dplyr::select(., all_of(nm.dtm)), na.rm = T)),
      dfs = do.call(pmax, c(dplyr::select(., all_of(nm.dts)), na.rm = T))
    )

  #now merge back into main data
  b <- lng_f(am %>% dplyr::select(-p), remvars = c("study", "var"))
  b <- b %>%
    dplyr::left_join(
      rng %>% dplyr::select(study, var, dfm),
      by = c("study", "var")
    )

  #work out which combinations have the biggest and smallest mean square difference
  b2 <- ms(b)

  #analyse separately note max = max diff ie min p and vice versa
  b_max <- dplyr::left_join(
    rng %>% dplyr::select(-all_of(nm_m)),
    b2 %>%
      dplyr::filter(stat == "max") %>%
      dplyr::select(study, var, all_of(nm_m)),
    by = c("study", "var")
  ) %>%
    dplyr::select(-c(p, all_of(nm.dtm), all_of(nm.dts), dpm:sfs)) %>%
    dplyr::mutate(across(all_of(nm_s), ~ .x - dfs))

  b_min <- dplyr::left_join(
    rng %>% dplyr::select(-all_of(nm_m)),
    b2 %>%
      dplyr::filter(stat == "min") %>%
      dplyr::select(study, var, all_of(nm_m)),
    by = c("study", "var")
  ) %>%
    dplyr::select(-c(p, all_of(nm.dtm), all_of(nm.dts), dpm:sfs)) %>%
    dplyr::mutate(across(all_of(nm_s), ~ .x + dfs))

  #get max and min pvals
  b_max1 <- lng_f(b_max, remvars = c("study", "var", "dfs", "dfm"))
  b_max1 <- pval(b_max1)

  b_min1 <- lng_f(b_min, remvars = c("study", "var", "dfs", "dfm"))
  b_min1 <- pval(b_min1)

  #range of pvals
  rng_data <- dplyr::bind_rows(
    dplyr::left_join(
      b_max,
      b_max1 %>%
        dplyr::rename(p_min = p.value) %>%
        dplyr::mutate(stat = "min_p-value"),
      by = c("study", "var")
    ),
    dplyr::left_join(
      b_min,
      b_min1 %>%
        dplyr::rename(p_max = p.value) %>%
        dplyr::mutate(stat = "max_p-value"),
      by = c("study", "var")
    )
  ) %>%
    dplyr::mutate(p_range = dplyr::coalesce(p_min, p_max)) %>%
    dplyr::select(
      study,
      var,
      stat,
      all_of(nm_m),
      all_of(nm_s),
      all_of(nm_n),
      dfm,
      dfs,
      p_range,
      p_max,
      p_min
    )

  return(rng_data)
}


#determine combinations of means with max and min weighted mean square error
ms <- function(data) {
  data %>%
    dplyr::group_by(study, var) %>%
    dplyr::group_modify(
      ~ {
        # Get values
        m_vals <- .x$m
        n_vals <- .x$n
        df_vals <- .x$dfm

        # Create all combinations of m+df and m-df for each observation
        # Each observation can be either m+df or m-df
        combs <- expand.grid(lapply(1:length(m_vals), function(i) {
          c(m_vals[i] + df_vals[i], m_vals[i] - df_vals[i])
        }))

        # Calculate mean square
        ms <- apply(combs, 1, function(m_comb) {
          tot_n <- sum(n_vals)
          ov_m <- sum(n_vals * m_comb) / tot_n
          ms <- sum(n_vals * (m_comb - ov_m)^2)
          return(ms)
        })

        m_cols <- combs
        colnames(m_cols) <- paste0("m", 1:ncol(combs))

        df <- data.frame(
          ms = ms,
          m_cols
        ) %>%
          dplyr::mutate(
            stat = dplyr::case_when(
              ms == max(ms) ~ "max",
              ms == min(ms) ~ "min",
              .default = ""
            )
          ) %>%
          dplyr::filter(stat %in% c("max", "min")) %>%
          dplyr::distinct(ms, stat, .keep_all = TRUE)
      }
    )
}

# work out how much to round values
max_round <- function(x, dpm, sig_figs) {
  dat <- data.frame(x = x, dpm = dpm, sig_figs = sig_figs) %>%
    dplyr::mutate(
      z = dplyr::case_when(
        is.na(x) ~ NA,
        x == 0 ~ 0.5,
        abs(x) >= 1 ~ {
          magnitude <- floor(log10(abs(x)))
          least_sig_position <- magnitude - sig_figs + 1
          max_val <- abs(x) + 0.5 * 10^least_sig_position
          max_val - abs(x)
        },
        abs(x) < 1 ~ 0.5 * 10^(-dpm),
        .default = x
      )
    )
  return(dat$z)
}

#check whether to flag pvalues
chk_pvals <- function(dat = a, am, nm, nm_m, nm_s, nm_n, nm.m, nm.s, nm.sf) {
  #find min and max pvalues
  rng_p_data <- rng_p_fn(dat = am, nm, nm_m, nm_s, nm_n, nm.m, nm.s, nm.sf)

  a1 <- dat %>%
    dplyr::left_join(
      rng_p_data %>%
        dplyr::filter(stat == "min_p-value") %>%
        dplyr::select(study, var, p_min),
      by = c("study", "var")
    ) %>%
    dplyr::left_join(
      rng_p_data %>%
        dplyr::filter(stat == "max_p-value") %>%
        dplyr::select(study, var, p_max),
      by = c("study", "var")
    )
  # Check whether p-values are possible
  a1 <- a1 %>%
    dplyr::mutate(
      p_dec.place = dplyr::coalesce(p_dec.place, 3),
      across(p_min:p_max, ~ round(.x, p_dec.place)),
      p_num_min = dplyr::case_when(
        grepl(">", p) ~ suppressWarnings(as.numeric(sub(".*>", "", p))),
        grepl("ns", tolower(p)) ~ 0.05,
        is.na(p_num) ~ NA,
        .default = p_min
      ),
      p_num_max = dplyr::case_when(
        grepl("<", p) ~ suppressWarnings(as.numeric(sub(".*<", "", p))),
        is.na(p_num) ~ NA,
        .default = p_max
      )
    ) %>%
    #now put flag for impossible results
    dplyr::mutate(
      p_flag = dplyr::case_when(
        is.na(p_num) & is.na(p_num_min) & is.na(p_num_max) ~ NA,
        p_num >= p_min & p_num <= p_max ~ "",
        p_num < p_min | p_num > p_max ~ "flag",
        p_num_max < p_min | p_num_min > p_max ~ "flag",
        .default = ""
      )
    ) %>%
    dplyr::select(-c(p_num_min, p_num_max))
  return(a1)
}


# shiny app ---------------------------------------------------------------

ui <- fluidPage(
  titlePanel(
    "Check p-values for baseline continuous variables from a randomised controlled trial"
  ),

  sidebarLayout(
    sidebarPanel(
      numericInput("groups", "Number of Groups:", value = 2, min = 1),
      numericInput("rows", "Number of Rows:", value = 10, min = 1),
      radioButtons(
        "nN_pattern",
        "Choose n and N arrangement:",
        choices = list(
          "n1, n2, n3 ... m1, m2, m3 ... s1, s2, s3 ..." = "nms",
          "m1, m2, .... s1, s2, ... n1, n2" = "msn",
          "n1, m1, s1, n2, m2, s2 ..." = "n1m1s1",
          "m1, s1, n1, m2, s2, n2 ..." = "m1s1n1"
        ),
        selected = "msn"
      ),
      actionButton("create_table", "Create Table"),
      br(),
      br(),
      actionButton("analyze", "Run Analysis"),
      br(),
      br(),
      downloadButton("download", "Download Results"),
      br(),
      br(),
      downloadButton("download_data", "Download Data"),
      br(),
      br(),
      actionButton("reset", "Reset App", width = "100%"),

      div(
        style = "background-color: #f8f9fa; padding: 10px; margin-top: 20px; border-radius: 5px; font-size: 12px;",
        strong("Feedback"),
        br(),
        "For bugs email ",
        a("m.bolland@auckland.ac.nz", href = "mailto:m.bolland@auckland.ac.nz"),
        br(),
        "I'd also love to hear if it has been useful for you."
      )
    ),

    mainPanel(
      rhandsontable::rHandsontableOutput("data_table"),
      conditionalPanel(
        condition = "output.table_created",
        p(
          style = "color: #666; margin-top: 10px;",
          "Fill in all cells with your data.",
          br(),
          "You can copy and paste from Excel or type directly into the cells. You can cut and paste multiple rows at once",
          br(),
          "Empty cells are allowed for Variable, and P only, otherwise an error will occur",
          br(),
          "Put the columns in your selected order and use the button to make sure the order in the table matches this order",
          br(),
          ""
        )
      ),
      br(),
      # Display results dataframe
      conditionalPanel(
        condition = "output.results_available",
        h3("Analysis Results"),
        tableOutput("results_table")
      )
    )
  )
)

server <- function(input, output, session) {
  # Reactive values to store data
  values <- reactiveValues(table_data = NULL, analysis_results = NULL)

  #activate once table created
  output$table_created <- reactive({
    !is.null(values$table_data)
  })
  outputOptions(output, "table_created", suspendWhenHidden = FALSE)

  get_column_order <- reactive({
    n_groups <- input$groups
    pattern <- input$nN_pattern

    col_names <- c("Study", "Variable")

    if (pattern == "nms") {
      col_names <- c(
        col_names,
        paste0("n", 1:n_groups),
        paste0("m", 1:n_groups),
        paste0("s", 1:n_groups),
        "p"
      )
    } else if (pattern == "msn") {
      col_names <- c(
        col_names,
        paste0("m", 1:n_groups),
        paste0("s", 1:n_groups),
        paste0("n", 1:n_groups),
        "p"
      )
    } else if (pattern == "n1m1s1") {
      col_names <- c(
        col_names,
        as.vector(sapply(1:n_groups, function(i) {
          c(paste0("n", i), paste0("m", i), paste0("s", i))
        })),
        "p"
      )
    } else if (pattern == "m1s1n1") {
      col_names <- c(
        col_names,
        as.vector(sapply(1:n_groups, function(i) {
          c(paste0("m", i), paste0("s", i), paste0("n", i))
        })),
        "p"
      )
    }
    return(col_names)
  })

  # Create dynamic table structure
  observeEvent(input$create_table, {
    req(get_column_order)
    columns <- get_column_order()

    # Create empty dataframe
    values$table_data <- data.frame(matrix(
      "",
      nrow = input$rows,
      ncol = length(columns)
    ))
    names(values$table_data) <- columns
  })

  # Render the editable table
  output$data_table <- rhandsontable::renderRHandsontable({
    if (!is.null(values$table_data)) {
      rhandsontable::rhandsontable(values$table_data, stretchH = "all")
    }
  })

  # Update data when table is edited
  observe({
    if (!is.null(input$data_table)) {
      values$table_data <- hot_to_r(input$data_table)
    }
  })

  # Run analysis
  observeEvent(input$analyze, {
    # Clear previous results first
    values$analysis_results <- NULL

    if (!is.null(values$table_data)) {
      # Check for missing values
      values$table_data[values$table_data == ""] <- NA

      tryCatch(
        {
          results <- base_p_cont(as.data.frame(values$table_data))
          values$analysis_results <- results[[1]]
          values$max_digits <- results[[2]]
        },

        error = function(e) {
          values$analysis_results <- data.frame(
            Error = paste("Error in analysis:", e$message)
          )
          values$max_digits <- 3
        }
      )
    }
  })

  # Display the results dataframe
  output$results_table <- renderTable(
    {
      values$analysis_results
    },
    digits = function() values$max_digits
  )

  output$results_available <- reactive({
    !is.null(values$analysis_results)
  })
  outputOptions(output, "results_available", suspendWhenHidden = FALSE)

  # Download handler
  output$download <- downloadHandler(
    filename = function() {
      paste("analysis_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(values$analysis_results, file, row.names = FALSE)
    }
  )

  output$download_data <- downloadHandler(
    filename = function() {
      paste("analysis_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(
        as.data.frame(values$table_data),
        file,
        row.names = FALSE,
        na = ""
      )
    }
  )

  observeEvent(input$reset, {
    session$reload()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
