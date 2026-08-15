library(shiny)
library(emdi)
library(dplyr)
library(DT)
library(bslib)
library(fontawesome)
library(bsicons)
library(ggplot2)

source("Bootstrap_emdi_new_function.R")

# --- Create custom emdi environment ---
emdi2_env <- new.env(parent = asNamespace("emdi"))

# Funzioni custom già create
emdi2_env$mse_estim2 <- mse_estim2
emdi2_env$mse_estim_wrapper2 <- mse_estim_wrapper2
emdi2_env$parametric_bootstrap2 <- parametric_bootstrap2
emdi2_env$ebp2 <- ebp2

# Nuove funzioni di trasformazione
emdi2_env$allowed_transformations_ebp2 <- allowed_transformations_ebp2
emdi2_env$check_unit_interval <- check_unit_interval
emdi2_env$clip_open01 <- clip_open01
emdi2_env$clip_closed01 <- clip_closed01

emdi2_env$logit_transform <- logit_transform
emdi2_env$logit_transform_back <- logit_transform_back
emdi2_env$probit_transform <- probit_transform
emdi2_env$probit_transform_back <- probit_transform_back
emdi2_env$arcsin_transform <- arcsin_transform
emdi2_env$arcsin_transform_back <- arcsin_transform_back

emdi2_env$data_transformation <- data_transformation
emdi2_env$std_data_transformation <- std_data_transformation
emdi2_env$back_transformation <- back_transformation
emdi2_env$optimal_parameter <- optimal_parameter
emdi2_env$ebp_check2 <- ebp_check2

# Funzioni interne originali di emdi che devono usare il nuovo ambiente
emdi2_env$point_estim <- getFromNamespace("point_estim", "emdi")
emdi2_env$monte_carlo <- getFromNamespace("monte_carlo", "emdi")
emdi2_env$prediction_y <- getFromNamespace("prediction_y", "emdi")
emdi2_env$errors_gen <- getFromNamespace("errors_gen", "emdi")
emdi2_env$model_par <- getFromNamespace("model_par", "emdi")
emdi2_env$gen_model <- getFromNamespace("gen_model", "emdi")

emdi2_env$superpopulation <- getFromNamespace("superpopulation", "emdi")
emdi2_env$superpopulation_wild <- getFromNamespace("superpopulation_wild", "emdi")
emdi2_env$bootstrap_par <- getFromNamespace("bootstrap_par", "emdi")
emdi2_env$bootstrap_par_wild <- getFromNamespace("bootstrap_par_wild", "emdi")

# Imposta l'ambiente custom per tutte le funzioni
for (nm in names(emdi2_env)) {
  if (is.function(emdi2_env[[nm]])) {
    environment(emdi2_env[[nm]]) <- emdi2_env
  }
}

# Funzione finale da usare nell'app
ebp2 <- emdi2_env$ebp2


# ============================================================
# Utility functions
# ============================================================

read_uploaded_data <- function(file) {
  ext <- tools::file_ext(file$name)
  
  if (ext == "csv") {
    read.csv(file$datapath, stringsAsFactors = TRUE)
  } else if (ext == "rds") {
    readRDS(file$datapath)
  } else {
    stop("Only .csv and .rds files are supported.")
  }
}

safe_prob <- function(p, B) {
  eps <- 1 / (2 * B)
  pmin(pmax(p, eps), 1 - eps)
}

method_suffix <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

indicator_support <- function(indicator) {
  if (indicator %in% c("Head_Count", "Gini")) {
    c(0, 1)
  } else {
    c(0, Inf)
  }
}

bootstrap_skewness <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(NA_real_)
  m <- mean(x)
  s <- sd(x)
  if (is.na(s) || s == 0) return(NA_real_)
  mean((x - m)^3) / s^3
}


# ============================================================
# BCa functions
# ============================================================

bca_ci <- function(theta_hat, boot_values, jack_values, conf = 0.95) {
  
  alpha <- 1 - conf
  
  boot_values <- boot_values[!is.na(boot_values)]
  jack_values <- jack_values[!is.na(jack_values)]
  
  B <- length(boot_values)
  
  if (B < 10 || length(jack_values) < 3) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  
  p_hat <- mean(boot_values < theta_hat)
  p_hat <- safe_prob(p_hat, B)
  z0 <- qnorm(p_hat)
  
  jack_mean <- mean(jack_values)
  num <- sum((jack_mean - jack_values)^3)
  den <- 6 * (sum((jack_mean - jack_values)^2))^(3 / 2)
  
  if (is.na(den) || den == 0) {
    a <- 0
  } else {
    a <- num / den
  }
  
  z_low <- qnorm(alpha / 2)
  z_high <- qnorm(1 - alpha / 2)
  
  alpha1 <- pnorm(z0 + (z0 + z_low) / (1 - a * (z0 + z_low)))
  alpha2 <- pnorm(z0 + (z0 + z_high) / (1 - a * (z0 + z_high)))
  
  alpha1 <- safe_prob(alpha1, B)
  alpha2 <- safe_prob(alpha2, B)
  
  probs <- sort(c(alpha1, alpha2))
  
  ci <- quantile(
    boot_values,
    probs = probs,
    na.rm = TRUE,
    names = FALSE,
    type = 6
  )
  
  c(lower = ci[1], upper = ci[2])
}


compute_area_jackknife <- function(
    fixed,
    pop_data,
    smp_data,
    pop_domains,
    smp_domains,
    transformation,
    L,
    na_rm,
    indicator
) {
  
  areas <- unique(smp_data[[smp_domains]])
  
  jack_list <- vector("list", length(areas))
  names(jack_list) <- areas
  
  for (a in areas) {
    
    smp_jack <- smp_data[smp_data[[smp_domains]] != a, ]
    
    fit_jack <- ebp2(
      fixed = fixed,
      pop_data = pop_data,
      pop_domains = pop_domains,
      smp_data = smp_jack,
      smp_domains = smp_domains,
      na.rm = na_rm,
      transformation = transformation,
      L = L,
      MSE = FALSE
    )
    
    jack_list[[a]] <- data.frame(
      Domain = fit_jack$ind$Domain,
      value = fit_jack$ind[[indicator]]
    )
  }
  
  domains <- jack_list[[1]]$Domain
  
  jack_mat <- sapply(jack_list, function(x) {
    x$value[match(domains, x$Domain)]
  })
  
  rownames(jack_mat) <- domains
  
  jack_mat
}


# ============================================================
# CI computation
# ============================================================

compute_ci <- function(
    model,
    indicator,
    ci_method,
    fixed = NULL,
    pop_data = NULL,
    smp_data = NULL,
    pop_domains = NULL,
    smp_domains = NULL,
    transformation = NULL,
    L = NULL,
    na_rm = TRUE
) {
  
  estimates <- model$ind[[indicator]]
  domains <- model$ind$Domain
  mse <- model$MSE[[indicator]]
  se <- sqrt(mse)
  
  out <- data.frame(
    Domain = domains,
    Estimate = estimates,
    MSE = mse,
    SE = se
  )
  
  if (ci_method == "Normal") {
    
    out <- out |>
      mutate(
        lower = pmax(Estimate - 1.96 * SE, 0),
        upper = Estimate + 1.96 * SE
      )
    
  } else if (ci_method == "Percentile") {
    
    boot_values <- model$boot_estimates[, indicator, ]
    
    ci <- t(apply(
      boot_values,
      1,
      quantile,
      probs = c(0.025, 0.975),
      na.rm = TRUE
    ))
    
    out <- out |>
      mutate(
        lower = ci[, 1],
        upper = ci[, 2]
      )
    
  } else if (ci_method == "Reverse Percentile") {
    
    boot_values <- model$boot_estimates[, indicator, ]
    
    ci_perc <- t(apply(
      boot_values,
      1,
      quantile,
      probs = c(0.025, 0.975),
      na.rm = TRUE
    ))
    
    out <- out |>
      mutate(
        lower = pmax(2 * Estimate - ci_perc[, 2], 0),
        upper = 2 * Estimate - ci_perc[, 1]
      )
    
  } else if (ci_method == "BCa area") {
    
    boot_values <- model$boot_estimates[, indicator, ]
    
    jack_mat <- compute_area_jackknife(
      fixed = fixed,
      pop_data = pop_data,
      smp_data = smp_data,
      pop_domains = pop_domains,
      smp_domains = smp_domains,
      transformation = transformation,
      L = L,
      na_rm = na_rm,
      indicator = indicator
    )
    
    ci <- matrix(NA_real_, nrow = length(domains), ncol = 2)
    colnames(ci) <- c("lower", "upper")
    
    for (i in seq_along(domains)) {
      
      domain_i <- domains[i]
      
      boot_i <- boot_values[i, ]
      jack_i <- jack_mat[domain_i, ]
      
      ci[i, ] <- bca_ci(
        theta_hat = estimates[i],
        boot_values = boot_i,
        jack_values = jack_i
      )
    }
    
    out <- out |>
      mutate(
        lower = pmax(ci[, "lower"], 0),
        upper = ci[, "upper"]
      )
  }
  
  out |>
    mutate(
      CI_length = upper - lower,
      Relative_length = if_else(Estimate != 0, CI_length / abs(Estimate), NA_real_),
      Contains_estimate = Estimate >= lower & Estimate <= upper
    )
}


compute_compare_ci <- function(
    model,
    indicator,
    methods,
    fixed,
    pop_data,
    smp_data,
    pop_domains,
    smp_domains,
    transformation,
    L,
    na_rm
) {
  
  base <- data.frame(
    Domain = model$ind$Domain,
    Estimate = model$ind[[indicator]]
  )
  
  method_tables <- lapply(methods, function(m) {
    
    tmp <- compute_ci(
      model = model,
      indicator = indicator,
      ci_method = m,
      fixed = fixed,
      pop_data = pop_data,
      smp_data = smp_data,
      pop_domains = pop_domains,
      smp_domains = smp_domains,
      transformation = transformation,
      L = L,
      na_rm = na_rm
    )
    
    suf <- method_suffix(m)
    
    tmp2 <- tmp |>
      select(Domain, lower, upper, CI_length, Relative_length, Contains_estimate)
    
    names(tmp2) <- c(
      "Domain",
      paste0("lower_", suf),
      paste0("upper_", suf),
      paste0("length_", suf),
      paste0("rel_length_", suf),
      paste0("contains_estimate_", suf)
    )
    
    tmp2
  })
  
  Reduce(function(x, y) left_join(x, y, by = "Domain"), c(list(base), method_tables))
}


# ============================================================
# UI
# ============================================================

ui <- page_sidebar(
  
  title = "Bootstrap Confidence Intervals for EBP Small Area Estimates",
  
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2C3E50",
    secondary = "#18BC9C",
    base_font = font_google("Inter"),
    heading_font = font_google("Inter")
  ),
  
  tags$head(
    tags$style(HTML("
      .nav-tabs {
        display: flex;
        flex-wrap: wrap;
        width: 100%;
      }

      .nav-tabs .nav-item {
        flex: 1 1 150px;
        text-align: center;
      }

      .nav-tabs .nav-link {
        width: 100%;
        text-align: center;
      }
    "))
  ),
  
  sidebar = sidebar(
    width = 400,
    
    h4(bsicons::bs_icon("database"), " Input data"),
    
    checkboxInput(
      "use_example",
      "Use example data included in the app",
      value = FALSE
    ),
    
    fileInput(
      "pop_file",
      "Upload population data",
      accept = c(".csv", ".rds")
    ),
    
    fileInput(
      "smp_file",
      "Upload sample data",
      accept = c(".csv", ".rds")
    ),
    
    tags$small(
      style = "color: #6c757d;",
      "If 'Use example data' is selected, the app reads pop.rds and smp.rds from the app folder."
    ),
    
    hr(),
    
    h4(bsicons::bs_icon("sliders"), " Model specification"),
    
    textInput(
      "formula_text",
      "Model formula",
      value = "eqIncome ~ gender + eqsize + cash + self_empl + unempl_ben + age_ben + surv_ben + sick_ben + dis_ben + rent + fam_allow + house_allow + cap_inv + tax_adj"
    ),
    
    textInput(
      "pop_domains",
      "Population domain column",
      value = "district"
    ),
    
    textInput(
      "smp_domains",
      "Sample domain column",
      value = "district"
    ),
    
    selectInput(
      "transformation",
      "Transformation",
      choices = c(
        "No transformation" = "no",
        "Log" = "log",
        "Box-Cox" = "box.cox",
        "Dual" = "dual",
        "Log-shift" = "log.shift",
        "Logit" = "logit",
        "Probit" = "probit",
        "Arcsin" = "arcsin"
      ),
      selected = "no"
    ),
    
    selectInput(
      "indicator",
      "Indicator",
      choices = c("Mean", "Head_Count", "Gini"),
      selected = "Head_Count"
    ),
    
    selectInput(
      "boot_type",
      "Bootstrap type",
      choices = c("parametric", "wild"),
      selected = "parametric"
    ),
    
    numericInput(
      "B",
      "Number of bootstrap replicates",
      value = 100,
      min = 5,
      step = 50
    ),
    
    numericInput(
      "L",
      "Number of Monte Carlo simulations (L)",
      value = 50,
      min = 1,
      step = 10
    ),
    
    selectInput(
      "ci_method",
      "Confidence interval method",
      choices = c(
        "Normal",
        "Percentile",
        "Reverse Percentile",
        "BCa area"
      ),
      selected = "Percentile"
    ),
    
    checkboxGroupInput(
      "compare_methods",
      "Methods to compare",
      choices = c(
        "Normal",
        "Percentile",
        "Reverse Percentile",
        "BCa area"
      ),
      selected = c("Normal", "Percentile")
    ),
    
    checkboxInput(
      "na_rm",
      "Remove missing values",
      value = TRUE
    ),
    
    hr(),
    
    h4(bsicons::bs_icon("bar-chart"), " Plot options"),
    
    numericInput(
      "n_plot",
      "Number of domains shown in the plot",
      value = 30,
      min = 5,
      step = 5
    ),
    
    selectInput(
      "plot_order",
      "Domains to show",
      choices = c(
        "Highest estimates",
        "Lowest estimates",
        "All domains"
      ),
      selected = "Highest estimates"
    ),
    
    selectInput(
      "diagnostic_domain",
      "Domain for bootstrap diagnostics",
      choices = NULL
    ),
    
    div(
      class = "d-grid gap-2 mt-3",
      actionButton(
        "run_model",
        "Run model",
        class = "btn-primary btn-lg"
      ),
      downloadButton(
        "download_results",
        "Download CI table",
        class = "btn-outline-secondary"
      ),
      downloadButton(
        "download_report",
        "Download HTML report",
        class = "btn-outline-secondary"
      )
    ),
    
    hr(),
    
    tags$small(
      style = "color: #6c757d;",
      "Tip: BCa area may take longer because the model is re-fitted once for each area."
    )
  ),
  
  navset_card_tab(
    
    title = div(
      style = "font-weight: 600;",
      bsicons::bs_icon("window"),
      " Outputs"
    ),
    
    nav_panel(
      title = "Home",
      br(),
      
      card(
        card_body(
          h2("Bootstrap Confidence Intervals for EBP Small Area Estimates"),
          
          p(
            "This Shiny app computes Empirical Best Prediction small area estimates ",
            "using modified emdi functions and constructs bootstrap-based confidence intervals."
          ),
          
          p(
            "The app allows the user to load population and sample data, specify an EBP model, ",
            "choose the transformation, compute bootstrap-based confidence intervals, ",
            "compare interval methods and inspect bootstrap diagnostics."
          ),
          
          tags$hr(),
          
          h4("Available tools"),
          
          tags$ul(
            tags$li("Upload population and sample datasets."),
            tags$li("Run unit-level EBP models using modified emdi functions."),
            tags$li("Use standard and additional transformations: logit, probit and arcsin."),
            tags$li("Compute Normal, Percentile, Reverse Percentile and BCa area confidence intervals."),
            tags$li("Compare confidence interval methods across domains."),
            tags$li("Inspect bootstrap distributions and diagnostic summaries."),
            tags$li("Download results and HTML reports.")
          ),
          
          tags$hr(),
          
          h4("References"),
          
          tags$ul(
            tags$li(
              tags$a(
                href = "https://www.tandfonline.com/doi/pdf/10.1080/00949650601141811",
                target = "_blank",
                " González-Manteiga et al. (2008)."
              )
            ),
            tags$li(
              tags$a(
                href = "https://www.jstatsoft.org/article/view/v091i07/0",
                target = "_blank",
                "Kreutzmann et al. (2019)."
              )
            )
          ),
          
          p(
            "To test the application without uploading external data, click the button below ",
            "to use the example population and sample datasets included in the app."
          ),
          
          actionButton(
            "load_example_home",
            "Load example dataset",
            class = "btn-primary btn-lg"
          )
        )
      )
    ),
    
    nav_panel(
      title = "Input summary",
      br(),
      uiOutput("input_warnings"),
      DTOutput("input_summary_table")
    ),
    
    nav_panel(
      title = "Results table",
      br(),
      uiOutput("quality_warnings"),
      DTOutput("ci_table")
    ),
    
    nav_panel(
      title = "Model diagnostic",
      br(),
      
      layout_column_wrap(
        width = 1,
        
        card(
          card_header("Checks"),
          card_body(
            DTOutput("checks_table")
          )
        ),
        
        card(
          card_header("Model status"),
          card_body(
            verbatimTextOutput("status")
          )
        )
      )
    ), 

    nav_panel(
      title = "Bootstrap diagnostics",
      br(),
      
      uiOutput("diagnostic_warning"),
      
      layout_column_wrap(
        width = 1,
        
        card(
          card_header("Bootstrap summary"),
          card_body(
            DTOutput("bootstrap_diag_table", height = "auto")
          )
        ),
        
        card(
          card_header("Bootstrap distribution"),
          card_body(
            plotOutput("bootstrap_hist", height = "500px")
          )
        )
      )
    ),
    
    nav_panel(
      title = "CI plot",
      br(),
      uiOutput("plot_warning"),
      plotOutput("ci_plot", height = "750px")
    ),
    
    nav_panel(
      title = "Compare CI methods",
      br(),
      uiOutput("compare_warning"),
      DTOutput("compare_table"),
      br()#,
      #plotOutput("compare_plot", height = "750px")
    ),
    
    nav_panel(
      title = "About",
      br(),
      
      h4("About this app"),
      p(
        "This Shiny app computes Empirical Best Prediction small area estimates using modified emdi functions and constructs bootstrap-based confidence intervals."
      ),
      p(
        "Available confidence interval methods currently include Normal, Percentile, Reverse Percentile and BCa area intervals."
      ),
      p(
        "The app includes tools for comparing confidence interval methods, inspecting bootstrap distributions, checking interval quality and downloading reproducible outputs."
      ),
      
      tags$hr(),
      
      h4("Confidence interval methods"),
      
      p(strong("Normal interval")),
      tags$p("The Normal interval is computed as: estimate ± 1.96 × sqrt(MSE)."),
      
      p(strong("Percentile interval")),
      tags$p("The Percentile interval is obtained from the 2.5th and 97.5th empirical quantiles of the bootstrap distribution."),
      
      p(strong("Reverse Percentile interval")),
      tags$p("The Reverse Percentile interval reflects the Percentile bounds around the original point estimate."),
      
      p(strong("BCa area interval")),
      tags$p("The BCa area interval adjusts the percentile bounds for bias and acceleration. The acceleration term is computed using an area-level jackknife."),
      
      tags$hr(),
      
      h4("References"),
      tags$ul(
        tags$li(
          tags$a(
            href = "https://www.tandfonline.com/doi/pdf/10.1080/00949650601141811",
            target = "_blank",
            "González-Manteiga et al. (2008)."
          )
        ),
        tags$li(
          tags$a(
            href = "https://www.jstatsoft.org/article/view/v091i07/0",
            target = "_blank",
            "Kreutzmann et al. (2019)."
          )
        )
      ),
      
      tags$hr(),
      
      h4("Contact"),
      p(
        "For questions, bug reports or problems with the app, please contact: ",
        tags$a(
          href = "mailto:lorenzo.mori7@unibo.it",
          "lorenzo.mori7@unibo.it"
        ),
        "."
      ),
      
      tags$hr(),
      
      tags$small(
        style = "color: #6c757d;",
        "Recommended for quick online testing: use small values of B and L first, then increase them for final analyses."
      )
    )
  )
)


# ============================================================
# Server
# ============================================================

server <- function(input, output, session) {
  
  pop_data <- reactive({
    if (isTRUE(input$use_example)) {
      validate(
        need(file.exists("pop.rds"), "Example file pop.rds was not found in the app folder.")
      )
      readRDS("pop.rds")
    } else {
      req(input$pop_file)
      read_uploaded_data(input$pop_file)
    }
  })
  
  smp_data <- reactive({
    if (isTRUE(input$use_example)) {
      validate(
        need(file.exists("smp.rds"), "Example file smp.rds was not found in the app folder.")
      )
      readRDS("smp.rds")
    } else {
      req(input$smp_file)
      read_uploaded_data(input$smp_file)
    }
  })
  
  observeEvent(input$load_example_home, {
    
    updateCheckboxInput(
      session,
      "use_example",
      value = TRUE
    )
    
    showNotification(
      "Example data loaded. You can now move to the Input summary tab or run the model.",
      type = "message"
    )
  })
  
  
  input_summary <- reactive({
    
    req(pop_data())
    req(smp_data())
    
    pop <- pop_data()
    smp <- smp_data()
    
    formula_vars <- tryCatch(all.vars(as.formula(input$formula_text)), error = function(e) character(0))
    outcome <- if (length(formula_vars) > 0) formula_vars[1] else NA_character_
    
    pop_domain_ok <- input$pop_domains %in% names(pop)
    smp_domain_ok <- input$smp_domains %in% names(smp)
    
    n_pop_domains <- if (pop_domain_ok) length(unique(pop[[input$pop_domains]])) else NA_integer_
    n_smp_domains <- if (smp_domain_ok) length(unique(smp[[input$smp_domains]])) else NA_integer_
    
    sample_sizes <- if (smp_domain_ok) table(smp[[input$smp_domains]]) else NA
    
    data.frame(
      Item = c(
        "Population rows",
        "Sample rows",
        "Population columns",
        "Sample columns",
        "Population domain column",
        "Sample domain column",
        "Number of population domains",
        "Number of sample domains",
        "Minimum sample size per domain",
        "Maximum sample size per domain",
        "Outcome variable",
        "Formula variables missing in population data",
        "Formula variables missing in sample data"
      ),
      Value = c(
        nrow(pop),
        nrow(smp),
        ncol(pop),
        ncol(smp),
        input$pop_domains,
        input$smp_domains,
        n_pop_domains,
        n_smp_domains,
        if (smp_domain_ok) min(sample_sizes) else NA,
        if (smp_domain_ok) max(sample_sizes) else NA,
        outcome,
        paste(setdiff(formula_vars, names(pop)), collapse = ", "),
        paste(setdiff(formula_vars, names(smp)), collapse = ", ")
      )
    )
  })
  
  output$input_summary_table <- renderDT({
    req(input_summary())
    
    datatable(
      input_summary(),
      rownames = FALSE,
      options = list(
        pageLength = 15,
        dom = "t",
        scrollX = TRUE
      )
    )
  })
  
  output$input_warnings <- renderUI({
    
    req(pop_data())
    req(smp_data())
    
    pop <- pop_data()
    smp <- smp_data()
    
    warnings <- character(0)
    
    if (!(input$pop_domains %in% names(pop))) {
      warnings <- c(warnings, paste0("Population domain column '", input$pop_domains, "' was not found."))
    }
    
    if (!(input$smp_domains %in% names(smp))) {
      warnings <- c(warnings, paste0("Sample domain column '", input$smp_domains, "' was not found."))
    }
    
    formula_vars <- tryCatch(all.vars(as.formula(input$formula_text)), error = function(e) character(0))
    
    missing_pop <- setdiff(formula_vars, names(pop))
    missing_smp <- setdiff(formula_vars, names(smp))
    
    if (length(missing_pop) > 0) {
      warnings <- c(warnings, paste0("Variables missing in population data: ", paste(missing_pop, collapse = ", ")))
    }
    
    if (length(missing_smp) > 0) {
      warnings <- c(warnings, paste0("Variables missing in sample data: ", paste(missing_smp, collapse = ", ")))
    }
    
    if (input$pop_domains %in% names(pop) && input$smp_domains %in% names(smp)) {
      
      pop_domains <- unique(pop[[input$pop_domains]])
      smp_domains <- unique(smp[[input$smp_domains]])
      
      missing_in_sample <- setdiff(pop_domains, smp_domains)
      extra_in_sample <- setdiff(smp_domains, pop_domains)
      
      if (length(missing_in_sample) > 0) {
        warnings <- c(
          warnings,
          paste0(length(missing_in_sample), " population domains are not present in the sample.")
        )
      }
      
      if (length(extra_in_sample) > 0) {
        warnings <- c(
          warnings,
          paste0(length(extra_in_sample), " sample domains are not present in the population.")
        )
      }
    }
    
    if (length(warnings) == 0) {
      div(
        class = "alert alert-success",
        strong("Input check passed: "),
        "No obvious input problems were detected."
      )
    } else {
      tagList(lapply(warnings, function(w) {
        div(class = "alert alert-warning", strong("Warning: "), w)
      }))
    }
  })
  
  model_result <- eventReactive(input$run_model, {
    
    req(pop_data())
    req(smp_data())
    
    fixed_formula <- as.formula(input$formula_text)
    
    transformation_error <- check_transformation_compatibility(
      smp_data = smp_data(),
      fixed_formula = fixed_formula,
      transformation = input$transformation
    )
    
    if (!is.null(transformation_error)) {
      showNotification(
        transformation_error,
        type = "error",
        duration = 10
      )
      
      validate(
        need(FALSE, transformation_error)
      )
    }
    
    withProgress(message = "Running EBP model...", value = 0, {
      
      incProgress(0.2, detail = "Preparing data")
      
      fit <- tryCatch(
        {
          ebp2(
            fixed = fixed_formula,
            pop_data = pop_data(),
            pop_domains = input$pop_domains,
            smp_data = smp_data(),
            smp_domains = input$smp_domains,
            na.rm = input$na_rm,
            transformation = input$transformation,
            L = input$L,
            MSE = TRUE,
            boot_type = input$boot_type,
            B = input$B
          )
        },
        error = function(e) {
          showNotification(
            paste("Model estimation failed:", conditionMessage(e)),
            type = "error",
            duration = 12
          )
          
          validate(
            need(FALSE, paste("Model estimation failed:", conditionMessage(e)))
          )
        }
      )
      
      incProgress(1, detail = "Done")
      
      fit
    })
  })  
  observeEvent(model_result(), {
    updateSelectInput(
      session,
      "diagnostic_domain",
      choices = model_result()$ind$Domain,
      selected = model_result()$ind$Domain[1]
    )
  })
  
  ci_result <- reactive({
    
    req(model_result())
    
    fixed_formula <- as.formula(input$formula_text)
    
    withProgress(message = "Computing confidence intervals...", value = 0, {
      
      incProgress(0.3, detail = paste("Method:", input$ci_method))
      
      ci <- compute_ci(
        model = model_result(),
        indicator = input$indicator,
        ci_method = input$ci_method,
        fixed = fixed_formula,
        pop_data = pop_data(),
        smp_data = smp_data(),
        pop_domains = input$pop_domains,
        smp_domains = input$smp_domains,
        transformation = input$transformation,
        L = input$L,
        na_rm = input$na_rm
      )
      
      incProgress(1, detail = "Done")
      
      ci
    })
  })
  
  compare_result <- reactive({
    
    req(model_result())
    req(input$compare_methods)
    
    fixed_formula <- as.formula(input$formula_text)
    
    withProgress(message = "Computing method comparison...", value = 0, {
      
      incProgress(0.2, detail = "Computing selected methods")
      
      out <- compute_compare_ci(
        model = model_result(),
        indicator = input$indicator,
        methods = input$compare_methods,
        fixed = fixed_formula,
        pop_data = pop_data(),
        smp_data = smp_data(),
        pop_domains = input$pop_domains,
        smp_domains = input$smp_domains,
        transformation = input$transformation,
        L = input$L,
        na_rm = input$na_rm
      )
      
      incProgress(1, detail = "Done")
      
      out
    })
  })
  
  output$quality_warnings <- renderUI({
    
    req(ci_result())
    
    df <- ci_result()
    support <- indicator_support(input$indicator)
    
    warnings <- character(0)
    
    if (input$B < 200) {
      warnings <- c(warnings, "The number of bootstrap replicates is below 200. Results should be interpreted with caution.")
    }
    
    n_not_contains <- sum(!df$Contains_estimate, na.rm = TRUE)
    if (n_not_contains > 0) {
      warnings <- c(
        warnings,
        paste0(n_not_contains, " confidence intervals do not contain the point estimate.")
      )
    }
    
    if (any(is.na(df$lower) | is.na(df$upper))) {
      warnings <- c(warnings, "Some confidence interval bounds are NA.")
    }
    
    if (any(df$lower > df$upper, na.rm = TRUE)) {
      warnings <- c(warnings, "Some intervals have lower bound greater than upper bound.")
    }
    
    if (is.finite(support[2]) && any(df$upper > support[2], na.rm = TRUE)) {
      warnings <- c(
        warnings,
        paste0("Some upper bounds exceed the natural support of ", input$indicator, ".")
      )
    }
    
    if (any(df$lower < support[1], na.rm = TRUE)) {
      warnings <- c(
        warnings,
        paste0("Some lower bounds are below the natural support of ", input$indicator, ".")
      )
    }
    
    if (length(warnings) == 0) {
      div(
        class = "alert alert-success",
        strong("CI check passed: "),
        "No obvious interval problems were detected."
      )
    } else {
      tagList(lapply(warnings, function(w) {
        div(class = "alert alert-warning", strong("Warning: "), w)
      }))
    }
  })
  
  output$plot_warning <- renderUI({
    
    req(ci_result())
    
    n_domains <- nrow(ci_result())
    
    if (input$plot_order != "All domains" && input$n_plot > n_domains) {
      div(
        class = "alert alert-warning",
        strong("Warning: "),
        paste0(
          "You requested ", input$n_plot,
          " domains, but only ", n_domains,
          " domains are available. The plot will show all available domains."
        )
      )
    }
  })
  
  output$compare_warning <- renderUI({
    
    if ("BCa area" %in% input$compare_methods) {
      div(
        class = "alert alert-warning",
        strong("Warning: "),
        "BCa area is selected in the comparison. This may be slow because the model is re-fitted once for each area."
      )
    }
  })
  
  plot_data <- reactive({
    
    req(ci_result())
    
    df <- ci_result() |>
      arrange(desc(Estimate))
    
    n_domains <- nrow(df)
    n_to_plot <- min(input$n_plot, n_domains)
    
    if (input$plot_order == "Highest estimates") {
      df <- df |>
        slice_head(n = n_to_plot)
    } else if (input$plot_order == "Lowest estimates") {
      df <- df |>
        arrange(Estimate) |>
        slice_head(n = n_to_plot)
    } else if (input$plot_order == "All domains") {
      df <- df
    }
    
    df |>
      arrange(Estimate) |>
      mutate(
        Domain = factor(Domain, levels = Domain)
      )
  })
  
  compare_plot_data <- reactive({
    
    req(model_result())
    req(input$compare_methods)
    
    fixed_formula <- as.formula(input$formula_text)
    
    all_ci <- lapply(input$compare_methods, function(m) {
      tmp <- compute_ci(
        model = model_result(),
        indicator = input$indicator,
        ci_method = m,
        fixed = fixed_formula,
        pop_data = pop_data(),
        smp_data = smp_data(),
        pop_domains = input$pop_domains,
        smp_domains = input$smp_domains,
        transformation = input$transformation,
        L = input$L,
        na_rm = input$na_rm
      )
      tmp$Method <- m
      tmp
    })
    
    df <- bind_rows(all_ci) |>
      arrange(desc(Estimate))
    
    n_domains <- length(unique(df$Domain))
    n_to_plot <- min(input$n_plot, n_domains)
    
    selected_domains <- df |>
      distinct(Domain, Estimate) |>
      arrange(desc(Estimate))
    
    if (input$plot_order == "Highest estimates") {
      selected_domains <- selected_domains |>
        slice_head(n = n_to_plot)
    } else if (input$plot_order == "Lowest estimates") {
      selected_domains <- selected_domains |>
        arrange(Estimate) |>
        slice_head(n = n_to_plot)
    }
    
    if (input$plot_order != "All domains") {
      df <- df |>
        filter(Domain %in% selected_domains$Domain)
    }
    
    df |>
      arrange(Estimate) |>
      mutate(
        Domain = factor(Domain, levels = unique(Domain))
      )
  })
  
  bootstrap_diag <- reactive({
    
    req(model_result())
    req(input$diagnostic_domain)
    req(input$indicator)
    
    domains <- model_result()$ind$Domain
    idx <- match(input$diagnostic_domain, domains)
    
    validate(
      need(!is.na(idx), "Selected domain not found.")
    )
    
    validate(
      need(input$indicator %in% dimnames(model_result()$boot_estimates)[[2]],
           paste0("Indicator ", input$indicator, " not found in bootstrap estimates."))
    )
    
    boot_values <- model_result()$boot_estimates[idx, input$indicator, , drop = TRUE]
    boot_values <- as.numeric(boot_values)
    boot_values <- boot_values[is.finite(boot_values)]
    
    theta_hat <- model_result()$ind[[input$indicator]][idx]
    
    validate(
      need(length(boot_values) > 0, "No valid bootstrap values are available for this domain and indicator.")
    )
    
    data.frame(
      Statistic = c(
        "Domain",
        "Indicator",
        "Point estimate",
        "Bootstrap mean",
        "Bootstrap median",
        "Bootstrap SD",
        "Bootstrap skewness",
        "2.5% quantile",
        "97.5% quantile",
        "Minimum bootstrap value",
        "Maximum bootstrap value",
        "Number of valid bootstrap replicates"
      ),
      Value = c(
        input$diagnostic_domain,
        input$indicator,
        theta_hat,
        mean(boot_values, na.rm = TRUE),
        median(boot_values, na.rm = TRUE),
        sd(boot_values, na.rm = TRUE),
        bootstrap_skewness(boot_values),
        quantile(boot_values, 0.025, na.rm = TRUE, names = FALSE),
        quantile(boot_values, 0.975, na.rm = TRUE, names = FALSE),
        min(boot_values, na.rm = TRUE),
        max(boot_values, na.rm = TRUE),
        length(boot_values)
      )
    )
  })  
  checks_result <- reactive({
    
    req(ci_result())
    
    df <- ci_result()
    support <- indicator_support(input$indicator)
    
    data.frame(
      Check = c(
        "Number of domains",
        "Number of bootstrap replicates",
        "Number of intervals with NA bounds",
        "Number of intervals not containing the point estimate",
        "Number of intervals with lower > upper",
        "Number of intervals with lower bound equal to 0",
        "Number of intervals below natural support",
        "Number of intervals above natural support",
        "Average CI length",
        "Median CI length",
        "Maximum CI length",
        "Average relative CI length",
        "Median relative CI length"
      ),
      Value = c(
        nrow(df),
        input$B,
        sum(is.na(df$lower) | is.na(df$upper)),
        sum(!df$Contains_estimate, na.rm = TRUE),
        sum(df$lower > df$upper, na.rm = TRUE),
        sum(df$lower == 0, na.rm = TRUE),
        sum(df$lower < support[1], na.rm = TRUE),
        if (is.finite(support[2])) sum(df$upper > support[2], na.rm = TRUE) else NA,
        mean(df$CI_length, na.rm = TRUE),
        median(df$CI_length, na.rm = TRUE),
        max(df$CI_length, na.rm = TRUE),
        mean(df$Relative_length, na.rm = TRUE),
        median(df$Relative_length, na.rm = TRUE)
      )
    )
  })
  
  output$ci_table <- renderDT({
    req(ci_result())
    
    datatable(
      ci_result(),
      options = list(
        pageLength = 15,
        scrollX = TRUE
      )
    ) |>
      formatRound(
        columns = c(
          "Estimate", "MSE", "SE", "lower", "upper",
          "CI_length", "Relative_length"
        ),
        digits = 4
      )
  })
  
  output$ci_plot <- renderPlot({
    
    req(plot_data())
    
    df <- plot_data()
    
    ggplot(df, aes(y = Domain)) +
      geom_segment(
        aes(x = lower, xend = upper, yend = Domain),
        linewidth = 0.8,
        alpha = 0.7
      ) +
      geom_point(
        aes(x = Estimate),
        size = 2.4
      ) +
      labs(
        title = paste("Confidence intervals -", input$indicator),
        subtitle = paste("CI method:", input$ci_method),
        x = "Estimate and confidence interval",
        y = NULL
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey30"),
        panel.grid.minor = element_blank()
      )
  })
  
  output$compare_table <- renderDT({
    req(compare_result())
    
    datatable(
      compare_result(),
      options = list(
        pageLength = 15,
        scrollX = TRUE
      )
    ) |>
      formatRound(
        columns = names(compare_result())[sapply(compare_result(), is.numeric)],
        digits = 4
      )
  })
  
  output$compare_plot <- renderPlot({
    
    req(compare_plot_data())
    
    df <- compare_plot_data()
    
    ggplot(df, aes(y = Domain)) +
      geom_segment(
        aes(x = lower, xend = upper, yend = Domain),
        linewidth = 0.75,
        alpha = 0.7
      ) +
      geom_point(
        aes(x = Estimate),
        size = 2
      ) +
      facet_wrap(~ Method, scales = "free_x") +
      labs(
        title = paste("Comparison of confidence interval methods -", input$indicator),
        x = "Estimate and confidence interval",
        y = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
  })
  
  output$diagnostic_warning <- renderUI({
    
    req(input$diagnostic_domain)
    
    div(
      class = "alert alert-info",
      strong("Bootstrap diagnostics: "),
      paste0(
        "The histogram shows the bootstrap distribution for domain ",
        input$diagnostic_domain,
        " and indicator ",
        input$indicator,
        "."
      )
    )
  })
  
  output$bootstrap_diag_table <- renderDT({
    req(bootstrap_diag())
    
    datatable(
      bootstrap_diag(),
      rownames = FALSE,
      options = list(
        pageLength = 15,
        scrollX = TRUE
      )
    )
  })  
  output$bootstrap_hist <- renderPlot({
    
    req(model_result())
    req(input$diagnostic_domain)
    req(input$indicator)
    
    domains <- model_result()$ind$Domain
    idx <- match(input$diagnostic_domain, domains)
    
    validate(
      need(!is.na(idx), "Selected domain not found.")
    )
    
    validate(
      need(input$indicator %in% dimnames(model_result()$boot_estimates)[[2]],
           paste0("Indicator ", input$indicator, " not found in bootstrap estimates."))
    )
    
    boot_values <- model_result()$boot_estimates[idx, input$indicator, , drop = TRUE]
    boot_values <- as.numeric(boot_values)
    boot_values <- boot_values[is.finite(boot_values)]
    
    validate(
      need(length(boot_values) > 0, "No valid bootstrap values are available for this domain and indicator.")
    )
    
    theta_hat <- model_result()$ind[[input$indicator]][idx]
    
    df <- data.frame(value = boot_values)
    
    ggplot(df, aes(x = value)) +
      geom_histogram(bins = min(30, max(5, length(unique(boot_values)))), alpha = 0.8) +
      geom_vline(xintercept = theta_hat, linewidth = 1) +
      labs(
        title = paste("Bootstrap distribution -", input$diagnostic_domain),
        subtitle = paste("Vertical line: original point estimate | Indicator:", input$indicator),
        x = "Bootstrap estimates",
        y = "Frequency"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
  })  
  output$checks_table <- renderDT({
    req(checks_result())
    
    datatable(
      checks_result(),
      rownames = FALSE,
      options = list(
        pageLength = 20,
        dom = "t",
        scrollX = TRUE
      )
    )
  })
  
  output$status <- renderPrint({
    req(model_result())
    
    
    cat("Model successfully estimated.\n")
    cat("Selected indicator:", input$indicator, "\n")
    cat("Selected CI method:", input$ci_method, "\n")
    cat("Bootstrap replicates:", input$B, "\n")
    cat("Monte Carlo simulations L:", input$L, "\n")
    cat("Transformation:", input$transformation, "\n")
    tp <- model_result()$transform_param
    
    if (!is.null(tp)) {
      
      if (input$transformation %in% c("box.cox", "dual")) {
        
        if (!is.null(tp$optimal_lambda) && length(tp$optimal_lambda) > 0) {
          cat("Transformation lambda:", round(tp$optimal_lambda, 4), "\n")
        }
        
        if (!is.null(tp$shift_par) && length(tp$shift_par) > 0) {
          cat("Transformation shift:", round(tp$shift_par, 4), "\n")
        }
      }
      
      if (input$transformation == "log.shift") {
        
        if (!is.null(tp$optimal_lambda) && length(tp$optimal_lambda) > 0) {
          cat("Log-shift parameter:", round(tp$optimal_lambda, 4), "\n")
        }
      }
    }
    cat("Bootstrap type:", input$boot_type, "\n")
    
    if (isTRUE(input$use_example)) {
      cat("Data source: example files pop.rds and smp.rds.\n")
    } else {
      cat("Data source: uploaded files.\n")
    }
    
    if (input$ci_method == "BCa area") {
      cat("Note: BCa area refits the model once for each area, so it may take longer.\n")
    }
  })
  
  output$download_results <- downloadHandler(
    filename = function() {
      paste0(
        "ci_results_",
        input$indicator,
        "_",
        gsub(" ", "_", input$ci_method),
        ".csv"
      )
    },
    content = function(file) {
      write.csv(ci_result(), file, row.names = FALSE)
    }
  )
  
  output$download_report <- downloadHandler(
    filename = function() {
      paste0(
        "ci_report_",
        input$indicator,
        "_",
        gsub(" ", "_", input$ci_method),
        ".html"
      )
    },
    content = function(file) {
      
      ci <- ci_result()
      checks <- checks_result()
      summary_tbl <- input_summary()
      
      html_table <- function(df) {
        paste0(
          "<table border='1' cellspacing='0' cellpadding='5'>",
          "<tr>",
          paste0("<th>", names(df), "</th>", collapse = ""),
          "</tr>",
          paste(
            apply(df, 1, function(row) {
              paste0(
                "<tr>",
                paste0("<td>", row, "</td>", collapse = ""),
                "</tr>"
              )
            }),
            collapse = "\n"
          ),
          "</table>"
        )
      }
      
      report <- paste0(
        "<html><head><title>CI report</title>",
        "<style>",
        "body { font-family: Arial, sans-serif; margin: 30px; }",
        "h1, h2 { color: #2C3E50; }",
        "table { border-collapse: collapse; font-size: 12px; }",
        "th { background-color: #f2f2f2; }",
        "</style>",
        "</head><body>",
        "<h1>Bootstrap Confidence Interval Report</h1>",
        "<h2>Model settings</h2>",
        "<p><b>Formula:</b> ", input$formula_text, "</p>",
        "<p><b>Indicator:</b> ", input$indicator, "</p>",
        "<p><b>CI method:</b> ", input$ci_method, "</p>",
        "<p><b>Transformation:</b> ", input$transformation, "</p>",
        "<p><b>Bootstrap type:</b> ", input$boot_type, "</p>",
        "<p><b>B:</b> ", input$B, "</p>",
        "<p><b>L:</b> ", input$L, "</p>",
        "<h2>Input summary</h2>",
        html_table(summary_tbl),
        "<h2>Checks</h2>",
        html_table(checks),
        "<h2>Confidence interval results</h2>",
        html_table(ci),
        "</body></html>"
      )
      
      writeLines(report, file)
    }
  )
}


shinyApp(ui = ui, server = server)