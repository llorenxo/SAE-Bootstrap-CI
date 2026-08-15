# ============================================================
# Additional transformations for unit-level EBP
# logit, probit, arcsin
# ============================================================

allowed_transformations_ebp2 <- c(
  "no", "log", "box.cox", "dual", "log.shift",
  "logit", "probit", "arcsin"
)

check_unit_interval <- function(y, transformation) {
  y_num <- as.numeric(y)
  
  if (any(y_num < 0 | y_num > 1, na.rm = TRUE)) {
    stop(
      paste0(
        "Transformation '", transformation,
        "' requires the response variable to be in [0, 1]."
      ),
      call. = FALSE
    )
  }
  
  y
}

clip_open01 <- function(y, eps = 1e-8) {
  pmin(pmax(y, eps), 1 - eps)
}

clip_closed01 <- function(y) {
  pmin(pmax(y, 0), 1)
}


# ------------------------------------------------------------
# Logit transformation:
# log(y / (1 - y))
# ------------------------------------------------------------

logit_transform <- function(y, shift = NULL, eps = 1e-8) {
  y <- check_unit_interval(y, transformation = "logit")
  y <- clip_open01(y, eps = eps)
  
  yt <- log(y / (1 - y))
  
  return(list(y = yt, shift = NULL))
}

logit_transform_back <- function(y) {
  plogis(y)
}


# ------------------------------------------------------------
# Probit transformation:
# Phi^{-1}(y)
# ------------------------------------------------------------

probit_transform <- function(y, shift = NULL, eps = 1e-8) {
  y <- check_unit_interval(y, transformation = "probit")
  y <- clip_open01(y, eps = eps)
  
  yt <- qnorm(y)
  
  return(list(y = yt, shift = NULL))
}

probit_transform_back <- function(y) {
  pnorm(y)
}


# ------------------------------------------------------------
# Arcsin transformation:
# arcsin(sqrt(y))
# ------------------------------------------------------------

arcsin_transform <- function(y, shift = NULL) {
  y <- check_unit_interval(y, transformation = "arcsin")
  y <- clip_closed01(y)
  
  yt <- asin(sqrt(y))
  
  return(list(y = yt, shift = NULL))
}

arcsin_transform_back <- function(y) {
  y <- pmin(pmax(y, 0), pi / 2)
  sin(y)^2
}

# ============================================================
# Modified emdi transformation functions
# ============================================================

data_transformation <- function(fixed, smp_data, transformation, lambda) {
  
  y_vector <- as.matrix(smp_data[paste(fixed[2])])
  
  transformed <- if (transformation == "no") {
    no_transform(y = y_vector, shift = NULL)
    
  } else if (transformation == "log") {
    log_transform(y = y_vector, shift = 0)
    
  } else if (transformation == "box.cox") {
    box_cox(y = y_vector, lambda = lambda, shift = 0)
    
  } else if (transformation == "dual") {
    dual(y = y_vector, lambda = lambda, shift = 0)
    
  } else if (transformation == "log.shift") {
    log_shift_opt(y = y_vector, lambda = lambda, shift = NULL)
    
  } else if (transformation == "logit") {
    logit_transform(y = y_vector)
    
  } else if (transformation == "probit") {
    probit_transform(y = y_vector)
    
  } else if (transformation == "arcsin") {
    arcsin_transform(y = y_vector)
    
  } else {
    stop("Unknown transformation: ", transformation, call. = FALSE)
  }
  
  smp_data[paste(fixed[2])] <- transformed$y
  
  return(list(
    transformed_data = smp_data,
    shift = transformed$shift
  ))
}


std_data_transformation <- function(fixed = fixed, smp_data, transformation, lambda) {
  
  y_vector <- as.matrix(smp_data[paste(fixed[2])])
  
  std_transformed <- if (transformation == "box.cox") {
    as.data.frame(box_cox_std(y = y_vector, lambda = lambda))
    
  } else if (transformation == "dual") {
    as.data.frame(dual_std(y = y_vector, lambda = lambda))
    
  } else if (transformation == "log.shift") {
    as.data.frame(log_shift_opt_std(y = y_vector, lambda = lambda))
    
  } else if (transformation == "log") {
    smp_data[paste(fixed[2])]
    
  } else if (transformation == "no") {
    smp_data[paste(fixed[2])]
    
  } else if (transformation == "logit") {
    as.data.frame(logit_transform(y = y_vector)$y)
    
  } else if (transformation == "probit") {
    as.data.frame(probit_transform(y = y_vector)$y)
    
  } else if (transformation == "arcsin") {
    as.data.frame(arcsin_transform(y = y_vector)$y)
    
  } else {
    stop("Unknown transformation: ", transformation, call. = FALSE)
  }
  
  smp_data[paste(fixed[2])] <- std_transformed
  
  return(transformed_data = smp_data)
}


back_transformation <- function(y, transformation, lambda, shift) {
  
  back_transformed <- if (transformation == "no") {
    no_transform_back(y = y)
    
  } else if (transformation == "log") {
    log_transform_back(y = y, shift = shift)
    
  } else if (transformation == "box.cox") {
    box_cox_back(y = y, lambda = lambda, shift = shift)
    
  } else if (transformation == "dual") {
    dual_back(y = y, lambda = lambda, shift = shift)
    
  } else if (transformation == "log.shift") {
    log_shift_opt_back(y = y, lambda = lambda)
    
  } else if (transformation == "logit") {
    logit_transform_back(y = y)
    
  } else if (transformation == "probit") {
    probit_transform_back(y = y)
    
  } else if (transformation == "arcsin") {
    arcsin_transform_back(y = y)
    
  } else {
    stop("Unknown transformation: ", transformation, call. = FALSE)
  }
  
  return(y = back_transformed)
}


optimal_parameter <- function(generic_opt, fixed, smp_data, smp_domains,
                              transformation, interval) {
  
  # logit, probit and arcsin have no lambda parameter
  no_lambda_transformations <- c("no", "log", "logit", "probit", "arcsin")
  
  if (!transformation %in% no_lambda_transformations) {
    
    if (transformation == "box.cox" && any(interval == "default")) {
      interval <- c(-1, 2)
      
    } else if (transformation == "dual" && any(interval == "default")) {
      interval <- c(0, 2)
      
    } else if (transformation == "log.shift" && any(interval == "default")) {
      
      span <- range(smp_data[paste(fixed[2])])
      
      if ((span[1] + 1) <= 1) {
        lower <- abs(span[1]) + 1
      } else {
        lower <- 0
      }
      
      upper <- diff(span) / 2
      interval <- c(lower, upper)
    }
    
    optimal_parameter <- optimize(
      generic_opt,
      fixed = fixed,
      smp_data = smp_data,
      smp_domains = smp_domains,
      transformation = transformation,
      interval = interval,
      maximum = FALSE
    )$minimum
    
  } else {
    optimal_parameter <- NULL
  }
  
  return(optimal_parameter)
}

# ============================================================
# Modified ebp_check2
# Allows logit, probit and arcsin
# ============================================================

ebp_check2 <- function(threshold, transformation, interval, MSE, boot_type,
                       B, custom_indicator, cpus, seed, na.rm,
                       weights, pop_weights) {
  
  if (!transformation %in% allowed_transformations_ebp2) {
    stop(
      paste0(
        "Argument 'transformation' must be one of: ",
        paste(allowed_transformations_ebp2, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Use original emdi checks for all other arguments.
  # For the new fixed transformations, pass "no" to the original checker.
  transformation_for_original_check <- if (
    transformation %in% c("logit", "probit", "arcsin")
  ) {
    "no"
  } else {
    transformation
  }
  
  old_ebp_check2 <- getFromNamespace("ebp_check2", "emdi")
  
  old_ebp_check2(
    threshold = threshold,
    transformation = transformation_for_original_check,
    interval = interval,
    MSE = MSE,
    boot_type = boot_type,
    B = B,
    custom_indicator = custom_indicator,
    cpus = cpus,
    seed = seed,
    na.rm = na.rm,
    weights = weights,
    pop_weights = pop_weights
  )
  
  invisible(TRUE)
}

mse_estim2 <- function (framework, lambda, shift, model_par, gen_model, res_s, 
          fitted_s, fixed, transformation, interval, L, boot_type) 
{
  if (boot_type == "wild") {
    superpop <- superpopulation_wild(framework = framework, 
                                     model_par = model_par, gen_model = gen_model, lambda = lambda, 
                                     shift = shift, transformation = transformation, 
                                     res_s = res_s, fitted_s = fitted_s)
  }
  else {
    superpop <- superpopulation(framework = framework, model_par = model_par, 
                                gen_model = gen_model, lambda = lambda, shift = shift, 
                                transformation = transformation)
  }
  pop_income_vector <- superpop$pop_income_vector
  if (inherits(framework$threshold, "function")) {
    framework$threshold <- framework$threshold(y = pop_income_vector)
  }
  if (!is.null(framework$aggregate_to_vec)) {
    N_dom_pop_tmp <- framework$N_dom_pop_agg
    pop_domains_vec_tmp <- framework$aggregate_to_vec
  }
  else {
    N_dom_pop_tmp <- framework$N_dom_pop
    pop_domains_vec_tmp <- framework$pop_domains_vec
  }
  if (!is.null(framework$pop_weights)) {
    pop_weights_vec <- framework$pop_data[[framework$pop_weights]]
  }
  else {
    pop_weights_vec <- rep(1, nrow(framework$pop_data))
  }
  true_indicators <- matrix(nrow = N_dom_pop_tmp, data = unlist(lapply(framework$indicator_list, 
                                                                       function(f, threshold) {
                                                                         matrix(nrow = N_dom_pop_tmp, data = unlist(mapply(y = split(pop_income_vector, 
                                                                                                                                     pop_domains_vec_tmp), pop_weights = split(pop_weights_vec, 
                                                                                                                                                                               pop_domains_vec_tmp), f, threshold = framework$threshold)), 
                                                                                byrow = TRUE)
                                                                       }, threshold = framework$threshold)))
  colnames(true_indicators) <- framework$indicator_names
  if (boot_type == "wild") {
    bootstrap_sample <- bootstrap_par_wild(fixed = fixed, 
                                           transformation = transformation, framework = framework, 
                                           model_par = model_par, lambda = lambda, shift = shift, 
                                           vu_tmp = superpop$vu_tmp, res_s = res_s, fitted_s = fitted_s)
  }
  else {
    bootstrap_sample <- bootstrap_par(fixed = fixed, transformation = transformation, 
                                      framework = framework, model_par = model_par, lambda = lambda, 
                                      shift = shift, vu_tmp = superpop$vu_tmp)
  }
  framework$smp_data <- bootstrap_sample
  bootstrap_point_estim <- as.matrix(point_estim(fixed = fixed, 
                                                 transformation = transformation, interval = interval, 
                                                 L = L, framework = framework)[[1]][, -1])
  return(list(
    mse_component = (bootstrap_point_estim - true_indicators)^2,
    boot_estimates = bootstrap_point_estim,
    true_indicators = true_indicators
  ))
  }


mse_estim_wrapper2 <- function(i, B, framework, lambda, shift, model_par,
                               gen_model, fixed, transformation, interval, L,
                               res_s, fitted_s, start_time, boot_type, seedvec = NULL) {
  
  tmp <- mse_estim2(
    framework = framework,
    lambda = lambda,
    shift = shift,
    model_par = model_par,
    gen_model = gen_model,
    res_s = res_s,
    fitted_s = fitted_s,
    fixed = fixed,
    transformation = transformation,
    interval = interval,
    L = L,
    boot_type = boot_type
  )
  
  if (i %% 10 == 0) {
    if (i != B) {
      delta <- difftime(Sys.time(), start_time, units = "secs")
      remaining <- (delta / i) * (B - i)
      remaining <- unclass(remaining)
      remaining <- sprintf(
        "%02d:%02d:%02d:%02d",
        remaining %/% 86400,
        remaining %% 86400 %/% 3600,
        remaining %% 3600 %/% 60,
        remaining %% 60 %/% 1
      )
      
      message(
        "\r", i, " of ", B,
        " Bootstrap iterations completed \t Approximately ",
        remaining, " remaining \n"
      )
      
      if (.Platform$OS.type == "windows") flush.console()
    }
  }
  
  return(tmp)
}

parametric_bootstrap2 <- function(framework, point_estim, fixed,
                                  transformation, interval = c(-1, 2),
                                  L, B, boot_type, parallel_mode, cpus) {
  
  message("\r", "Bootstrap started ")
  
  if (boot_type == "wild") {
    res_s <- residuals(point_estim$model)
    fitted_s <- fitted(point_estim$model, level = 1)
  } else {
    res_s <- NULL
    fitted_s <- NULL
  }
  
  start_time <- Sys.time()
  
  if (cpus > 1) {
    cpus <- min(cpus, parallel::detectCores())
    
    parallelMap::parallelStart(
      mode = parallel_mode,
      cpus = cpus,
      show.info = FALSE
    )
    
    if (parallel_mode == "socket") {
      parallel::clusterSetRNGStream()
    }
    
    parallelMap::parallelLibrary("nlme")
    
    boot_out <- parallelMap::parallelLapply(
      xs = seq_len(B),
      fun = mse_estim_wrapper2,
      B = B,
      framework = framework,
      lambda = point_estim$optimal_lambda,
      shift = point_estim$shift_par,
      model_par = point_estim$model_par,
      gen_model = point_estim$gen_model,
      fixed = fixed,
      transformation = transformation,
      interval = interval,
      L = L,
      res_s = res_s,
      fitted_s = fitted_s,
      start_time = start_time,
      boot_type = boot_type
    )
    
    parallelMap::parallelStop()
    
  } else {
    
    boot_out <- lapply(
      X = seq_len(B),
      FUN = mse_estim_wrapper2,
      B = B,
      framework = framework,
      lambda = point_estim$optimal_lambda,
      shift = point_estim$shift_par,
      model_par = point_estim$model_par,
      gen_model = point_estim$gen_model,
      fixed = fixed,
      transformation = transformation,
      interval = interval,
      L = L,
      res_s = res_s,
      fitted_s = fitted_s,
      start_time = start_time,
      boot_type = boot_type
    )
  }
  
  message("\r", "Bootstrap completed", "\n")
  if (.Platform$OS.type == "windows") flush.console()
  
  mse_array <- simplify2array(
    lapply(boot_out, `[[`, "mse_component")
  )
  
  boot_est_array <- simplify2array(
    lapply(boot_out, `[[`, "boot_estimates")
  )
  
  true_ind_array <- simplify2array(
    lapply(boot_out, `[[`, "true_indicators")
  )
  
  mses <- apply(mse_array, c(1, 2), mean)
  
  if (is.null(framework$aggregate_to_vec)) {
    domains <- unique(framework$pop_domains_vec)
  } else {
    domains <- unique(framework$aggregate_to_vec)
  }
  
  mses <- data.frame(Domain = domains, mses)
  
  return(list(
    MSE = mses,
    boot_estimates = boot_est_array,
    true_indicators = true_ind_array
  ))
}

ebp2 <- function(fixed, pop_data, pop_domains, smp_data, smp_domains, 
                 L = 50, threshold = NULL, transformation = "box.cox", 
                 interval = "default", MSE = FALSE, B = 50, seed = 123, 
                 boot_type = "parametric", 
                 parallel_mode = ifelse(grepl("windows", .Platform$OS.type), 
                                        "socket", "multicore"), 
                 cpus = 1, custom_indicator = NULL, 
                 na.rm = FALSE, weights = NULL, pop_weights = NULL, 
                 aggregate_to = NULL) {
  
  ebp_check1(
    fixed = fixed, 
    pop_data = pop_data, 
    pop_domains = pop_domains, 
    smp_data = smp_data, 
    smp_domains = smp_domains, 
    L = L
  )
  
  ebp_check2(
    threshold = threshold, 
    transformation = transformation, 
    interval = interval, 
    MSE = MSE, 
    boot_type = boot_type, 
    B = B, 
    custom_indicator = custom_indicator, 
    cpus = cpus, 
    seed = seed, 
    na.rm = na.rm, 
    weights = weights, 
    pop_weights = pop_weights
  )
  
  call <- match.call()
  
  if (inherits(call$fixed, "name")) {
    call$fixed <- fixed
  }
  
  if (!is.null(seed)) {
    if (cpus > 1 && parallel_mode != "socket") {
      RNG_kind <- RNGkind()
      set.seed(seed, kind = "L'Ecuyer")
    } else {
      set.seed(seed)
    }
  }
  
  framework <- framework_ebp(
    pop_data = pop_data, 
    pop_domains = pop_domains, 
    smp_data = smp_data, 
    smp_domains = smp_domains, 
    aggregate_to = aggregate_to, 
    custom_indicator = custom_indicator, 
    fixed = fixed, 
    threshold = threshold, 
    na.rm = na.rm, 
    weights = weights, 
    pop_weights = pop_weights
  )
  
  point_estim <- point_estim(
    framework = framework, 
    fixed = fixed, 
    transformation = transformation, 
    interval = interval, 
    L = L, 
    keep_data = TRUE
  )
  
  if (MSE == TRUE) {
    
    mse_estimates <- parametric_bootstrap2(
      framework = framework, 
      point_estim = point_estim, 
      fixed = fixed, 
      transformation = transformation, 
      interval = interval, 
      L = L, 
      B = B, 
      boot_type = boot_type, 
      parallel_mode = parallel_mode, 
      cpus = cpus
    )
    
    ebp_out <- list(
      ind = point_estim$ind,
      
      # MSE classico, uguale come struttura a quello originale
      MSE = mse_estimates$MSE,
      
      # Nuovi oggetti salvati
      boot_estimates = mse_estimates$boot_estimates,
      true_indicators_boot = mse_estimates$true_indicators,
      
      transform_param = point_estim[c(
        "optimal_lambda", 
        "shift_par"
      )],
      
      model = point_estim$model,
      
      framework = framework[c(
        "N_dom_unobs", 
        "N_dom_smp", 
        "N_smp", 
        "N_pop", 
        "smp_domains", 
        "smp_data", 
        "smp_domains_vec", 
        "pop_domains_vec", 
        "response"
      )],
      
      transformation = transformation,
      method = "reml",
      fixed = fixed,
      call = call,
      successful_bootstraps = NULL
    )
    
  } else {
    
    ebp_out <- list(
      ind = point_estim$ind,
      MSE = NULL,
      
      # Li metto comunque, così l'oggetto ha sempre la stessa struttura
      boot_estimates = NULL,
      true_indicators_boot = NULL,
      
      transform_param = point_estim[c(
        "optimal_lambda", 
        "shift_par"
      )],
      
      model = point_estim$model,
      
      framework = framework[c(
        "N_dom_unobs", 
        "N_dom_smp", 
        "N_smp", 
        "N_pop", 
        "smp_domains", 
        "smp_data", 
        "smp_domains_vec", 
        "pop_domains_vec", 
        "response"
      )],
      
      transformation = transformation,
      method = "reml",
      fixed = fixed,
      call = call,
      successful_bootstraps = NULL
    )
  }
  
  if (cpus > 1 && parallel_mode != "socket") {
    RNGkind(RNG_kind[1])
  }
  
  class(ebp_out) <- c("ebp", "emdi")
  
  return(ebp_out)
}

###### Controll variable

check_transformation_compatibility <- function(smp_data, fixed_formula, transformation) {
  
  response <- all.vars(fixed_formula)[1]
  y <- as.numeric(smp_data[[response]])
  
  if (transformation %in% c("logit", "probit", "arcsin")) {
    
    if (any(!is.finite(y), na.rm = TRUE)) {
      return(
        paste0(
          "Transformation '", transformation,
          "' requires finite response values in [0, 1]. ",
          "The response variable contains non-finite values."
        )
      )
    }
    
    if (any(y < 0 | y > 1, na.rm = TRUE)) {
      return(
        paste0(
          "Transformation '", transformation,
          "' is not suitable for these data. ",
          "The response variable must be in [0, 1]. ",
          "Observed range: [",
          round(min(y, na.rm = TRUE), 4), ", ",
          round(max(y, na.rm = TRUE), 4), "]."
        )
      )
    }
  }
  
  if (transformation == "log") {
    
    if (any(y <= 0, na.rm = TRUE)) {
      return(
        paste0(
          "Transformation 'log' is not suitable for these data. ",
          "The response variable must be strictly positive. ",
          "Observed minimum: ",
          round(min(y, na.rm = TRUE), 4), "."
        )
      )
    }
  }
  
  NULL
}