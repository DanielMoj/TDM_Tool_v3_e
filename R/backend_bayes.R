# R/backend_bayes.R
# Backends: Laplace (optim), Stan (cmdstanr/rstan), JAGS (rjags)
suppressPackageStartupMessages({
  library(numDeriv)
})

backend_status <- function() {
  has_cmdstan <- requireNamespace("cmdstanr", quietly = TRUE)
  has_rstan   <- requireNamespace("rstan", quietly = TRUE)
  has_rjags   <- requireNamespace("rjags", quietly = TRUE)
  glue::glue("cmdstanr: {has_cmdstan}, rstan: {has_rstan}, rjags: {has_rjags}")
}

# Prior-Parsing: erwartet priors$theta_log (mu/sd auf log-Skala)
draw_from_priors <- function(priors) {
  th <- priors$theta
  exp(setNames(rnorm(length(th), priors$theta_log$mu[names(th)], priors$theta_log$sd[names(th)]), names(th)))
}

# Helper function for covariate adjustment
apply_covariates <- function(theta, covariates, drug_name = NULL) {
  # Simple allometric scaling for CL and V
  if (!is.null(covariates$weight) && !is.null(theta[["CL"]])) {
    theta[["CL"]] <- theta[["CL"]] * (covariates$weight / 70)^0.75
  }
  if (!is.null(covariates$weight) && !is.null(theta[["Vc"]])) {
    theta[["Vc"]] <- theta[["Vc"]] * (covariates$weight / 70)^1.0
  }
  theta
}

# Helper function for creatinine penalty
cl_creatinine_penalty <- function(CL, age, weight, sex, creatinine_data) {
  # Simplified - returns 0 if no creatinine data
  if (is.null(creatinine_data) || nrow(creatinine_data) == 0) return(0)
  # Placeholder for actual implementation
  0
}

# Helper function for pediatric prior adjustment
adjust_priors_for_covariates <- function(priors, covariates) {
  # Placeholder - returns unchanged priors
  priors
}

# Helper to select appropriate Stan model file
stan_file_for_model <- function(model_type) {
  if (model_type == "MM-1C") {
    return("models/stan/pk_mm_onecpt_ode.stan")
  } else if (model_type == "TMDD-QSS-1C") {
    return("models/stan/pk_tmdd_qss_onecpt_ode.stan")
  } else {
    return("models/stan/pk_multicpt_ode.stan")
  }
}

# Helper to build Stan data list
stan_data_list2 <- function(obs, regimen, priors, model_type, error_model, estimate_sigma, sigma_init, blq_lloq, is_blq) {
  # Map error model string to numeric code
  error_code <- switch(error_model,
    "additiv" = 1L,
    "proportional" = 2L,
    "kombiniert" = 3L,
    "t-additiv" = 4L,
    "t-proportional" = 5L,
    "mixture" = 6L,
    3L  # default to combined
  )
  
  # Build infusion schedule
  n_inf <- regimen$n_doses
  t0 <- regimen$start_time + (0:(n_inf-1)) * regimen$tau
  tinf <- rep(regimen$tinf, n_inf)
  rate <- rep(regimen$dose / regimen$tinf, n_inf)
  
  # BLQ handling
  if (is.null(is_blq)) is_blq <- rep(0L, nrow(obs))
  if (is.na(blq_lloq)) blq_lloq <- 0.0
  
  # Parameter indices for priors
  param_names <- names(priors$theta)
  idx_CL <- which(param_names == "CL")
  idx_Vc <- which(param_names == "Vc")
  idx_Q1 <- which(param_names == "Q1")
  idx_Vp1 <- which(param_names == "Vp1")
  idx_Q2 <- which(param_names == "Q2")
  idx_Vp2 <- which(param_names == "Vp2")
  idx_Vmax <- which(param_names == "Vmax")
  idx_Km <- which(param_names == "Km")
  idx_kint <- which(param_names == "kint")
  idx_Rtot <- which(param_names == "Rtot")
  idx_Kss <- which(param_names == "Kss")
  
  # Set n_cmt based on model
  if (model_type == "1C" || model_type == "MM-1C" || model_type == "TMDD-QSS-1C") {
    n_cmt <- 1L
  } else if (model_type == "2C") {
    n_cmt <- 2L
  } else if (model_type == "3C") {
    n_cmt <- 3L
  } else {
    n_cmt <- 1L
  }
  
  data_list <- list(
    N = nrow(obs),
    t_obs = obs$time,
    y = obs$conc,
    is_blq = as.integer(is_blq),
    lloq = blq_lloq,
    n_cmt = n_cmt,
    n_inf = n_inf,
    t0 = t0,
    tinf = tinf,
    rate = rate,
    error_model = error_code,
    est_sigma = as.integer(estimate_sigma),
    nu = 4.0,  # df for t-distribution
    mix_w = 0.9,  # mixture weight
    mix_scale = 3.0,  # mixture scale factor
    K = length(param_names),
    mu = unname(priors$theta_log$mu[param_names]),
    sd = unname(priors$theta_log$sd[param_names]),
    sigma_add_init = sigma_init$add %||% 1.0,
    sigma_prop_init = sigma_init$prop %||% 0.1
  )
  
  # Add indices (use 1 if not found, Stan will ignore if not needed)
  data_list$idx_CL <- if(length(idx_CL) > 0) idx_CL else 1L
  data_list$idx_Vc <- if(length(idx_Vc) > 0) idx_Vc else 1L
  data_list$idx_Q1 <- if(length(idx_Q1) > 0) idx_Q1 else 1L
  data_list$idx_Vp1 <- if(length(idx_Vp1) > 0) idx_Vp1 else 1L
  data_list$idx_Q2 <- if(length(idx_Q2) > 0) idx_Q2 else 1L
  data_list$idx_Vp2 <- if(length(idx_Vp2) > 0) idx_Vp2 else 1L
  data_list$idx_Vmax <- if(length(idx_Vmax) > 0) idx_Vmax else 1L
  data_list$idx_Km <- if(length(idx_Km) > 0) idx_Km else 1L
  data_list$idx_kint <- if(length(idx_kint) > 0) idx_kint else 1L
  data_list$idx_Rtot <- if(length(idx_Rtot) > 0) idx_Rtot else 1L
  data_list$idx_Kss <- if(length(idx_Kss) > 0) idx_Kss else 1L
  
  data_list
}

neg_log_post_map <- function(par_log, obs, regimen, priors, model_type, error_model, sigma_add, sigma_prop, covariates, blq_lloq = NA_real_, is_blq = NULL, creatinine_data = NULL) {
  # Parameter auf natürlicher Skala
  th_names <- names(priors$theta)
  th <- exp(setNames(par_log, th_names))
  th <- apply_covariates(th, covariates, priors$name)
  # Vorhersage
  pred <- predict_conc_grid(obs$time, regimen, th, model_type)
  # Fehler
  sigfun <- make_sigma_fun(error_model, sigma_add, sigma_prop)
  sig <- sigfun(pred)
  ll <- loglik_residuals_vec(obs$conc, pred, error_model, sigma_add, sigma_prop, lloq = blq_lloq, is_blq = is_blq)
  # joint-like penalty linking CL to creatinine-derived expectation
  ll <- ll + tryCatch(cl_creatinine_penalty(exp(par_log[["logCL"]]), covariates$age %||% 60, covariates$weight %||% 70, covariates$sex %||% "male", creatinine_data), error = function(e) 0)
  # Priors (lognormal)
  lp <- sum(dnorm(par_log, mean = priors$theta_log$mu[th_names], sd = priors$theta_log$sd[th_names], log = TRUE))
  -(ll + lp)
}

run_fit_laplace <- function(obs, regimen, priors, model_type, error_model, covariates, estimate_sigma, sigma_init, blq_lloq = NA_real_, is_blq = NULL) {
  validate_inputs_units(regimen, obs)
  th0 <- log(priors$theta) # Start bei prior-Mean
  sigma_add <- sigma_init[["add"]]; sigma_prop <- sigma_init[["prop"]]
  creatinine_data <- NULL  # placeholder
  obj <- function(p) neg_log_post_map(p, obs, regimen, priors, model_type, error_model, sigma_add, sigma_prop, covariates, blq_lloq, is_blq, creatinine_data)
  opt <- optim(th0, obj, method = "BFGS", hessian = TRUE, control = list(maxit = 1000))
  cov <- tryCatch(solve(opt$hessian), error = function(e) diag(rep(0.05^2, length(th0))))
  draws <- MASS::mvrnorm(n = 800, mu = opt$par, Sigma = cov)
  draws_nat <- exp(draws)
  colnames(draws_nat) <- names(priors$theta)
  list(draws = draws_nat)
}

run_fit_stan <- function(obs, regimen, priors, model_type, error_model, covariates,
                         estimate_sigma, sigma_init, blq_lloq = NA_real_, is_blq = NULL) {
  # --- HMC Controls (from options or defaults) ---
  .hmc <- getOption("tdmx_hmc", default = list(
    chains = 4L, iter_warmup = 1000L, iter_sampling = 1000L,
    parallel_chains = NULL, adapt_delta = 0.9, max_treedepth = 12L, seed = 1234L
  ))

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    warning("cmdstanr nicht verfügbar, fallback auf Laplace.")
    return(run_fit_laplace(obs, regimen, priors, model_type, error_model, covariates, estimate_sigma, sigma_init, blq_lloq, is_blq))
  }
  stan_file <- stan_file_for_model(model_type)
  data_list <- stan_data_list2(obs, regimen, priors, model_type, error_model, estimate_sigma, sigma_init, blq_lloq, is_blq)
  mod <- cmdstanr::cmdstan_model(stan_file)
  fit <- mod$sample(
    data = data_list,
    seed = .hmc$seed,
    chains = .hmc$chains,
    parallel_chains = if (is.null(.hmc$parallel_chains)) .hmc$chains else .hmc$parallel_chains,
    iter_warmup = .hmc$iter_warmup,
    iter_sampling = .hmc$iter_sampling,
    adapt_delta = .hmc$adapt_delta,
    max_treedepth = .hmc$max_treedepth
  )
  draws <- as.data.frame(fit$draws(variables = c("CL_out","Vc_out","Q1_out","Vp1_out","Q2_out","Vp2_out","sigma_add","sigma_prop","nu"), format = "df"))
  # rename outputs to expected names
  names(draws) <- sub("_out$", "", names(draws))
  diagnostics <- NULL
  try({
    if (requireNamespace("posterior", quietly = TRUE)) {
      summ <- posterior::summarise_draws(fit$draws())
      keep <- intersect(c("CL_out","Vc_out","Q1_out","Vp1_out","Q2_out","Vp2_out","sigma_add","sigma_prop","nu"), summ$variable)
      summ <- summ[summ$variable %in% keep, c("variable","rhat","ess_bulk","ess_tail"), drop = FALSE]
    } else { summ <- NULL }
    sdiag <- try(fit$diagnostic_summary(), silent = TRUE)
    div <- try(sdiag$num_divergent[1], silent = TRUE)
    treedepth <- try(sdiag$num_max_treedepth[1], silent = TRUE)
    stepsize <- try(as.numeric(fit$metadata()$step_size_adaptation), silent = TRUE)
    diagnostics <- list(summary = summ, divergences = div, treedepth_hits = treedepth, stepsize = stepsize)
  }, silent = TRUE)
  list(draws = draws, diagnostics = diagnostics)
}

run_fit_stan_advi <- function(obs, regimen, priors, model_type, error_model, covariates, estimate_sigma, sigma_init, blq_lloq = NA_real_, is_blq = NULL) {
  if (!requireNamespace("cmdstanr", quietly = TRUE) && !requireNamespace("rstan", quietly = TRUE)) {
    warning("Stan-ADVI nicht verfügbar, fallback auf Laplace.")
    return(run_fit_laplace(obs, regimen, priors, model_type, error_model, covariates, estimate_sigma, sigma_init, blq_lloq, is_blq))
  }
  stan_file <- stan_file_for_model(model_type)
  data <- stan_data_list2(obs, regimen, priors, model_type, error_model, estimate_sigma, sigma_init, blq_lloq, is_blq)
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    mod <- cmdstanr::cmdstan_model(stan_file)
    fit <- mod$variational(data = data, output_samples = 1000, seed = 123)
    dr <- as.data.frame(fit$draws(variables = c("CL","Vc","Q1","Vp1","Q2","Vp2"), format = "df"))
  } else {
    stan_text <- readChar(stan_file, file.info(stan_file)$size)
    sm <- rstan::stan_model(model_code = stan_text)
    fit <- rstan::vb(sm, data = data, output_samples = 1000, seed = 123)
    dr <- as.data.frame(rstan::extract(fit, pars = c("CL","Vc","Q1","Vp1","Q2","Vp2")))
  }
  keep <- intersect(colnames(dr), names(priors$theta))
  dr <- dr[, keep, drop = FALSE]
  list(draws = dr, diagnostics = NULL)
}

run_fit_jags <- function(obs, regimen, priors, model_type, error_model, covariates,
                         estimate_sigma, sigma_init, blq_lloq = NA_real_, is_blq = NULL) {
  if (!requireNamespace("rjags", quietly = TRUE)) stop("rjags nicht installiert")
  # Data prep
  obs_times <- obs$time; y <- obs$conc
  is_blq <- if (is.null(is_blq)) as.integer(is.na(y)) else as.integer(is_blq)
  y_obs <- ifelse(is.na(y), 0.0, y)
  lloq <- ifelse(is.na(y), blq_lloq, blq_lloq)  # constant lloq for all (if provided)
  grid <- build_time_grid_adaptive(regimen, obs_times, dt_min = 0.025, dt_base = 0.25, refine_window = 0.5)
  # Initial concentration
  C0 <- 0.0
  # Select JAGS model file
  jags_file <- switch(model_type,
    "MM-1C" = "models/jags/pk_mm_onecpt_disc.jags",
    "TMDD-QSS-1C" = "models/jags/pk_tmdd_qss_onecpt_disc.jags",
    # Fallback to linear 1C infusion model if exists
    "1C" = "models/jags/pk_onecpt_inf.jags",
    "models/jags/pk_onecpt_inf.jags"
  )
  data_list <- list(
    N = length(y_obs),
    y = y_obs,
    is_blq = is_blq,
    lloq = rep(blq_lloq, length(y_obs)),
    NT = length(grid$t),
    dt = grid$dt,
    rate = grid$rate,
    idx = grid$idx,
    C0 = C0
  )
  # Inits (log-scale priors)
  inits <- function() {
    lst <- list(
      logCL = log(priors$theta$CL %||% 5),
      logVc = log(priors$theta$Vc %||% 30),
      sigma = sigma_init %||% 2
    )
    if (model_type == "MM-1C") {
      lst$logVmax <- log(priors$theta$Vmax %||% 500)
      lst$logKm <- log(priors$theta$Km %||% 10)
    }
    if (model_type == "TMDD-QSS-1C") {
      lst$logKint <- log(priors$theta$kint %||% 0.1)
      lst$logRtot <- log(priors$theta$Rtot %||% 50)
      lst$logKss <- log(priors$theta$Kss %||% 10)
    }
    lst
  }
  # Parameters to monitor
  params <- c("CL","Vc","sigma")
  if (model_type == "MM-1C") params <- c(params, "Vmax","Km")
  if (model_type == "TMDD-QSS-1C") params <- c(params, "kint","Rtot","Kss")
  j <- rjags::jags.model(jags_file, data = data_list, inits = inits, n.chains = 3, quiet = TRUE)
  rjags::update(j, n.iter = 1000, progress.bar = "none")
  m <- rjags::coda.samples(j, variable.names = params, n.iter = 2000, thin = 2, progress.bar = "none")
  dr <- as.data.frame(do.call(rbind, m))
  names(dr) <- gsub("^\.", "", names(dr))
  list(draws = dr, diagnostics = NULL)
}

run_fit <- function(obs, regimen, priors, model_type, error_model, covariates, backend, estimate_sigma, sigma_init, blq_lloq = NA_real_, is_blq = NULL, use_cache = TRUE, creatinine_data = NULL) {
  backend <- sub(" .*","", backend) # "Laplace", "Stan", "JAGS"
  # cache key
  key <- cache_key_for_fit(obs, regimen, priors, model_type, error_model, covariates, backend, estimate_sigma, sigma_init)
  if (use_cache) {
    cval <- cache_get(key)
    if (!is.null(cval)) return(cval)
  }
  # adjust priors for covariates (pediatric & weight)
  pri_adj <- tryCatch(adjust_priors_for_covariates(priors, covariates), error = function(e) priors)
  priors <- pri_adj
  res <- switch(backend,
    "Laplace" = run_fit_laplace(obs, regimen, priors, model_type, error_model, covariates, estimate_sigma, sigma_init, blq_lloq, is_blq),
    "Stan"    = run_fit_stan(obs, regimen, priors, model_type, error_model, covariates, estimate_sigma, sigma_init, blq_lloq, is_blq),
    "JAGS"    = run_fit_jags2(obs, regimen, priors, model_type, error_model, covariates, estimate_sigma, sigma_init, blq_lloq, is_blq),
    "Stan-ADVI" = run_fit_stan_advi(obs, regimen, priors, model_type, error_model, covariates, estimate_sigma, sigma_init, blq_lloq, is_blq),
    "Stan-Pathfinder" = run_fit_stan_pathfinder(obs, regimen, priors, model_type, error_model, covariates, estimate_sigma, sigma_init, blq_lloq, is_blq),
    run_fit_laplace(obs, regimen, priors, model_type, error_model, covariates, estimate_sigma, sigma_init)
  )
  draws <- as.data.frame(res$draws)
  # Posterior-Zusammenfassung
  summ <- lapply(draws, function(x) c(median = stats::median(x), q2.5 = quantile(x, 0.025), q97.5 = quantile(x, 0.975)))
  median <- sapply(summ, function(z) z["median"])
  q2.5 <- sapply(summ, function(z) z["q2.5"])
  q97.5 <- sapply(summ, function(z) z["q97.5"])
  out <- list(draws = draws, posterior_summary = list(median = median, q2.5 = q2.5, q97.5 = q97.5), diagnostics = res$diagnostics %||% NULL)
  if (use_cache) cache_put(key, out)
  out
}

# --- ADVI / Variational path (requires cmdstanr) ---
fit_advi <- function(data_list, stan_file, seed = 1234, iter = 20000, output_samples = 1000, cache_key = NULL) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) stop("cmdstanr nicht installiert")
  if (!file.exists(stan_file)) stop("Stan file fehlt: ", stan_file)
  cache_file <- if (!is.null(cache_key)) file.path("cache", paste0("advi_", gsub("[^A-Za-z0-9_]+","_", cache_key), ".rds")) else NULL
  if (!is.null(cache_file) && file.exists(cache_file)) {
    res <- readRDS(cache_file); return(res)
  }
  mod <- cmdstanr::cmdstan_model(stan_file)
  fit <- mod$variational(data = data_list, seed = seed, iter = iter, output_samples = output_samples)
  draws <- as.data.frame(fit$draws())
  diagnostics <- NULL
  try({
    if (requireNamespace("posterior", quietly = TRUE)) {
      summ <- posterior::summarise_draws(fit$draws())
      # Keep key PK params if exist
      keep <- intersect(c("CL","Vc","Q1","Vp1","Q2","Vp2","sigma","nu"), summ$variable)
      summ <- summ[summ$variable %in% keep, c("variable","rhat","ess_bulk","ess_tail"), drop = FALSE]
    } else { summ <- NULL }
    # Sampler diagnostics
    sdiag <- try(fit$diagnostic_summary(), silent = TRUE)
    div <- try(fit$diagnostic_summary()$num_divergent[1], silent = TRUE)
    treedepth <- try(fit$diagnostic_summary()$num_max_treedepth[1], silent = TRUE)
    stepsize <- try(as.numeric(fit$metadata()$step_size_adaptation), silent = TRUE)
    diagnostics <- list(summary = summ, divergences = div, treedepth_hits = treedepth, stepsize = stepsize)
  }, silent = TRUE)
  res <- list(draws = draws, method = "advi", seed = seed)
  dir.create("cache", showWarnings = FALSE)
  if (!is.null(cache_file)) saveRDS(res, cache_file)
  res
}

run_fit_stan_pathfinder <- function(obs, regimen, priors, model_type, error_model, covariates, estimate_sigma, sigma_init, blq_lloq = NA_real_, is_blq = NULL) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    warning("cmdstanr nicht verfügbar, fallback auf ADVI.")
    return(run_fit_laplace(obs, regimen, priors, model_type, error_model, covariates, estimate_sigma, sigma_init, blq_lloq, is_blq))
  }
  stan_file <- stan_file_for_model(model_type)
  data <- stan_data_list2(obs, regimen, priors, model_type, error_model, estimate_sigma, sigma_init, blq_lloq, is_blq)
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    mod <- cmdstanr::cmdstan_model(stan_file)
    fit <- mod$variational(data = data, output_samples = 1000, seed = 123)
    dr <- as.data.frame(fit$draws(variables = c("CL","Vc","Q1","Vp1","Q2","Vp2"), format = "df"))
  } else {
    stan_text <- readChar(stan_file, file.info(stan_file)$size)
    sm <- rstan::stan_model(model_code = stan_text)
    fit <- rstan::vb(sm, data = data, output_samples = 1000, seed = 123)
    dr <- as.data.frame(rstan::extract(fit, pars = c("CL","Vc","Q1","Vp1","Q2","Vp2")))
  }
  keep <- intersect(colnames(dr), names(priors$theta))
  dr <- dr[, keep, drop = FALSE]
  list(draws = dr, diagnostics = NULL)
}

run_fit_jags2 <- function(obs, regimen, priors, model_type, error_model, covariates,
                         estimate_sigma, sigma_init, blq_lloq = NA_real_, is_blq = NULL) {
  if (!requireNamespace("rjags", quietly = TRUE)) stop("rjags nicht installiert")
  # Data prep
  obs_times <- obs$time; y <- obs$conc
  is_blq <- if (is.null(is_blq)) as.integer(is.na(y)) else as.integer(is_blq)
  y_obs <- ifelse(is.na(y), 0.0, y)
  lloq <- ifelse(is.na(y), blq_lloq, blq_lloq)  # constant lloq for all (if provided)
  grid <- build_time_grid_adaptive(regimen, obs_times, dt_min = 0.025, dt_base = 0.25, refine_window = 0.5)
  # Initial concentration
  C0 <- 0.0
  # Select JAGS model file
  jags_file <- switch(model_type,
    "MM-1C" = "models/jags/pk_mm_onecpt_disc.jags",
    "TMDD-QSS-1C" = "models/jags/pk_tmdd_qss_onecpt_disc.jags",
    # Fallback to linear 1C infusion model if exists
    "1C" = "models/jags/pk_onecpt_inf.jags",
    "models/jags/pk_onecpt_inf.jags"
  )
  data_list <- list(
    N = length(y_obs),
    y = y_obs,
    is_blq = is_blq,
    lloq = rep(blq_lloq, length(y_obs)),
    NT = length(grid$t),
    dt = grid$dt,
    rate = grid$rate,
    idx = grid$idx,
    C0 = C0
  )
  # Inits (log-scale priors)
  inits <- function() {
    lst <- list(
      logCL = log(priors$theta$CL %||% 5),
      logVc = log(priors$theta$Vc %||% 30),
      sigma = sigma_init %||% 2
    )
    if (model_type == "MM-1C") {
      lst$logVmax <- log(priors$theta$Vmax %||% 500)
      lst$logKm <- log(priors$theta$Km %||% 10)
    }
    if (model_type == "TMDD-QSS-1C") {
      lst$logKint <- log(priors$theta$kint %||% 0.1)
      lst$logRtot <- log(priors$theta$Rtot %||% 50)
      lst$logKss <- log(priors$theta$Kss %||% 10)
    }
    lst
  }
  # Parameters to monitor
  params <- c("CL","Vc","sigma")
  if (model_type == "MM-1C") params <- c(params, "Vmax","Km")
  if (model_type == "TMDD-QSS-1C") params <- c(params, "kint","Rtot","Kss")
  j <- rjags::jags.model(jags_file, data = data_list, inits = inits, n.chains = 3, quiet = TRUE)
  rjags::update(j, n.iter = 1000, progress.bar = "none")
  m <- rjags::coda.samples(j, variable.names = params, n.iter = 2000, thin = 2, progress.bar = "none")
  dr <- as.data.frame(do.call(rbind, m))
  names(dr) <- gsub("^\.", "", names(dr))
  list(draws = dr, diagnostics = NULL)
}