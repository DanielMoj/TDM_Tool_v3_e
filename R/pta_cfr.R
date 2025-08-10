# R/pta_cfr.R
# PTA/CFR computation

# Simulate steady-state over a long horizon and then compute last interval metrics
ss_window <- function(regimen, n_intervals = 20L, dt = 0.05) {
  list(
    times = seq(0, n_intervals * regimen$tau, by = dt),
    t_end  = n_intervals * regimen$tau
  )
}

compute_metrics_for_draw <- function(theta, regimen, model_type, MIC) {
  # FIX: Added missing variable definition
  n_intervals <- 20L
  
  # simulate steady-state and take last interval
  win <- ss_window(regimen)
  conc <- predict_conc_grid(win$times, list(dose = regimen$dose, tau = regimen$tau, tinf = regimen$tinf,
                                            n_doses = n_intervals, start_time = 0),
                            theta, model_type)
  # apply site penetration factor
  drg <- getOption("current_drug_name", default = "Drug")
  st  <- getOption("current_site_name", default = "Plasma")
  conc <- apply_site_penetration(conc, drg, st, load_tissue_cfg("config/tissue.json"))
  # last interval [t_end - tau, t_end]
  t0 <- tail(win$times, 1) - regimen$tau
  idx <- which(win$times >= t0 - 1e-9)
  t_iv <- win$times[idx]
  c_iv <- conc[idx]
  # fT>MIC
  ft <- mean(c_iv > MIC)
  # AUC over last interval (tau), trapezoid
  auc_tau <- sum(diff(t_iv) * zoo::rollmean(c_iv, 2))
  # scale to 24h
  auc24 <- auc_tau * (24 / regimen$tau)
  # Cmax on last interval
  cmax <- max(c_iv)
  list(
    ft_gt_mic = ft,
    auc_tau = auc_tau,
    auc24 = auc24,
    auc24_over_mic = ifelse(MIC > 0, auc24 / MIC, Inf),
    cmax = cmax,
    cmax_over_mic = ifelse(MIC > 0, cmax / MIC, Inf)
  )
}

pta_for_regimen <- function(draws, regimen, model_type, target_def, MIC) {
  if (nrow(draws) == 0) return(NA_real_)
  ok <- logical(nrow(draws))
  for (i in seq_len(nrow(draws))) {
    th <- as.list(draws[i, , drop = FALSE])
    m <- compute_metrics_for_draw(th, regimen, model_type, MIC)
    ok[i] <- meets_target(m, target_def)
  }
  mean(ok)
}

parse_mic_distribution <- function(txt) {
  # format: "0.25:0.1, 0.5:0.2, 1:0.4, 2:0.3"
  if (is.null(txt) || !nzchar(txt)) return(NULL)
  parts <- unlist(strsplit(txt, ","))
  vals <- sapply(parts, function(p) {
    kv <- trimws(unlist(strsplit(p, ":")))
    if (length(kv) != 2) return(c(NA, NA))
    c(as.numeric(kv[1]), as.numeric(kv[2]))
  })
  if (is.null(dim(vals))) return(NULL)
  df <- data.frame(mic = as.numeric(vals[1,]), p = as.numeric(vals[2,]))
  df <- df[is.finite(df$mic) & is.finite(df$p) & df$p >= 0, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df$p <- df$p / sum(df$p)
  df
}

cfr_for_regimen <- function(draws, regimen, model_type, target_def, mic_dist) {
  if (is.null(mic_dist) || nrow(mic_dist) == 0) return(NA_real_)
  pta_vals <- sapply(seq_len(nrow(mic_dist)), function(i) {
    MIC <- mic_dist$mic[i]
    pta_for_regimen(draws, regimen, model_type, target_def, MIC)
  })
  sum(pta_vals * mic_dist$p)
}

# Explore PTA vs dose grid (for Dosing-Studio basic)
pta_vs_dose_grid <- function(draws, regimen, model_type, target_def, MIC, doses) {
  sapply(doses, function(d) {
    reg <- regimen; reg$dose <- d
    pta_for_regimen(draws, reg, model_type, target_def, MIC)
  })
}


# Backwards-compatibility shim for older tests
cfr_from_micdist <- function(draws, regimen, model_type, target, mic_df) {
  .Deprecated("cfr_for_regimen", package = "tdmx")
  cfr_for_regimen(draws, regimen, model_type, target, mic_df)
}