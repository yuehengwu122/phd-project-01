#' Extract Key Diagnostics from Stan Fit
#' @param fit Stan fit object
#' @return List containing convergence diagnostics
extract_diagnostics <- function(fit) {
  summary_df <- rstan::summary(fit)$summary
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  
  # Calculate mean accept_stat across all chains
  accept_stat_values <- sapply(sampler_params, function(chain) {
    mean(chain[, "accept_stat__"], na.rm = TRUE)
  })
  
  list(
    rhat = summary_df[, "Rhat"],
    n_eff_ratio = summary_df[, "n_eff"] / nrow(rstan::extract(fit)$lp__),
    n_divergent = sum(rstan::get_num_divergent(fit)),
    n_max_treedepth = sum(rstan::get_num_max_treedepth(fit)),
    mean_accept_stat = mean(accept_stat_values, na.rm = TRUE)
  )
}


# ============================================================
# Summary Function
# ============================================================

#' Print Quick Summary of Fit Results
#' @param result Result object from run_analysis()
#' @param key_params Character vector of key parameter names to display
#' @param unstandardize Logical, if TRUE show unstandardized linear coefficients 
#'   in a separate table (beta1_unstd = theta * sigma / x_range). Default: FALSE.
#' @param predictor_names Optional names for predictors (used in unstandardized table)
#' @param X_orig Original (unscaled) predictor matrix for computing unstandardized 
#'   coefficients. If NULL (default), uses result$data$X (which may already be rescaled).
print_quick_summary <- function(result, 
                                key_params = c("delta", "lambda", "sigma", "beta0", "beta1"),
                                unstandardize = FALSE,
                                predictor_names = NULL,
                                X_orig = NULL) {
  
  fit <- result$fit
  elapsed <- rstan::get_elapsed_time(fit)
  
  cat("\n=== MODEL:", result$model_id, "===\n")
  cat("Sample size:", result$data$N, "\n")
  if (!is.null(result$data$P)) {
  cat("Number of predictors:", result$data$P, "\n")
  }
  cat("Fitting time:", round((elapsed[1,1] + elapsed[1,2]) / 3600, 2), "hours\n")
  cat("\n")

  
  # Extract diagnostics
  summary_df <- rstan::summary(fit)$summary
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  
  accept_stat_values <- sapply(sampler_params, function(chain) {
    mean(chain[, "accept_stat__"], na.rm = TRUE)
  })
  
  cat("Convergence diagnostics:\n")
  diag_table <- data.frame(
    "Max Rhat" = sprintf("%.4f", max(summary_df[, "Rhat"], na.rm = TRUE)),
    "Min n_eff/N" = sprintf("%.4f", min(summary_df[, "n_eff"] / nrow(rstan::extract(fit)$lp__), na.rm = TRUE)),
    "Divergences" = sprintf("%d", sum(rstan::get_num_divergent(fit))),
    "Max treedepth" = sprintf("%d", sum(rstan::get_num_max_treedepth(fit))),
    "Mean accept_stat" = sprintf("%.4f", mean(accept_stat_values, na.rm = TRUE)),
    row.names = "Value",
    check.names = FALSE
  )
  print(t(diag_table))
  cat("\n")
  
  # Posterior summary for key parameters
  cat("Posterior summary (key parameters):\n")
  param_rows <- rownames(summary_df)
  key_mask <- sapply(key_params, function(p) {
    grepl(sprintf("^%s($|\\[)", p), param_rows)
  })
  key_mask <- apply(key_mask, 1, any)
  
  if (any(key_mask)) {
    post_summary_rounded <- round(summary_df[key_mask, , drop = FALSE], 2)
    print(post_summary_rounded)
  } else {
    cat("No key parameters found.\n")
  }
  
  # Unstandardized linear coefficients table
 if (unstandardize) {
    cat("\n")
    cat("Unstandardized linear coefficients (beta1_unstd = theta * sigma / x_range):\n")
    
    # Extract posterior draws
    post <- rstan::extract(fit)
    
    if (is.null(post$theta) || is.null(post$sigma)) {
      cat("  Cannot compute: theta or sigma not found in fit.\n")
    } else {
      theta_draws <- post$theta
      sigma_draws <- post$sigma
      
      # Get X for computing x_range: use X_orig if provided, otherwise result$data$X
      if (!is.null(X_orig)) {
        X <- X_orig
      } else {
        X <- result$data$X
      }
      
      if (is.null(X)) {
        cat("  Cannot compute: X not found. Provide X_orig argument.\n")
      } else {
        if (is.vector(X)) X <- matrix(X, ncol = 1)
        P <- ncol(X)
        
        # Compute x_range for each predictor (range of original X)
        x_ranges <- apply(X, 2, function(x) diff(range(x)))
        
        # Set predictor names
        if (is.null(predictor_names)) {
          predictor_names <- colnames(X)
          if (is.null(predictor_names)) {
            predictor_names <- paste0("X", 1:P)
          }
        }
        
        # Compute unstandardized beta1 for each draw: beta1_unstd = theta * sigma / x_range
        # theta_draws is (n_draws x P), sigma_draws is (n_draws)
        if (is.matrix(theta_draws)) {
          beta1_unstd <- sweep(theta_draws, 1, sigma_draws, "*")
          beta1_unstd <- sweep(beta1_unstd, 2, x_ranges, "/")
        } else {
          # Single predictor case
          beta1_unstd <- matrix(theta_draws * sigma_draws / x_ranges[1], ncol = 1)
        }
        
        # Summarize: mean, sd, 2.5%, 50%, 97.5%
        unstd_summary <- data.frame(
          predictor = predictor_names,
          mean = round(colMeans(beta1_unstd), 4),
          sd = round(apply(beta1_unstd, 2, sd), 4),
          `2.5%` = round(apply(beta1_unstd, 2, quantile, 0.025), 4),
          `50%` = round(apply(beta1_unstd, 2, quantile, 0.50), 4),
          `97.5%` = round(apply(beta1_unstd, 2, quantile, 0.975), 4),
          check.names = FALSE
        )
        rownames(unstd_summary) <- NULL
        print(unstd_summary, row.names = FALSE)
      }
    }
  }
}


#' Summary of GP Model Fit
#'
#' Provides a clean, publication-ready summary of a fitted GP model,
#' including linear effects, nonlinear effects, and GP hyperparameters.
#' Similar in spirit to brms output but tailored for GP regression.
#'
#' @param result Result object from fit_model()
#' @param predictor_names Optional names for predictors
#' @param compute_sdr Logical, if TRUE compute Savage-Dickey ratios (default: TRUE)
#' @param digits Number of digits to display (default: 4)
#'
#' @return Invisibly returns a list containing the summary tables
#'
#' @examples
#' \dontrun{
#' summary_gp_fit(fit_result)
#' }
summary_gp_fit <- function(result,
                           predictor_names = NULL,
                           compute_sdr = TRUE,
                           digits = 4) {
  
  fit <- result$fit
  post <- rstan::extract(fit)
  summary_df <- rstan::summary(fit)$summary
  elapsed <- rstan::get_elapsed_time(fit)
  n_draws <- nrow(post$lp__)
  
  # ============================================================
  # Header: Model Info
  # ============================================================
  cat("\n")
  cat("Gaussian Process Regression Model\n")
  cat(paste0(rep("=", 50), collapse = ""), "\n")
  cat("Observations:", result$data$N, "\n")
  if (!is.null(result$data$P)) {
    cat("Predictors:  ", result$data$P, "\n")
  }
  cat("Draws:       ", n_draws, "post-warmup samples\n")
  cat("Fit time:    ", round((elapsed[1,1] + elapsed[1,2]) / 3600, 2), "hours\n")
  cat("\n")
  
  # ============================================================
  # Convergence Summary (compact)
  # ============================================================
  max_rhat <- max(summary_df[, "Rhat"], na.rm = TRUE)
  min_neff <- min(summary_df[, "n_eff"], na.rm = TRUE)
  n_div <- sum(rstan::get_num_divergent(fit))
  
  if (max_rhat > 1.01 || n_div > 0) {
    cat("Convergence: Rhat_max =", sprintf("%.3f", max_rhat), 
        "| ESS_min =", round(min_neff),
        "| Divergences =", n_div, "\n")
    if (max_rhat > 1.01) cat("  Warning: Some Rhat > 1.01\n")
    if (n_div > 0) cat("  Warning:", n_div, "divergent transitions\n")
    cat("\n")
  }
  
  # ============================================================
  # Setup: Get X_orig, delta_prior, and predictor names from result
  # ============================================================
  X_scaled <- result$data$X
  if (is.vector(X_scaled)) X_scaled <- matrix(X_scaled, ncol = 1)
  P <- ncol(X_scaled)
  
  # Get X_orig from result (stored by fit_model)
  X_orig <- result$orig$X
  if (is.null(X_orig)) {
    X_orig <- X_scaled
    warning("X_orig not found in result; using scaled X. Unstandardized coefficients may be incorrect.")
  }
  if (is.vector(X_orig)) X_orig <- matrix(X_orig, ncol = 1)
  
  # Get delta_prior from result (stored by fit_model)
  delta_prior <- result$delta_prior
  if (is.null(delta_prior)) delta_prior <- "p1"
  
  if (is.null(predictor_names)) {
    predictor_names <- colnames(X_orig)
    if (is.null(predictor_names)) {
      predictor_names <- paste0("X", 1:P)
    }
  }
  
  # ============================================================
  # Compute SDR if requested
  # ============================================================
  sdr_linear <- rep(NA, P)
  sdr_nonlinear <- rep(NA, P)
  
  if (compute_sdr) {
    # Try to compute SDR using the bayes_factor.R functions
    sdr_result <- tryCatch({
      compute_savage_dickey_ratios(result, delta_prior = delta_prior)
    }, error = function(e) {
      warning("Could not compute SDR: ", e$message)
      NULL
    })
    
    if (!is.null(sdr_result)) {
      sdr_linear <- as.numeric(sdr_result$SDR_linear)
      sdr_nonlinear <- as.numeric(sdr_result$SDR_nonlinear)
    }
  }
  
  # ============================================================
  # Linear Effects (unstandardized beta from pre-computed coef_orig)
  # ============================================================
  cat("Linear Effects:\n")
  cat(paste0(rep("-", 70), collapse = ""), "\n")
  
  # Get pre-computed coefficients on original scale from fit_model
  if (!is.null(result$coef_orig)) {
    beta0_draws_orig <- result$coef_orig$beta0_draws
    beta1_draws_orig <- result$coef_orig$beta1_draws
  } else {
    # Fallback: compute on the fly if coef_orig not available (legacy results)
    beta1_draws <- post$beta1
    y_bar <- mean(result$data$y)
    x_bar <- colMeans(X_orig)
    x_min <- apply(X_orig, 2, min)
    x_max <- apply(X_orig, 2, max)
    rng <- x_max - x_min
    
    if (is.matrix(beta1_draws)) {
      beta1_draws_orig <- sweep(beta1_draws, 2, rng, "/")
    } else {
      beta1_draws_orig <- matrix(beta1_draws / rng[1], ncol = 1)
    }
    
    if (is.matrix(beta1_draws_orig)) {
      beta0_draws_orig <- as.vector(y_bar - beta1_draws_orig %*% x_bar)
    } else {
      beta0_draws_orig <- y_bar - beta1_draws_orig * x_bar
    }
  }
  
  # Get n_eff and Rhat from beta1 (use these for original scale beta)
  beta1_rows <- grep("^beta1\\[", rownames(summary_df))
  if (length(beta1_rows) == 0) {
    beta1_rows <- grep("^beta1$", rownames(summary_df))
  }
  
  if (length(beta1_rows) > 0) {
    beta1_neff <- summary_df[beta1_rows, "n_eff"]
    beta1_rhat <- summary_df[beta1_rows, "Rhat"]
  } else {
    beta1_neff <- rep(NA, P)
    beta1_rhat <- rep(NA, P)
  }
  
  # Compute log10(SDR) for linear effects
  log_sdr_linear <- ifelse(is.na(sdr_linear) | sdr_linear <= 0, NA, log10(sdr_linear))
  
  # Build intercept row (beta0 is derived, so no individual Rhat/ESS)
  intercept_row <- data.frame(
    Predictor = "Intercept",
    Estimate = round(mean(beta0_draws_orig), digits),
    Est.Error = round(sd(beta0_draws_orig), digits),
    `l-95% CI` = round(quantile(beta0_draws_orig, 0.025), digits),
    `u-95% CI` = round(quantile(beta0_draws_orig, 0.975), digits),
    Rhat = "—",
    ESS = "—",
    logSDR10 = "—",
    check.names = FALSE
  )
  
  # Build predictor rows
  predictor_rows <- data.frame(
    Predictor = predictor_names,
    Estimate = round(colMeans(beta1_draws_orig), digits),
    Est.Error = round(apply(beta1_draws_orig, 2, sd), digits),
    `l-95% CI` = round(apply(beta1_draws_orig, 2, quantile, 0.025), digits),
    `u-95% CI` = round(apply(beta1_draws_orig, 2, quantile, 0.975), digits),
    Rhat = ifelse(is.na(beta1_rhat), "—", sprintf("%.3f", beta1_rhat)),
    ESS = ifelse(is.na(beta1_neff), "—", as.character(round(beta1_neff))),
    logSDR10 = ifelse(is.na(log_sdr_linear), "—", sprintf("%.2f", log_sdr_linear)),
    check.names = FALSE
  )
  
  # Combine intercept and predictors
  linear_table <- rbind(intercept_row, predictor_rows)
  
  print(linear_table, row.names = FALSE, right = FALSE)
  cat("\n")
  
  # ============================================================
  # Gaussian Process Terms
  # ============================================================
  cat("Gaussian Process Terms:\n")
  
  # --- Standardized Magnitude (delta) ---
  cat("  Standardized Magnitude (delta):\n")
  cat(paste0(rep("-", 70), collapse = ""), "\n")
  
  # Get delta draws for display (always use delta for table display)
  delta_draws <- post$delta
  delta_param <- "delta"
  
  # Get n_eff and Rhat from delta
  # Try indexed form first, then scalar form
  delta_rows <- grep("^delta\\[", rownames(summary_df))
  if (length(delta_rows) == 0) {
    delta_rows <- grep("^delta$", rownames(summary_df))
  }
  
  # Extract n_eff and Rhat, with fallback to NA if not found
  if (length(delta_rows) > 0) {
    delta_neff <- summary_df[delta_rows, "n_eff"]
    delta_rhat <- summary_df[delta_rows, "Rhat"]
  } else {
    delta_neff <- rep(NA, P)
    delta_rhat <- rep(NA, P)
  }
  
  if (is.matrix(delta_draws)) {
    delta_est <- colMeans(delta_draws)
    delta_se <- apply(delta_draws, 2, sd)
    delta_l95 <- apply(delta_draws, 2, quantile, 0.025)
    delta_u95 <- apply(delta_draws, 2, quantile, 0.975)
  } else {
    delta_est <- mean(delta_draws)
    delta_se <- sd(delta_draws)
    delta_l95 <- quantile(delta_draws, 0.025)
    delta_u95 <- quantile(delta_draws, 0.975)
  }
  
  # Compute log10(SDR) for nonlinear effects
  log_sdr_nonlinear <- ifelse(is.na(sdr_nonlinear) | sdr_nonlinear <= 0, NA, log10(sdr_nonlinear))
  
  nonlinear_table <- data.frame(
    Predictor = predictor_names,
    Estimate = round(delta_est, digits),
    Est.Error = round(delta_se, digits),
    `l-95% CI` = round(delta_l95, digits),
    `u-95% CI` = round(delta_u95, digits),
    Rhat = ifelse(is.na(delta_rhat), "—", sprintf("%.3f", delta_rhat)),
    ESS = ifelse(is.na(delta_neff), "—", as.character(round(delta_neff))),
    logSDR10 = ifelse(is.na(log_sdr_nonlinear), "—", sprintf("%.2f", log_sdr_nonlinear)),
    check.names = FALSE
  )
  
  print(nonlinear_table, row.names = FALSE, right = FALSE)
  cat("\n")
  
  # --- Length Scale (lambda) ---
  cat("  Length Scale (lambda):\n")
  cat(paste0(rep("-", 70), collapse = ""), "\n")
  
  # Lambda (length-scale)
  lambda_draws <- post$lambda
  lambda_rows <- grep("^lambda\\[", rownames(summary_df))
  if (length(lambda_rows) == 0) {
    lambda_rows <- grep("^lambda$", rownames(summary_df))
  }
  
  lambda_neff <- summary_df[lambda_rows, "n_eff"]
  lambda_rhat <- summary_df[lambda_rows, "Rhat"]
  
  if (is.matrix(lambda_draws)) {
    lambda_est <- colMeans(lambda_draws)
    lambda_se <- apply(lambda_draws, 2, sd)
    lambda_l95 <- apply(lambda_draws, 2, quantile, 0.025)
    lambda_u95 <- apply(lambda_draws, 2, quantile, 0.975)
    
    lambda_table <- data.frame(
      Predictor = predictor_names,
      Estimate = round(lambda_est, digits),
      Est.Error = round(lambda_se, digits),
      `l-95% CI` = round(lambda_l95, digits),
      `u-95% CI` = round(lambda_u95, digits),
      Rhat = round(lambda_rhat, 3),
      ESS = round(lambda_neff),
      check.names = FALSE
    )
  } else {
    lambda_table <- data.frame(
      Predictor = "—",
      Estimate = round(mean(lambda_draws), digits),
      Est.Error = round(sd(lambda_draws), digits),
      `l-95% CI` = round(quantile(lambda_draws, 0.025), digits),
      `u-95% CI` = round(quantile(lambda_draws, 0.975), digits),
      Rhat = round(lambda_rhat, 3),
      ESS = round(lambda_neff),
      check.names = FALSE
    )
  }
  
  print(lambda_table, row.names = FALSE, right = FALSE)
  cat("\n")
  
  # ============================================================
  # Residual
  # ============================================================
  cat("Residual:\n")
  cat(paste0(rep("-", 70), collapse = ""), "\n")
  
  # Sigma (residual SD)
  sigma_draws <- post$sigma
  sigma_rows <- grep("^sigma$", rownames(summary_df))
  sigma_est <- mean(sigma_draws)
  sigma_se <- sd(sigma_draws)
  sigma_l95 <- quantile(sigma_draws, 0.025)
  sigma_u95 <- quantile(sigma_draws, 0.975)
  sigma_neff <- summary_df[sigma_rows, "n_eff"]
  sigma_rhat <- summary_df[sigma_rows, "Rhat"]
  
  sigma_table <- data.frame(
    Parameter = "sigma",
    Estimate = round(sigma_est, digits),
    Est.Error = round(sigma_se, digits),
    `l-95% CI` = round(sigma_l95, digits),
    `u-95% CI` = round(sigma_u95, digits),
    Rhat = round(sigma_rhat, 3),
    ESS = round(sigma_neff),
    check.names = FALSE
  )
  
  print(sigma_table, row.names = FALSE, right = FALSE)
  cat("\n")
  
  # ============================================================
  # Footer: SDR interpretation guide
  # ============================================================
  if (compute_sdr) {
    cat("logSDR10: log10 of Savage-Dickey Ratio (evidence for effect vs. no effect)\n")
  }
  cat(paste0(rep("=", 50), collapse = ""), "\n")
  
 # Return summary tables invisibly
  invisible(list(
    linear = linear_table,
    gp_delta = nonlinear_table,
    gp_lambda = lambda_table,
    residual = sigma_table
  ))
}


# ============================================================
# Plot Functions
# ============================================================

#' Plot Posterior Distribution of Delta
#' @param result Result object from run_analysis()
#' @param param Parameter name to plot (default: "delta")
#' @param xlim Custom x-axis limits (default: NULL, auto)
#' @param predictor_names Optional names for predictors
#' @param bins Number of histogram bins (default: 30)
#' @param trim_quantile Quantile threshold for trimming heavy tails (default: 1, no trimming)
#' @param model_name Optional model name to display in title (default: NULL, uses result$model_name)
#' @param plot_prior Logical, if TRUE plot the prior density curve (default: FALSE)
#' @param prior_scale Scale parameter for half-normal prior (default: 1)
#' @param ylim Custom y-axis limits (default: NULL, auto)
plot_posterior_delta <- function(result, 
                          param = "delta", 
                          xlim = NULL,
                          ylim = NULL,
                          predictor_names = NULL,
                          bins = 30,
                          trim_quantile = 1,
                          model_name = NULL,
                          plot_prior = FALSE,
                          prior_scale = 1) {
  
  fit <- result$fit
  draws <- rstan::extract(fit)[[param]]
  
  if (is.null(draws)) {
    stop(sprintf("Parameter '%s' not found in fit object.", param))
  }

  # Trim data based on quantile threshold
  if (trim_quantile < 1) {
    if (is.matrix(draws)) {
      draws <- apply(draws, 2, function(x) x[x <= quantile(x, trim_quantile)])
    } else {
      draws <- draws[draws <= quantile(draws, trim_quantile)]
    }
  }
  
  # Create plot title with model name
  if (is.null(model_name)) {
    model_name <- result$model_name
  }
  
  if (!is.null(model_name)) {
    title_line1 <- paste0(model_name, ": Posterior of ", param)
  } else {
    title_line1 <- paste0("Posterior of ", param)
  }
  
  if (trim_quantile < 1) {
    plot_title <- paste0(title_line1, "\n(trimmed at ", sprintf("%.0f%%", trim_quantile * 100), ")")
  } else {
    plot_title <- title_line1
  }
  
  # Common theme settings
  common_theme <- theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  # Helper function for half-t prior density
  half_t_density <- function(x) {
    (1/2.7) * dt(x/2.7, df = 4) * 2
  }
  
  # Handle vector vs scalar parameter
  if (is.matrix(draws)) {
    P <- ncol(draws)
    
    # Set predictor names from data
    if (is.null(predictor_names)) {
      predictor_names <- colnames(result$data$X)
      if (is.null(predictor_names)) {
        predictor_names <- paste0(param, "[", 1:P, "]")
      }
    }
    
    # Convert to long format
    draws_df <- data.frame(
      value = as.vector(draws),
      predictor = rep(predictor_names, each = length(draws[, 1]))
    )
    draws_df$predictor <- factor(draws_df$predictor, levels = predictor_names)
    
    p <- ggplot(draws_df, aes(x = value)) +
      geom_histogram(aes(y = after_stat(density)), bins = bins, fill = "black", alpha=0.6, color = "white") +
      facet_wrap(~ predictor, scales = "free", ncol = 3) +
      labs(title = plot_title, x = param, y = "Density") +
      common_theme
    
    # Add prior density curve if requested
    if (plot_prior) {
      prior_df <- do.call(rbind, lapply(predictor_names, function(pn) {
        x_range <- range(draws_df$value[draws_df$predictor == pn], na.rm = TRUE)
        x_seq <- seq(0, x_range[2], length.out = 200)
        data.frame(x = x_seq, y = half_t_density(x_seq), 
                   predictor = factor(pn, levels = predictor_names))
      }))
      
      p <- p + geom_line(data = prior_df, aes(x = x, y = y), 
                         color = "red", linewidth = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4)
    }
    
  } else {
    # Scalar parameter
    draws_df <- data.frame(value = draws)
    
    p <- ggplot(draws_df, aes(x = value)) +
      geom_histogram(aes(y = after_stat(density)), bins = bins, fill = "black", alpha=0.6, color = "white") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      labs(title = plot_title, x = param, y = "Density") +
      common_theme
    
    # Add prior density curve if requested
    if (plot_prior) {
      x_seq <- seq(0, max(draws_df$value, na.rm = TRUE), length.out = 200)
      prior_df <- data.frame(x = x_seq, y = half_t_density(x_seq))
      p <- p + geom_line(data = prior_df, aes(x = x, y = y), 
                         color = "red", linewidth = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4)
    }
  }
  
  # Apply custom xlim and ylim if provided
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  }
  
  p
}


#' Plot Posterior Distribution of Theta (Linear Coefficients)
#' @param result Result object from run_analysis()
#' @param param Parameter name to plot (default: "theta")
#' @param xlim Custom x-axis limits (default: NULL, auto)
#' @param predictor_names Optional names for predictors
#' @param bins Number of histogram bins (default: 30)
#' @param trim_quantile Quantile threshold for trimming heavy tails (default: 1, no trimming)
#' @param model_name Optional model name to display in title (default: NULL, uses result$model_name)
#' @param plot_prior Logical, if TRUE plot the prior density curve (default: FALSE)
#' @param ylim Custom y-axis limits (default: NULL, auto)
plot_posterior_theta <- function(result, 
                          param = "theta", 
                          xlim = NULL,
                          ylim = NULL,
                          predictor_names = NULL,
                          bins = 30,
                          trim_quantile = 1,
                          model_name = NULL,
                          plot_prior = FALSE) {
  
  fit <- result$fit
  draws <- rstan::extract(fit)[[param]]
  
  if (is.null(draws)) {
    stop(sprintf("Parameter '%s' not found in fit object.", param))
  }

  # Trim data based on quantile threshold (symmetric trimming for theta)
  if (trim_quantile < 1) {
    lower_q <- (1 - trim_quantile) / 2
    upper_q <- 1 - lower_q
    if (is.matrix(draws)) {
      draws <- apply(draws, 2, function(x) {
        x[x >= quantile(x, lower_q) & x <= quantile(x, upper_q)]
      })
    } else {
      draws <- draws[draws >= quantile(draws, lower_q) & draws <= quantile(draws, upper_q)]
    }
  }
  
  # Create plot title with model name
  if (is.null(model_name)) {
    model_name <- result$model_name
  }
  
  if (!is.null(model_name)) {
    title_line1 <- paste0(model_name, ": Posterior of ", param)
  } else {
    title_line1 <- paste0("Posterior of ", param)
  }
  
  if (trim_quantile < 1) {
    plot_title <- paste0(title_line1, "\n(trimmed at ", sprintf("%.0f%%", trim_quantile * 100), ")")
  } else {
    plot_title <- title_line1
  }
  
  # Common theme settings
  common_theme <- theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  # Handle vector vs scalar parameter
  if (is.matrix(draws)) {
    P <- ncol(draws)
    
    # Set predictor names from data
    if (is.null(predictor_names)) {
      predictor_names <- colnames(result$data$X)
      if (is.null(predictor_names)) {
        predictor_names <- paste0(param, "[", 1:P, "]")
      }
    }
    
    # Convert to long format
    draws_df <- data.frame(
      value = as.vector(draws),
      predictor = rep(predictor_names, each = length(draws[, 1]))
    )
    draws_df$predictor <- factor(draws_df$predictor, levels = predictor_names)
    
    p <- ggplot(draws_df, aes(x = value)) +
      geom_histogram(aes(y = after_stat(density)), bins = bins, fill = "black", alpha = 0.6, color = "white") +
      facet_wrap(~ predictor, scales = "free", ncol = 3) +
      labs(title = plot_title, x = param, y = "Density") +
      common_theme
    
    # Add prior density curve if requested
    if (plot_prior) {
      # Get prior scales using compute_prior_density_theta logic
      X <- result$data$X
      N <- nrow(X)
      X_c <- scale(X, center = TRUE, scale = FALSE)
      scales_theta <- apply(X_c, 2, function(x_col) {
        sqrt(N) / sqrt(sum(x_col^2))
      })
      
      # Create prior data for each predictor (Cauchy prior)
      prior_df <- do.call(rbind, lapply(seq_along(predictor_names), function(j) {
        pn <- predictor_names[j]
        x_range <- range(draws_df$value[draws_df$predictor == pn], na.rm = TRUE)
        # Ensure x_seq includes 0
        x_min <- min(x_range[1], 0)
        x_max <- max(x_range[2], 0)
        x_seq <- seq(x_min, x_max, length.out = 200)
        # Cauchy density centered at 0 with scale from JZS prior
        y_seq <- dcauchy(x_seq, location = 0, scale = scales_theta[j])
        data.frame(x = x_seq, y = y_seq, predictor = factor(pn, levels = predictor_names))
      }))
      
      p <- p + geom_line(data = prior_df, aes(x = x, y = y), 
                         color = "red", linewidth = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4)
    }
    
  } else {
    # Scalar parameter
    draws_df <- data.frame(value = draws)
    
    p <- ggplot(draws_df, aes(x = value)) +
      geom_histogram(aes(y = after_stat(density)), bins = bins, fill = "black", alpha = 0.6, color = "white") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      labs(title = plot_title, x = param, y = "Density") +
      common_theme
    
    # Add prior density curve if requested
    if (plot_prior) {
      X <- result$data$X
      N <- nrow(X)
      X_c <- scale(X, center = TRUE, scale = FALSE)
      scale_theta <- sqrt(N) / sqrt(sum(X_c^2))
      
      x_range <- range(draws_df$value, na.rm = TRUE)
      # Ensure x_seq includes 0
      x_min <- min(x_range[1], 0)
      x_max <- max(x_range[2], 0)
      x_seq <- seq(x_min, x_max, length.out = 200)
      prior_df <- data.frame(x = x_seq, y = dcauchy(x_seq, location = 0, scale = scale_theta))
      p <- p + geom_line(data = prior_df, aes(x = x, y = y), 
                         color = "red", linewidth = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4)
    }
  }
  
  # Apply custom xlim and ylim if provided
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  }
  
  p
}



#' Plot Estimated GP Trends
#' @param result Result object from run_analysis()
#' @param X Predictor matrix (N x P), should match the data used for fitting
#' @param predictor_names Optional names for predictors
#' @param nrow Number of rows in plot layout (default: auto)
#' @param ncol Number of columns in plot layout (default: auto)
plot_gp_trends_old <- function(result, 
                          X = NULL,
                          predictor_names = NULL,
                          nrow = NULL, 
                          ncol = NULL) {
  
  fit <- result$fit
  y <- result$data$y
  
  # Use X from data if not provided
  if (is.null(X)) {
    X <- result$data$X
  }
  
  # Convert to matrix if single predictor
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }
  
  P <- ncol(X)
  
  # Set predictor names
  if (is.null(predictor_names)) {
    predictor_names <- colnames(X)
    if (is.null(predictor_names)) {
      predictor_names <- paste0("X", 1:P)
    }
  }
  
  # Extract posterior samples
  fitGP_post <- rstan::extract(fit)
  
  # Auto layout
  if (is.null(nrow) && is.null(ncol)) {
    if (P == 1) {
      nrow <- 1; ncol <- 1
    } else if (P <= 3) {
      nrow <- 1; ncol <- P
    } else if (P <= 6) {
      nrow <- 2; ncol <- 3
    } else {
      nrow <- ceiling(sqrt(P))
      ncol <- ceiling(P / nrow)
    }
  }
  
  # Set up plot layout
  old_par <- par(mfrow = c(nrow, ncol), mar = c(3, 3, 2, 1), mgp = c(2, 0.7, 0))
  on.exit(par(old_par))
  
  # Plot each predictor
  for (i in 1:P) {
    # Extract samples for predictor i
    if (P == 1) {
      mu_samples <- fitGP_post$full_est_grid
      lin_samples <- fitGP_post$linear_est_grid
    } else {
      mu_samples <- fitGP_post$full_est_grid[, , i]
      lin_samples <- fitGP_post$linear_est_grid[, , i]
    }
    
    # Calculate summaries
    full_mean <- apply(mu_samples, 2, mean)
    full_lower <- apply(mu_samples, 2, quantile, 0.025)
    full_upper <- apply(mu_samples, 2, quantile, 0.975)
    lin_mean <- apply(lin_samples, 2, mean)
    
    # Grid for plotting
    x_i <- X[, i]
    x_grid <- seq(from = min(x_i), to = max(x_i), length.out = 50)
    
    # Create plot
    plot(x_i, y, pch = 16, col = 1,
         main = predictor_names[i],
         xlab = predictor_names[i],
         ylab = "y")
    
    # Add credible interval
    polygon(
      c(x_grid, rev(x_grid)), 
      c(full_lower, rev(full_upper)),
      col = adjustcolor("blue", alpha.f = 0.2),
      border = NA
    )
    
    # Add lines
    lines(x_grid, lin_mean, col = 1, lty = 2, lwd = 1.5)
    lines(x_grid, full_mean, col = "blue", lwd = 1.5)
    
    # Add legend on first plot only
    if (i == 1) {
      legend("topleft", 
             c("Linear", "Linear + GP"), 
             col = c(1, "blue"), 
             lwd = 1.5, 
             lty = c(2, 1),
             bty = "n",
             cex = 0.8)
    }
  }
}


#' Plot GP Trends for Exact GP Model
#' @param result Result object from fit_model()
#' @param p Optional vector of predictor indices to plot (default: all)
#' @param predictor_names Optional names for predictors
#' @param N_grid Number of grid points for prediction (default: 80)
#' @param ndraws Number of posterior draws to use (default: 400)
#' @param seed Random seed for drawing samples (default: 1)
#' @param model_name Optional model name to display in title (default: NULL, uses result$model_name)
#' @param ylim Optional y-axis limits as a numeric vector of length 2 (default: NULL, auto)
plot_gp_trends <- function(result,
                           p = NULL,
                           predictor_names = NULL,
                           N_grid = 80,
                           ndraws = 400,
                           seed = 1,
                           model_name = NULL,
                           ylim = NULL) {
  stopifnot(!is.null(result$fit))
  stopifnot(!is.null(result$data$X), !is.null(result$data$y))
  stopifnot(!is.null(result$orig$X))  # original X (same columns as Stan X)

  fit    <- result$fit
  X_s    <- as.matrix(result$data$X)  # rescaled X used for Stan input (before centering)
  y      <- as.numeric(result$data$y)
  X_orig <- as.matrix(result$orig$X)  # original X for x-axis

  if (is.vector(X_s))    X_s    <- matrix(X_s, ncol = 1)
  if (is.vector(X_orig)) X_orig <- matrix(X_orig, ncol = 1)

  N <- nrow(X_s)
  P <- ncol(X_s)

  stopifnot(nrow(X_orig) == N, ncol(X_orig) == P)

  if (is.null(p)) p_idx <- seq_len(P) else p_idx <- as.integer(p)

  if (is.null(predictor_names)) {
    predictor_names <- colnames(X_orig)
    if (is.null(predictor_names)) predictor_names <- paste0("X", seq_len(P))
  }

  # ---- Centering exactly as in Stan ----
  y_bar <- mean(y)
  y_c   <- y - y_bar

  x_bar_s <- colMeans(X_s)
  X_c     <- sweep(X_s, 2, x_bar_s, "-")  # N x P

  # ---- Helpers: kernel + projection ----
  se_cov_from_D2 <- function(D2, alpha, l) {
    K <- (alpha^2) * exp(-0.5 * D2 / (l^2))
    diag(K) <- diag(K) + 1e-12
    K
  }

  ktilde_from_K <- function(K, Q, R) {
    U  <- K %*% Q
    PK <- Q %*% (R %*% t(U))
    M  <- t(Q) %*% U
    PKP <- Q %*% (R %*% M %*% R) %*% t(Q)
    Kt <- K - PK - t(PK) + PKP
    0.5 * (Kt + t(Kt))
  }

  ktilde_cross <- function(Kc, Qn, Rn, Qt, Rt) {
    # (I - Pn) Kc (I - Pt) where Kc is G x N
    T    <- Kc %*% Qt
    KcPt <- (T %*% Rt) %*% t(Qt)

    S    <- t(Qn) %*% Kc
    PnKc <- Qn %*% (Rn %*% S)

    S2     <- t(Qn) %*% KcPt
    PnKcPt <- Qn %*% (Rn %*% S2)

    Kc - PnKc - KcPt + PnKcPt
  }

  # ---- Precompute training projectors + D2 on centered-rescaled scale ----
  Q_tr <- vector("list", P)
  R_tr <- vector("list", P)
  D2_tr <- vector("list", P)

  for (j in 1:P) {
    Q <- cbind(1, X_c[, j])
    R <- solve(crossprod(Q))
    Q_tr[[j]] <- Q
    R_tr[[j]] <- R

    xj <- X_c[, j]
    D2 <- outer(xj, xj, "-"); D2 <- D2 * D2; diag(D2) <- 0
    D2_tr[[j]] <- D2
  }

  # ---- Draws ----
  post <- rstan::extract(fit, pars = c("sigma","theta","delta","lambda"))
  S_all <- length(post$sigma)

  set.seed(seed)
  idx <- if (ndraws < S_all) sample.int(S_all, ndraws) else seq_len(S_all)
  nd <- length(idx)

  sigma_draw  <- post$sigma[idx]
  theta_draw  <- post$theta[idx, , drop = FALSE]
  delta_draw  <- post$delta[idx, , drop = FALSE]
  lambda_draw <- post$lambda[idx, , drop = FALSE]

  # scaling map used in preprocessing: x_s = (x - min) / range - 0.5
  x_min <- apply(X_orig, 2, min)
  x_rng <- apply(X_orig, 2, max) - x_min

  out_list <- vector("list", length(p_idx))
  raw_list <- vector("list", length(p_idx))

  for (kk in seq_along(p_idx)) {
    i <- p_idx[kk]
    x_name <- predictor_names[i]

    # grid on ORIGINAL scale for x-axis
    xg_orig <- seq(min(X_orig[, i]), max(X_orig[, i]), length.out = N_grid)

    # transform grid to model scale (rescaled then centered like Stan)
    xg_s <- (xg_orig - x_min[i]) / x_rng[i]
    xg_c <- xg_s - x_bar_s[i]

    # grid-side projector for predictor i (span{1, xg_c})
    Qn <- cbind(1, xg_c)
    Rn <- solve(crossprod(Qn))

    # cross distances on centered-rescaled scale
    D2_cross <- outer(xg_c, X_c[, i], "-"); D2_cross <- D2_cross * D2_cross

    full_mat <- matrix(NA_real_, nd, N_grid)
    lin_mat  <- matrix(NA_real_, nd, N_grid)

    Qt <- Q_tr[[i]]
    Rt <- R_tr[[i]]

    for (s in 1:nd) {
      sigma <- sigma_draw[s]
      beta1 <- sigma * theta_draw[s, ]          # slopes on centered-rescaled scale
      alpha <- sigma * delta_draw[s, ]
      lam   <- lambda_draw[s, ]

      # build training covariance
      Sigma <- matrix(0, N, N)
      for (j in 1:P) {
        Kj <- se_cov_from_D2(D2_tr[[j]], alpha[j], lam[j])
        Sigma <- Sigma + ktilde_from_K(Kj, Q_tr[[j]], R_tr[[j]])
      }
      diag(Sigma) <- diag(Sigma) + sigma^2

      # residual on centered-y scale (NO intercept in centered model)
      resid <- y_c - as.numeric(X_c %*% beta1)

      # a = Sigma^{-1} resid
      U <- chol(Sigma)                      # upper
      v <- forwardsolve(t(U), resid)
      a_vec <- backsolve(U, v)

      # predictor-i GP mean at grid (conditional on data)
      Kc <- (alpha[i]^2) * exp(-0.5 * D2_cross / (lam[i]^2))
      Kc_tilde <- ktilde_cross(Kc, Qn, Rn, Qt, Rt)
      f_mean <- as.numeric(Kc_tilde %*% a_vec)

      # trend on ORIGINAL y scale:
      # y = y_bar + (xg_c * beta1_i) + f(xg_c)
      lin  <- y_bar + beta1[i] * xg_c
      full <- lin + f_mean

      lin_mat[s, ]  <- lin
      full_mat[s, ] <- full
    }

    out_list[[kk]] <- data.frame(
      predictor = factor(x_name, levels = predictor_names[p_idx]),
      x = xg_orig,
      full_mean = colMeans(full_mat),
      full_lo   = apply(full_mat, 2, quantile, 0.025),
      full_hi   = apply(full_mat, 2, quantile, 0.975),
      lin_mean  = colMeans(lin_mat)
    )

    raw_list[[kk]] <- data.frame(
      predictor = factor(x_name, levels = predictor_names[p_idx]),
      x = X_orig[, i],
      y = y
    )
  }

  df  <- do.call(rbind, out_list)
  raw <- do.call(rbind, raw_list)

  if (is.null(model_name)) model_name <- result$model_name
  plot_title <- if (!is.null(model_name)) paste("GP Trends:", model_name) else "GP Trends"

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = raw, ggplot2::aes(x = x, y = y), color = "black", size = 1) +
    ggplot2::geom_ribbon(data = df, ggplot2::aes(x = x, ymin = full_lo, ymax = full_hi),
                         fill = "lightskyblue", alpha = 0.30) +
    ggplot2::geom_line(data = df, ggplot2::aes(x = x, y = lin_mean),
                       color = "black", linetype = 2, linewidth = 0.9) +
    ggplot2::geom_line(data = df, ggplot2::aes(x = x, y = full_mean),
                       color = "blue", linewidth = 1.1) +
    ggplot2::facet_wrap(~ predictor, scales = "free_x") +
    ggplot2::labs(title = plot_title, x = NULL, y = "y") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "none"
    )

  if (!is.null(ylim)) p <- p + ggplot2::coord_cartesian(ylim = ylim)
  p
}


#' Helper function for HSGP basis computation
#' @param x_c Centered predictor values
#' @param M Number of basis functions
#' @param L Boundary factor
phi_grid_ortho <- function(x_c, M, L) {
  N <- length(x_c)
  PHI <- matrix(NA_real_, N, M)
  for (m in 1:M) {
    arg <- m * pi * (x_c + L) / (2 * L)
    PHI[, m] <- sin(arg) / sqrt(L)
  }
  PHI
}


#' Plot HSGP Trends
#' @param object Result object from run_analysis()
#' @param X_orig Original (unscaled) predictor matrix for plotting x-axis (default: NULL, uses X)
#' @param p Optional vector of predictor indices to plot (default: all)
#' @param N_grid Number of grid points for prediction (default: 50)
#' @param ndraws Number of posterior draws to use (default: 1000)
#' @param seed Random seed for drawing samples (default: 1)
#' @param predictor_names Optional names for predictors
#' @param model_name Optional model name to display in title (default: NULL, uses object$model_name)
#' @param ylim Optional y-axis limits as a numeric vector of length 2 (default: NULL, auto)
plot_hsgp_trends <- function(object, X_orig = NULL, p = NULL, N_grid = 50, ndraws = 1000, seed = 1,
                             predictor_names = NULL, model_name = NULL, ylim = NULL) {

  fit <- object$fit
  X   <- as.matrix(object$data$X)
  y   <- object$data$y
  M   <- object$data$M

  N <- nrow(X); P <- ncol(X)
  
  # Use X_orig for plotting if provided, otherwise use X
  if (is.null(X_orig)) {
    X_plot <- X
  } else {
    if (is.vector(X_orig)) X_orig <- matrix(X_orig, ncol = 1)
    X_plot <- as.matrix(X_orig)
  }

  if (is.null(p)) p_idx <- seq_len(P) else p_idx <- as.integer(p)

  if (is.null(predictor_names)) {
    predictor_names <- colnames(X)
    if (is.null(predictor_names)) predictor_names <- paste0("X", seq_len(P))
  }

  # centering like Stan
  y_bar <- mean(y)
  x_bar <- colMeans(X)
  X_c <- sweep(X, 2, x_bar, "-")

  post <- rstan::extract(fit, pars = c("beta0","sigma","theta","delta","lambda","z_basis"))
  S <- length(post$beta0)
  if (S > ndraws) {
    set.seed(seed)
    idx <- sample.int(S, ndraws)
  } else idx <- seq_len(S)

  mu_draw    <- post$beta0[idx]
  sigma_draw <- post$sigma[idx]

  out_list <- vector("list", length(p_idx))
  raw_list <- vector("list", length(p_idx))

  for (kk in seq_along(p_idx)) {
    pp <- p_idx[kk]
    x_name <- predictor_names[pp]

    xcp <- X_c[, pp]
    half_range <- 0.5 * (max(xcp) - min(xcp))
    Lp <- max(1e-6, 1.5 * half_range)

    # grid in rescaled scale (for model calculations)
    x_grid_raw <- seq(min(X[, pp]), max(X[, pp]), length.out = N_grid)
    x_grid_c   <- x_grid_raw - x_bar[pp]
    PHI_g      <- phi_grid_ortho(x_grid_c, M, Lp)
    
    # grid in original scale (for plotting)
    x_grid_plot <- seq(min(X_plot[, pp]), max(X_plot[, pp]), length.out = N_grid)

    theta_p  <- post$theta[idx, pp]
    delta_p  <- post$delta[idx, pp]
    lambda_p <- post$lambda[idx, pp]
    ZB <- post$z_basis[idx, pp, 1:M, drop = FALSE]  # [draw, 1, m]

    m_vec <- 1:M
    omega <- m_vec * pi / (2 * Lp)

    full_mat <- matrix(NA_real_, length(idx), N_grid)
    lin_mat  <- matrix(NA_real_, length(idx), N_grid)

    for (i in seq_along(idx)) {
      beta0_orig <- mu_draw[i] + y_bar
      beta1_p <- sigma_draw[i] * theta_p[i]
      a <- sigma_draw[i] * delta_p[i]
      ell <- lambda_p[i]

      # Matérn-3/2 weights (match Stan)
      s <- (a^2) * ell / (1 + (ell^2) * (omega^2))^2
      w <- sqrt(s) * as.numeric(ZB[i, 1, ])

      f_g    <- as.vector(PHI_g %*% w)
      lin_g  <- beta0_orig + beta1_p * x_grid_c
      full_g <- lin_g + f_g

      lin_mat[i, ]  <- lin_g
      full_mat[i, ] <- full_g
    }

    out_list[[kk]] <- data.frame(
      predictor = factor(x_name, levels = predictor_names[p_idx]),
      x = x_grid_plot,
      full_mean = colMeans(full_mat),
      full_lo = apply(full_mat, 2, quantile, 0.025),
      full_hi = apply(full_mat, 2, quantile, 0.975),
      lin_mean = colMeans(lin_mat)
    )

    raw_list[[kk]] <- data.frame(
      predictor = factor(x_name, levels = predictor_names[p_idx]),
      x = X_plot[, pp],
      y = y
    )
  }

  df  <- do.call(rbind, out_list)
  raw <- do.call(rbind, raw_list)

  # Create plot title with model name if available
  if (is.null(model_name)) {
    model_name <- object$model_name
  }
  
  plot_title <- if (!is.null(model_name)) {
    paste("HSGP Trends:", model_name)
  } else {
    "HSGP Trends"
  }

  p <- ggplot() +
    geom_point(data = raw, aes(x = x, y = y), color = "black", size = 1) +
    geom_ribbon(data = df, aes(x = x, ymin = full_lo, ymax = full_hi),
                fill = "lightskyblue", alpha = 0.30) +
    geom_line(data = df, aes(x = x, y = lin_mean),
              color = "black", linetype = 2, linewidth = 0.9) +
    geom_line(data = df, aes(x = x, y = full_mean),
              color = "blue", linewidth = 1.1) +
    facet_wrap(~ predictor, scales = "free_x") +
    labs(title = plot_title, x = NULL, y = "y") +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(size = 11),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "none", 
      plot.margin = margin(10, 10, 10, 10),
      panel.spacing = unit(1.2, "lines")
    )
  
  # Apply custom ylim if provided
  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  p
}


