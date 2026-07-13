# ============================================================
# Diagnostics, Summary, and Plotting for gp_fit Objects
# ============================================================


# ============================================================
# Convergence Diagnostics
# ============================================================

#' Extract convergence diagnostics from a Stan fit
#' @param fit Stan fit object
#' @return List of diagnostic quantities
extract_diagnostics <- function(fit) {
  summary_df <- rstan::summary(fit)$summary
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)

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
# Print / Summary Methods
# ============================================================

#' Print method for gp_fit objects
#' @param x A gp_fit object
#' @param ... Ignored
#' @export
print.gp_fit <- function(x, ...) {
  cat("GP Regression Fit\n")
  cat("  Model:      ", x$model_type, "\n")
  cat("  d_order:    ", x$d_order, "\n")
  cat("  N:          ", x$stan_data$N, "\n")
  cat("  P:          ", x$stan_data$P, "\n")
  cat("  Predictors: ", paste(x$predictor_names, collapse = ", "), "\n")

  elapsed <- rstan::get_elapsed_time(x$fit)
  total_sec <- max(rowSums(elapsed))
  if (total_sec > 3600) {
    cat("  Fit time:   ", round(total_sec / 3600, 2), "hours\n")
  } else {
    cat("  Fit time:   ", round(total_sec / 60, 1), "minutes\n")
  }

  diag <- extract_diagnostics(x$fit)
  if (diag$n_divergent > 0)
    cat("  WARNING:    ", diag$n_divergent, "divergent transitions\n")
  if (max(diag$rhat, na.rm = TRUE) > 1.01)
    cat("  WARNING:     Max Rhat =", sprintf("%.3f", max(diag$rhat, na.rm = TRUE)), "\n")

  invisible(x)
}


#' Publication-ready summary of a GP fit
#'
#' Prints tables for linear effects (on original scale), GP effects
#' (delta, lambda), residual sigma, and optionally log10(BF10).
#'
#' @param result A gp_fit object from fit_model()
#' @param compute_bf If TRUE, compute Savage-Dickey / interval-null BFs
#' @param eps Epsilon for nonlinear interval-null BF
#'   (default: BF_EPSILON_DEFAULT = 0.2)
#' @param digits Number of digits for rounding
#' @return Invisibly returns a list of the summary tables
summary_gp_fit <- function(
  result,
  compute_bf = TRUE,
  eps = BF_EPSILON_DEFAULT,
  digits = 4
) {
  stopifnot(inherits(result, "gp_fit"))

  fit <- result$fit
  post <- rstan::extract(fit)
  summary_df <- rstan::summary(fit)$summary
  elapsed <- rstan::get_elapsed_time(fit)
  n_draws <- nrow(post$lp__)

  P <- result$stan_data$P
  pred_names <- result$predictor_names
  masks <- get_gp_fit_masks(result)

  # ---- Header ----
  cat("\n")
  cat("Gaussian Process Regression Model\n")
  cat(strrep("=", 60), "\n")
  cat("Observations:", result$stan_data$N, "\n")
  cat("Predictors:  ", P, "\n")
  cat("Model type:  ", result$model_type, "\n")
  cat("Prior order: ", result$d_order, "\n")
  cat("Draws:       ", n_draws, "post-warmup samples\n")
  total_sec <- max(rowSums(elapsed))
  if (total_sec > 3600) {
    cat("Fit time:    ", round(total_sec / 3600, 2), "hours\n")
  } else {
    cat("Fit time:    ", round(total_sec / 60, 1), "minutes\n")
  }
  cat("\n")

  # ---- Convergence ----
  max_rhat <- max(summary_df[, "Rhat"], na.rm = TRUE)
  min_neff <- min(summary_df[, "n_eff"], na.rm = TRUE)
  n_div <- sum(rstan::get_num_divergent(fit))

  if (max_rhat > 1.01 || n_div > 0) {
    cat("Convergence: Rhat_max =",
        sprintf("%.3f", max_rhat), "| ESS_min =",
        round(min_neff), "| Divergences =", n_div, "\n")
    if (max_rhat > 1.01) cat("  WARNING: Some Rhat > 1.01\n")
    if (n_div > 0) cat("  WARNING:", n_div, "divergent transitions\n")
    cat("\n")
  }

  # ---- BF computation (if requested) ----
  log_bf10_linear <- rep(NA_real_, P)
  log_bf10_nonlinear <- rep(NA_real_, P)

  if (compute_bf && is_gp_full_fit(result)) {
    bf_result <- tryCatch(
      compute_bf_table(result, eps = eps),
      error = function(e) {
        warning("BF computation failed: ", e$message)
        NULL
      }
    )
    if (!is.null(bf_result)) {
      log_bf10_linear <- bf_result$log_bf10_linear
      log_bf10_nonlinear <- bf_result$log_bf10_nonlinear
    }
  }

  # ---- Linear Effects (original scale) ----
  cat("Linear Effects (original scale):\n")
  cat(strrep("-", 60), "\n")

  beta0_draws <- result$coef_orig$beta0_draws
  beta1_draws <- result$coef_orig$beta1_draws

  # Rhat / ESS for beta1
  beta1_rows <- grep("^beta1(\\[|$)", rownames(summary_df))
  beta1_neff <- if (length(beta1_rows) > 0) summary_df[beta1_rows, "n_eff"] else rep(NA, P)
  beta1_rhat <- if (length(beta1_rows) > 0) summary_df[beta1_rows, "Rhat"] else rep(NA, P)

  fmt_bf <- function(v) ifelse(is.na(v), "---", sprintf("%.2f", v))

  # Intercept row
  intercept_row <- data.frame(
    Predictor = "Intercept",
    Estimate = round(mean(beta0_draws), digits),
    Est.Error = round(sd(beta0_draws), digits),
    `l-95% CI` = round(quantile(beta0_draws, 0.025), digits),
    `u-95% CI` = round(quantile(beta0_draws, 0.975), digits),
    Rhat = "---",
    ESS = "---",
    `log BF10` = "---",
    check.names = FALSE
  )

  # Predictor rows
  predictor_rows <- data.frame(
    Predictor = pred_names,
    Estimate = round(colMeans(beta1_draws), digits),
    Est.Error = round(apply(beta1_draws, 2, sd), digits),
    `l-95% CI` = round(apply(beta1_draws, 2, quantile, 0.025), digits),
    `u-95% CI` = round(apply(beta1_draws, 2, quantile, 0.975), digits),
    Rhat = ifelse(is.na(beta1_rhat), "---", sprintf("%.3f", beta1_rhat)),
    ESS = ifelse(is.na(beta1_neff), "---", as.character(round(beta1_neff))),
    `log BF10` = fmt_bf(log_bf10_linear),
    check.names = FALSE
  )
  predictor_rows$Rhat[masks$linear == 0L] <- "---"
  predictor_rows$ESS[masks$linear == 0L] <- "---"

  linear_table <- rbind(intercept_row, predictor_rows)
  print(linear_table, row.names = FALSE, right = FALSE)
  cat("\n")

  # ---- GP Terms (full model only) ----
  gp_delta_table <- NULL
  gp_lambda_table <- NULL

  if (is_gp_full_fit(result)) {

    # --- Delta (nonlinear magnitude) ---
    cat("GP Nonlinear Magnitude (delta):\n")
    cat(strrep("-", 60), "\n")

    delta_draws <- extract_gp_draws_full(result, post = post, param = "delta")
    delta_row_pattern <- if (identical(result$model_type, "full_select")) "^delta_full(\\[|$)" else "^delta(\\[|$)"
    delta_rows <- grep(delta_row_pattern, rownames(summary_df))
    delta_neff <- if (length(delta_rows) > 0) summary_df[delta_rows, "n_eff"] else rep(NA, P)
    delta_rhat <- if (length(delta_rows) > 0) summary_df[delta_rows, "Rhat"] else rep(NA, P)

    d_est <- colMeans(delta_draws)
    d_se  <- apply(delta_draws, 2, sd)
    d_lo  <- apply(delta_draws, 2, quantile, 0.025)
    d_hi  <- apply(delta_draws, 2, quantile, 0.975)

    gp_delta_table <- data.frame(
      Predictor = pred_names,
      Estimate = round(d_est, digits),
      Est.Error = round(d_se, digits),
      `l-95% CI` = round(d_lo, digits),
      `u-95% CI` = round(d_hi, digits),
      Rhat = ifelse(is.na(delta_rhat), "---", sprintf("%.3f", delta_rhat)),
      ESS = ifelse(is.na(delta_neff), "---", as.character(round(delta_neff))),
      `log BF10` = fmt_bf(log_bf10_nonlinear),
      check.names = FALSE
    )
    gp_delta_table$Rhat[masks$nonlinear == 0L] <- "---"
    gp_delta_table$ESS[masks$nonlinear == 0L] <- "---"
    print(gp_delta_table, row.names = FALSE, right = FALSE)
    cat("\n")

    # --- Lambda (lengthscale) ---
    cat("GP Lengthscale (lambda):\n")
    cat(strrep("-", 60), "\n")

    lambda_draws <- extract_gp_draws_full(result, post = post, param = "lambda")
    lambda_row_pattern <- if (identical(result$model_type, "full_select")) "^lambda_full(\\[|$)" else "^lambda(\\[|$)"
    lambda_rows <- grep(lambda_row_pattern, rownames(summary_df))
    lambda_neff <- if (length(lambda_rows) > 0) summary_df[lambda_rows, "n_eff"] else rep(NA, P)
    lambda_rhat <- if (length(lambda_rows) > 0) summary_df[lambda_rows, "Rhat"] else rep(NA, P)

    l_est <- colMeans(lambda_draws)
    l_se  <- apply(lambda_draws, 2, sd)
    l_lo  <- apply(lambda_draws, 2, quantile, 0.025)
    l_hi  <- apply(lambda_draws, 2, quantile, 0.975)

    gp_lambda_table <- data.frame(
      Predictor = pred_names,
      Estimate = round(l_est, digits),
      Est.Error = round(l_se, digits),
      `l-95% CI` = round(l_lo, digits),
      `u-95% CI` = round(l_hi, digits),
      Rhat = sprintf("%.3f", lambda_rhat),
      ESS = as.character(round(lambda_neff)),
      check.names = FALSE
    )
    gp_lambda_table$Rhat[masks$nonlinear == 0L] <- "---"
    gp_lambda_table$ESS[masks$nonlinear == 0L] <- "---"
    print(gp_lambda_table, row.names = FALSE, right = FALSE)
    cat("\n")
  }

  # ---- Residual (sigma) ----
  cat("Residual:\n")
  cat(strrep("-", 60), "\n")

  sigma_draws <- post$sigma
  sigma_row <- grep("^sigma$", rownames(summary_df))

  sigma_table <- data.frame(
    Parameter = "sigma",
    Estimate = round(mean(sigma_draws), digits),
    Est.Error = round(sd(sigma_draws), digits),
    `l-95% CI` = round(quantile(sigma_draws, 0.025), digits),
    `u-95% CI` = round(quantile(sigma_draws, 0.975), digits),
    Rhat = sprintf("%.3f", summary_df[sigma_row, "Rhat"]),
    ESS = as.character(round(summary_df[sigma_row, "n_eff"])),
    check.names = FALSE
  )
  print(sigma_table, row.names = FALSE, right = FALSE)
  cat("\n")

  # ---- Footer ----
  if (compute_bf && is_gp_full_fit(result)) {
    cat("log BF10: log Bayes factor in favour of the effect (linear or nonlinear).\n")
    cat("  Nonlinear BF uses interval-null with eps =", eps,
        "and moment-prior order d =", result$d_order, ".\n")
  }
  cat(strrep("=", 60), "\n")

  invisible(list(
    linear = linear_table,
    gp_delta = gp_delta_table,
    gp_lambda = gp_lambda_table,
    residual = sigma_table
  ))
}


# ============================================================
# Plotting: Posterior Distributions
# ============================================================

#' Common ggplot theme for posterior plots
.posterior_theme <- function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11, face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "black", fill = NA, linewidth = 0.5
      )
    )
}


#' Plot posterior distribution of delta (nonlinear magnitude)
#'
#' @param result A gp_fit object from fit_model()
#' @param xlim Custom x-axis limits
#' @param ylim Custom y-axis limits
#' @param bins Number of histogram bins
#' @param plot_prior If TRUE, overlay the moment prior density
#' @param title Optional custom title
#' @return A ggplot object
plot_posterior_delta <- function(
  result,
  xlim = NULL,
  ylim = NULL,
  bins = 30,
  plot_prior = FALSE,
  title = NULL
) {
  stopifnot(inherits(result, "gp_fit"), is_gp_full_fit(result))

  post <- rstan::extract(result$fit)
  draws <- extract_gp_draws_full(result, post = post, param = "delta")
  pred_names <- result$predictor_names

  if (is.null(title)) title <- paste0("Posterior of delta (d_order = ", result$d_order, ")")

  draws_df <- data.frame(
    value = as.vector(draws),
    predictor = rep(pred_names, each = nrow(draws))
  )
  draws_df$predictor <- factor(draws_df$predictor, levels = pred_names)

  p <- ggplot2::ggplot(draws_df, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = bins, fill = "black", alpha = 0.6, color = "white"
    ) +
    ggplot2::facet_wrap(~predictor, scales = "free", ncol = 3) +
    ggplot2::labs(title = title, x = "delta", y = "Density") +
    .posterior_theme()

  if (plot_prior) {
    kap <- get_kappa(result$d_order)
    prior_df <- do.call(rbind, lapply(pred_names, function(pn) {
      x_range <- range(draws_df$value[draws_df$predictor == pn], na.rm = TRUE)
      x_seq <- seq(0, x_range[2], length.out = 200)
      data.frame(
        x = x_seq,
        y = prior_moment_density(x_seq, d = result$d_order, kappa = kap),
        predictor = factor(pn, levels = pred_names)
      )
    }))
    p <- p +
      ggplot2::geom_line(
        data = prior_df, ggplot2::aes(x = x, y = y),
        color = "red", linewidth = 0.8
      )
  }

  if (!is.null(xlim) || !is.null(ylim))
    p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)

  p
}


#' Plot posterior distribution of theta (standardized linear coefficients)
#'
#' @param result A gp_fit object from fit_model()
#' @param xlim Custom x-axis limits
#' @param ylim Custom y-axis limits
#' @param bins Number of histogram bins
#' @param plot_prior If TRUE, overlay the implied prior density under g ~ IG(2, N / 2)
#' @param title Optional custom title
#' @return A ggplot object
plot_posterior_theta <- function(
  result,
  xlim = NULL,
  ylim = NULL,
  bins = 30,
  plot_prior = FALSE,
  title = NULL
) {
  stopifnot(inherits(result, "gp_fit"))

  post <- rstan::extract(result$fit)
  draws <- extract_gp_draws_full(result, post = post, param = "theta")
  pred_names <- result$predictor_names

  if (is.null(title)) title <- "Posterior of theta"

  draws_df <- data.frame(
    value = as.vector(draws),
    predictor = rep(pred_names, each = nrow(draws))
  )
  draws_df$predictor <- factor(draws_df$predictor, levels = pred_names)

  p <- ggplot2::ggplot(draws_df, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = bins, fill = "black", alpha = 0.6, color = "white"
    ) +
    ggplot2::facet_wrap(~predictor, scales = "free", ncol = 3) +
    ggplot2::labs(title = title, x = "theta", y = "Density") +
    .posterior_theme()

  if (plot_prior) {
    X <- result$stan_data$X
    N <- nrow(X)
    X_c <- scale(X, center = TRUE, scale = FALSE)
    XtX_inv <- solve(crossprod(X_c))
    v_theta <- diag(XtX_inv)
    ig_a <- 2
    ig_b <- 0.5 * N
    scales_theta <- sqrt((ig_b * v_theta) / ig_a)

    prior_df <- do.call(rbind, lapply(seq_along(pred_names), function(j) {
      pn <- pred_names[j]
      x_range <- range(draws_df$value[draws_df$predictor == pn], na.rm = TRUE)
      x_seq <- seq(min(x_range[1], 0), max(x_range[2], 0), length.out = 200)
      data.frame(
        x = x_seq,
        y = stats::dt(x_seq / scales_theta[j], df = 2 * ig_a) / scales_theta[j],
        predictor = factor(pn, levels = pred_names)
      )
    }))
    p <- p +
      ggplot2::geom_line(
        data = prior_df, ggplot2::aes(x = x, y = y),
        color = "red", linewidth = 0.8
      ) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4)
  }

  if (!is.null(xlim) || !is.null(ylim))
    p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)

  p
}


# ============================================================
# Plotting: GP Trend Functions
# ============================================================

#' Plot estimated GP trends (exact GP conditional posterior)
#'
#' For each selected predictor, computes the posterior mean and 95% CI
#' of the GP function f_j(x_j) + linear_j, based on the conditional
#' posterior of the GP given the fitted hyperparameters.
#'
#' @param result A gp_fit object from fit_model(model = "full")
#' @param p Optional vector of predictor indices to plot (default: all)
#' @param N_grid Number of grid points for prediction
#' @param ndraws Number of posterior draws to use
#' @param seed Random seed for draw subsampling
#' @param title Optional plot title
#' @param ylim Optional y-axis limits
#' @return A ggplot object
plot_gp_trends <- function(
  result,
  p = NULL,
  N_grid = 80,
  ndraws = 400,
  seed = 1,
  title = NULL,
  ylim = NULL,
  x_label = NULL,
  y_label = "y",
  predictor_labels = NULL,
  show_legend = NULL
) {
  stopifnot(inherits(result, "gp_fit"), is_gp_full_fit(result))
  stopifnot(!is.null(result$orig$X))

  fit <- result$fit
  X_s <- as.matrix(result$stan_data$X)
  y <- as.numeric(result$orig$y)
  X_orig <- as.matrix(result$orig$X)

  N <- nrow(X_s)
  P <- ncol(X_s)
  pred_names <- result$predictor_names
  masks <- get_gp_fit_masks(result)
  active_nl <- which(masks$nonlinear == 1L)

  if (is.null(p)) p_idx <- seq_len(P) else p_idx <- as.integer(p)

  # ---- Centering (mirrors Stan transformed data) ----
  y_bar <- mean(y)
  y_c <- y - y_bar
  x_bar_s <- colMeans(X_s)
  X_c <- sweep(X_s, 2, x_bar_s, "-")

  # ---- Kernel helpers ----
  se_cov_from_D2 <- function(D2, alpha, l) {
    K <- (alpha^2) * exp(-0.5 * D2 / (l^2))
    diag(K) <- diag(K) + 1e-12
    K
  }
  ktilde_from_K <- function(K, Q) {
    U <- K %*% Q
    M <- t(Q) %*% U
    Kt <- K - U %*% t(Q) - Q %*% t(U) + Q %*% M %*% t(Q)
    0.5 * (Kt + t(Kt))
  }
  ktilde_cross <- function(Kc, Qn, Qt) {
    KcPt <- Kc %*% Qt %*% t(Qt)
    PnKc <- Qn %*% (t(Qn) %*% Kc)
    PnKcPt <- Qn %*% (t(Qn) %*% KcPt)
    Kc - PnKc - KcPt + PnKcPt
  }
  chol_with_jitter <- function(Sigma, jitter_seq = c(0, 1e-10, 1e-8, 1e-6, 1e-4)) {
    Sigma <- 0.5 * (Sigma + t(Sigma))
    for (jit in jitter_seq) {
      Sigma_try <- Sigma
      if (jit > 0) diag(Sigma_try) <- diag(Sigma_try) + jit
      U <- try(chol(Sigma_try), silent = TRUE)
      if (!inherits(U, "try-error")) {
        return(U)
      }
    }
    NULL
  }
  projector_basis <- function(Z) {
    qz <- qr(Z)
    qr.Q(qz, complete = FALSE)[, seq_len(qz$rank), drop = FALSE]
  }

  # ---- Precompute training projector + D2 ----
  Q_tr <- projector_basis(cbind(1, X_c))
  D2_tr <- vector("list", P)
  for (j in seq_len(P)) {
    xj <- X_c[, j]
    D2 <- outer(xj, xj, "-")^2
    diag(D2) <- 0
    D2_tr[[j]] <- D2
  }

  # ---- Posterior draws ----
  post <- rstan::extract(fit)
  S_all <- length(post$sigma)
  set.seed(seed)
  idx <- if (ndraws < S_all) sample.int(S_all, ndraws) else seq_len(S_all)
  nd <- length(idx)

  sigma_draw  <- post$sigma[idx]
  beta1_draw  <- extract_gp_draws_full(result, post = post, param = "beta1")[idx, , drop = FALSE]
  delta_draw  <- extract_gp_draws_full(result, post = post, param = "delta")[idx, , drop = FALSE]
  lambda_draw <- extract_gp_draws_full(result, post = post, param = "lambda")[idx, , drop = FALSE]

  beta1_orig_draw <- result$coef_orig$beta1_draws[idx, , drop = FALSE]
  x_bar_orig <- colMeans(X_orig)

  # Scaling: x_s = (x_orig - min) / range - 0.5
  x_min <- apply(X_orig, 2, min)
  x_rng <- apply(X_orig, 2, max) - x_min

  # ---- Compute conditional posterior for each predictor ----
  out_list <- raw_list <- vector("list", length(p_idx))

  for (kk in seq_along(p_idx)) {
    i <- p_idx[kk]
    x_name <- pred_names[i]

    # Grid on original scale -> scaled -> centered
    xg_orig <- seq(min(X_orig[, i]), max(X_orig[, i]), length.out = N_grid)
    xg_s <- (xg_orig - x_min[i]) / x_rng[i] - 0.5
    xg_c <- xg_s - x_bar_s[i]

    Xg_c <- matrix(0, nrow = N_grid, ncol = P)
    Xg_c[, i] <- xg_c
    Qn <- projector_basis(cbind(1, Xg_c))
    D2_cross <- outer(xg_c, X_c[, i], "-")^2

    full_mat <- lin_mat <- matrix(NA_real_, nd, N_grid)
    ok_draw <- rep(FALSE, nd)

    for (s in seq_len(nd)) {
      sig <- sigma_draw[s]
      beta1 <- beta1_draw[s, ]

      if (!is.finite(sig) || sig <= 0) {
        next
      }

      # Build full Sigma
      Sigma <- matrix(0, N, N)
      for (j in active_nl) {
        alpha_j <- sig * delta_draw[s, j]
        lam_j <- lambda_draw[s, j]
        if (!is.finite(alpha_j) || !is.finite(lam_j) || alpha_j < 0 || lam_j <= 0) {
          next
        }
        Kj <- se_cov_from_D2(D2_tr[[j]], alpha_j, lam_j)
        Sigma <- Sigma + ktilde_from_K(Kj, Q_tr)
      }
      diag(Sigma) <- diag(Sigma) + sig^2

      resid <- y_c - as.numeric(X_c %*% beta1)
      U <- chol_with_jitter(Sigma)
      if (is.null(U)) {
        next
      }
      a_vec <- backsolve(U, forwardsolve(t(U), resid))

      # Cross-covariance for predictor i
      if (masks$nonlinear[i] == 1L && delta_draw[s, i] > 0 && lambda_draw[s, i] > 0) {
        alpha_i <- sig * delta_draw[s, i]
        lam_i <- lambda_draw[s, i]
        Kc <- (alpha_i^2) * exp(-0.5 * D2_cross / (lam_i^2))
        Kc_tilde <- ktilde_cross(Kc, Qn, Q_tr)
        f_mean <- as.numeric(Kc_tilde %*% a_vec)
      } else {
        f_mean <- rep(0, N_grid)
      }

      lin <- y_bar + beta1_orig_draw[s, i] * (xg_orig - x_bar_orig[i])
      full_mat[s, ] <- lin + f_mean
      lin_mat[s, ] <- lin
      ok_draw[s] <- TRUE
    }

    if (!any(ok_draw)) {
      stop(
        "plot_gp_trends() could not obtain a numerically stable conditional ",
        "posterior for predictor '", x_name, "'."
      )
    }

    out_list[[kk]] <- data.frame(
      predictor = x_name,
      x = xg_orig,
      full_mean = colMeans(full_mat[ok_draw, , drop = FALSE]),
      full_lo = apply(full_mat[ok_draw, , drop = FALSE], 2, quantile, 0.025),
      full_hi = apply(full_mat[ok_draw, , drop = FALSE], 2, quantile, 0.975),
      lin_mean = colMeans(lin_mat[ok_draw, , drop = FALSE])
    )
    raw_list[[kk]] <- data.frame(predictor = x_name, x = X_orig[, i], y = y)
  }

  df <- do.call(rbind, out_list)
  raw <- do.call(rbind, raw_list)
  df$predictor <- factor(df$predictor, levels = pred_names[p_idx])
  raw$predictor <- factor(raw$predictor, levels = pred_names[p_idx])

  if (is.null(title)) title <- NULL
  if (is.null(x_label)) x_label <- if (length(p_idx) == 1L) pred_names[p_idx] else NULL
  if (is.null(show_legend)) show_legend <- length(p_idx) == 1L

  predictor_map <- setNames(pred_names[p_idx], pred_names[p_idx])
  if (!is.null(predictor_labels)) {
    predictor_map[names(predictor_labels)] <- predictor_labels
  }

  df$predictor_label <- factor(
    predictor_map[as.character(df$predictor)],
    levels = unname(predictor_map)
  )
  raw$predictor_label <- factor(
    predictor_map[as.character(raw$predictor)],
    levels = unname(predictor_map)
  )

  line_df <- rbind(
    data.frame(
      predictor_label = df$predictor_label,
      x = df$x,
      mean = df$lin_mean,
      curve = "Linear mean"
    ),
    data.frame(
      predictor_label = df$predictor_label,
      x = df$x,
      mean = df$full_mean,
      curve = "Linear + GP mean"
    )
  )

  plt <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = raw,
      ggplot2::aes(x = x, y = y),
      color = "grey20",
      alpha = 0.85,
      size = 1.7
    ) +
    ggplot2::geom_ribbon(
      data = df,
      ggplot2::aes(x = x, ymin = full_lo, ymax = full_hi),
      fill = "#7AA6C2",
      alpha = 0.22
    ) +
    ggplot2::geom_line(
      data = line_df,
      ggplot2::aes(x = x, y = mean, color = curve, linetype = curve),
      linewidth = 0.9
    ) +
    ggplot2::scale_color_manual(
      values = c("Linear mean" = "grey20", "Linear + GP mean" = "#0B5A8F"),
      breaks = c("Linear mean", "Linear + GP mean")
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Linear mean" = "22", "Linear + GP mean" = "solid"),
      breaks = c("Linear mean", "Linear + GP mean")
    ) +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = y_label,
      color = NULL,
      linetype = NULL
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0),
      axis.title = ggplot2::element_text(size = 11),
      axis.text = ggplot2::element_text(size = 10, color = "black"),
      legend.position = if (isTRUE(show_legend)) "top" else "none",
      legend.box.just = "left",
      legend.text = ggplot2::element_text(size = 10),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      panel.border = ggplot2::element_blank()
    )

  if (length(p_idx) > 1L) {
    plt <- plt + ggplot2::facet_wrap(~predictor_label, scales = "free_x")
  }

  if (!is.null(ylim)) plt <- plt + ggplot2::coord_cartesian(ylim = ylim)

  plt
}
