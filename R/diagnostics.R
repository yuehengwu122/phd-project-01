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
print_quick_summary <- function(result, key_params = c("delta", "lambda", "sigma", "beta0", "beta1")) {
  
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
}


# ============================================================
# Plot Functions
# ============================================================

#' Plot Posterior Distribution
#' @param result Result object from run_analysis()
#' @param param Parameter name to plot (default: "delta")
#' @param xlim Custom x-axis limits (default: NULL, auto)
#' @param predictor_names Optional names for predictors
#' @param bins Number of histogram bins (default: 30)
#' @param trim_quantile Quantile threshold for trimming heavy tails (default: 1, no trimming)
#' @param model_name Optional model name to display in title (default: NULL, uses result$model_name)
plot_posterior <- function(result, 
                          param = "delta", 
                          xlim = NULL,
                          predictor_names = NULL,
                          bins = 30,
                          trim_quantile = 1,
                          model_name = NULL) {
  
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
  
  # Handle vector vs scalar parameter
  if (is.matrix(draws)) {
    P <- ncol(draws)
    
    # Set predictor names
    if (is.null(predictor_names)) {
      predictor_names <- paste0(param, "[", 1:P, "]")
    }
    
    # Convert to long format
    draws_df <- data.frame(
      value = as.vector(draws),
      predictor = rep(predictor_names, each = length(draws[, 1]))
    )
    draws_df$predictor <- factor(draws_df$predictor, levels = predictor_names)
    
    p <- ggplot(draws_df, aes(x = value)) +
      geom_histogram(aes(y = after_stat(density)), bins = bins, fill = "steelblue", alpha = 0.6, color = "black") +
      facet_wrap(~ predictor, scales = "free") +  # Each predictor gets its own axis
      labs(title = plot_title,
           x = param, y = "Density") +
      theme_minimal() +
      theme(plot.title = element_text(size = 11))
    
  } else {
    # Scalar parameter
    draws_df <- data.frame(value = draws)
    
    p <- ggplot(draws_df, aes(x = value)) +
      geom_histogram(aes(y = after_stat(density)), bins = bins, fill = "steelblue", alpha = 0.6, color = "black") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      labs(title = plot_title,
           x = param, y = "Density") +
      theme_minimal() +
      theme(plot.title = element_text(size = 11))
  }
  
  # Apply custom xlim if provided
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim)
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
#' @param result Result object from run_analysis()
#' @param X Predictor matrix (N x P), should match the data used for fitting
#' @param predictor_names Optional names for predictors
#' @param p Optional vector of predictor indices to plot (default: all)
#' @param N_grid Number of grid points for prediction (default: 50)
#' @param ndraws Number of posterior draws to use (default: 400)
#' @param seed Random seed for drawing samples (default: 1)
#' @param model_name Optional model name to display in title (default: NULL, uses result$model_name)
#' @param ylim Optional y-axis limits as a numeric vector of length 2 (default: NULL, auto)
plot_gp_trends <- function(result,
                           X = NULL,
                           predictor_names = NULL,
                           p = NULL,
                           N_grid = 50,
                           ndraws = 400,
                           seed = 1,
                           model_name = NULL,
                           ylim = NULL) {
  fit <- result$fit
  y   <- result$data$y
  if (is.null(X)) X <- result$data$X
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  X <- as.matrix(X)

  N <- nrow(X)
  P <- ncol(X)

  if (is.null(p)) p_idx <- seq_len(P) else p_idx <- as.integer(p)

  if (is.null(predictor_names)) {
    predictor_names <- colnames(X)
    if (is.null(predictor_names)) predictor_names <- paste0("X", 1:P)
  }

  # ---- centering like Stan ----
  y_bar <- mean(y)
  y_c   <- y - y_bar
  x_bar <- colMeans(X)
  X_c   <- sweep(X, 2, x_bar, "-")

  # ---- precompute Q_tr, R_tr for each predictor (span{1, x_cj}) ----
  Q_tr <- vector("list", P)  # each N x 2
  R_tr <- vector("list", P)  # each 2 x 2 (inverse of Q'Q)
  xjxj <- numeric(P)

  for (j in 1:P) {
    Q <- cbind(rep(1, N), X_c[, j])
    R <- solve(crossprod(Q))
    Q_tr[[j]] <- Q
    R_tr[[j]] <- R
    xjxj[j]   <- sum(X_c[, j]^2)
  }

  # ---- precompute D2 for each predictor (as in Stan transformed data) ----
  D2_list <- vector("list", P)
  for (j in 1:P) {
    xj <- X_c[, j]
    # squared distance matrix
    D2 <- outer(xj, xj, "-")
    D2 <- D2 * D2
    diag(D2) <- 0
    D2_list[[j]] <- D2
  }

  # ---- helper: SE kernel from D2 ----
  se_cov_from_D2_R <- function(D2, alpha, l) {
    a2 <- alpha^2
    inv_l2 <- 1 / (l^2)
    K <- a2 * exp(-0.5 * D2 * inv_l2)
    diag(K) <- diag(K) + 1e-12
    K
  }

  # ---- helper: Ktilde = K - PK - KP + PKP with P = Q R Q' ----
  ktilde_from_K <- function(K, Q, R) {
    # U = K Q
    U <- K %*% Q              # N x 2
    # PK = Q (R U')
    PK <- Q %*% (R %*% t(U))  # N x N
    # M = Q' K Q = Q' U
    M  <- t(Q) %*% U          # 2 x 2
    # PKP = Q (R M R) Q'
    PKP <- Q %*% (R %*% M %*% R) %*% t(Q)
    Kt <- K - PK - t(PK) + PKP
    0.5 * (Kt + t(Kt))
  }

  # ---- helper: orthogonalized cross-cov (grid x train) ----
  ktilde_cross <- function(Kc, Qn, Rn, Qt, Rt) {
    # Kc is G x N

    # Kc * Pt = (Kc*Qt) Rt Qt'
    T  <- Kc %*% Qt                   # G x 2
    KcPt <- (T %*% Rt) %*% t(Qt)      # G x N

    # Pn * Kc = Qn Rn (Qn' Kc)
    S  <- t(Qn) %*% Kc                # 2 x N
    PnKc <- Qn %*% (Rn %*% S)         # G x N

    # Pn * KcPt
    S2 <- t(Qn) %*% KcPt
    PnKcPt <- Qn %*% (Rn %*% S2)

    Kc - PnKc - KcPt + PnKcPt
  }

  # ---- extract posterior draws needed ----
  post <- rstan::extract(fit, pars = c("mu","sigma","theta","delta","lambda"))
  S <- length(post$mu)

  set.seed(seed)
  if (ndraws < S) idx <- sample.int(S, ndraws) else idx <- seq_len(S)
  nd <- length(idx)

  mu_draw     <- post$mu[idx]
  sigma_draw  <- post$sigma[idx]
  theta_draw  <- post$theta[idx, , drop = FALSE]
  delta_draw  <- post$delta[idx, , drop = FALSE]
  lambda_draw <- post$lambda[idx, , drop = FALSE]

  out_list <- vector("list", length(p_idx))
  raw_list <- vector("list", length(p_idx))

  # ---- plotting loop ----
  for (kk in seq_along(p_idx)) {
    i <- p_idx[kk]
    x_name <- predictor_names[i]
    
    # grid in raw scale
    x_i <- X[, i]
    xg_raw <- seq(min(x_i), max(x_i), length.out = N_grid)
    xg_c   <- xg_raw - x_bar[i]

    # grid-side Qn, Rn for projector span{1, xg_c}
    Qn <- cbind(rep(1, N_grid), xg_c)
    Rn <- solve(crossprod(Qn))

    # store per-draw curves
    full_mat <- matrix(NA_real_, nd, N_grid)
    lin_mat  <- matrix(NA_real_, nd, N_grid)

    # precompute grid->train squared distances for predictor i (depends only on X_c and xg_c)
    x_train_c <- X_c[, i]
    D2_cross <- outer(xg_c, x_train_c, "-")
    D2_cross <- D2_cross * D2_cross

    for (s in 1:nd) {
      # parameters
      sigma <- sigma_draw[s]
      mu    <- mu_draw[s]
      beta1 <- sigma * theta_draw[s, ]
      alpha <- sigma * delta_draw[s, ]  # alpha_j = delta_j * sigma
      lam   <- lambda_draw[s, ]

      # build Sigma = sum_j Ktilde_j + sigma^2 I
      Sigma <- matrix(0, N, N)
      for (j in 1:P) {
        Kj <- se_cov_from_D2_R(D2_list[[j]], alpha[j], lam[j])
        Q  <- Q_tr[[j]]
        R  <- R_tr[[j]]
        Sigma <- Sigma + ktilde_from_K(Kj, Q, R)
      }
      diag(Sigma) <- diag(Sigma) + sigma^2

      # resid = y_c - (mu + X_c beta1)
      mean_vec <- as.numeric(mu + X_c %*% beta1)
      resid <- y_c - mean_vec

      # a = Sigma^{-1} resid via Cholesky
      L <- chol(Sigma)
      v <- forwardsolve(t(L), resid)  # because chol() gives upper in R
      a_vec <- backsolve(L, v)

      # K_cross for predictor i (grid x train) using SE from D2_cross
      Kc <- (alpha[i]^2) * exp(-0.5 * D2_cross / (lam[i]^2))

      # orthogonalize cross-cov: (I-Pn) Kc (I-Pt)
      Qt <- Q_tr[[i]]
      Rt <- R_tr[[i]]
      Kc_tilde <- ktilde_cross(Kc, Qn, Rn, Qt, Rt)

      f_mean <- as.numeric(Kc_tilde %*% a_vec)  # grid vector

      beta0_orig <- mu + y_bar
      lin <- beta0_orig + beta1[i] * xg_c
      full <- lin + f_mean

      lin_mat[s, ]  <- lin
      full_mat[s, ] <- full
    }

    # summarize
    full_mean  <- colMeans(full_mat)
    full_lower <- apply(full_mat, 2, quantile, 0.025)
    full_upper <- apply(full_mat, 2, quantile, 0.975)
    lin_mean   <- colMeans(lin_mat)

    out_list[[kk]] <- data.frame(
      predictor = factor(x_name, levels = predictor_names[p_idx]),
      x = xg_raw,
      full_mean = full_mean,
      full_lo = full_lower,
      full_hi = full_upper,
      lin_mean = lin_mean
    )

    raw_list[[kk]] <- data.frame(
      predictor = factor(x_name, levels = predictor_names[p_idx]),
      x = X[, i],
      y = y
    )
  }

  df  <- do.call(rbind, out_list)
  raw <- do.call(rbind, raw_list)

  # Create plot title with model name if available
  if (is.null(model_name)) {
    model_name <- result$model_name
  }
  
  plot_title <- if (!is.null(model_name)) {
    paste("GP Trends:", model_name)
  } else {
    "GP Trends"
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
#' @param p Optional vector of predictor indices to plot (default: all)
#' @param N_grid Number of grid points for prediction (default: 50)
#' @param ndraws Number of posterior draws to use (default: 1000)
#' @param seed Random seed for drawing samples (default: 1)
#' @param predictor_names Optional names for predictors
#' @param model_name Optional model name to display in title (default: NULL, uses object$model_name)
#' @param ylim Optional y-axis limits as a numeric vector of length 2 (default: NULL, auto)
plot_hsgp_trends <- function(object, p = NULL, N_grid = 50, ndraws = 1000, seed = 1,
                             predictor_names = NULL, model_name = NULL, ylim = NULL) {

  fit <- object$fit
  X   <- as.matrix(object$data$X)
  y   <- object$data$y
  M   <- object$data$M

  N <- nrow(X); P <- ncol(X)

  if (is.null(p)) p_idx <- seq_len(P) else p_idx <- as.integer(p)

  if (is.null(predictor_names)) {
    predictor_names <- colnames(X)
    if (is.null(predictor_names)) predictor_names <- paste0("X", seq_len(P))
  }

  # centering like Stan
  y_bar <- mean(y)
  x_bar <- colMeans(X)
  X_c <- sweep(X, 2, x_bar, "-")

  post <- rstan::extract(fit, pars = c("mu","sigma","theta","delta","lambda","z_basis"))
  S <- length(post$mu)
  if (S > ndraws) {
    set.seed(seed)
    idx <- sample.int(S, ndraws)
  } else idx <- seq_len(S)

  mu_draw    <- post$mu[idx]
  sigma_draw <- post$sigma[idx]

  out_list <- vector("list", length(p_idx))
  raw_list <- vector("list", length(p_idx))

  for (kk in seq_along(p_idx)) {
    pp <- p_idx[kk]
    x_name <- predictor_names[pp]

    xcp <- X_c[, pp]
    half_range <- 0.5 * (max(xcp) - min(xcp))
    Lp <- max(1e-6, 1.5 * half_range)

    # grid
    x_grid_raw <- seq(min(X[, pp]), max(X[, pp]), length.out = N_grid)
    x_grid_c   <- x_grid_raw - x_bar[pp]
    PHI_g      <- phi_grid_ortho(x_grid_c, M, Lp)

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
      x = x_grid_raw,
      full_mean = colMeans(full_mat),
      full_lo = apply(full_mat, 2, quantile, 0.025),
      full_hi = apply(full_mat, 2, quantile, 0.975),
      lin_mean = colMeans(lin_mat)
    )

    raw_list[[kk]] <- data.frame(
      predictor = factor(x_name, levels = predictor_names[p_idx]),
      x = X[, pp],
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
