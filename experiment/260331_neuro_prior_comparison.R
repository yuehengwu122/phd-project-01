# ============================================================
# Experiment: Neuroscience — comparing all four prior orders
# Dataset:    Lock5Data::FacebookFriends
#   Predictor: sqrt(FBfriends)
#   Response:  GMdensity (gray matter density)
#
# Steps:
#   1. Fit full GP (d_order 0–3) + linear null
#   2. Capture estimation summaries → theta / delta tables
#   3. Bridge sampling BF (full vs linear)
#   4. Interval-null BFs (logspline / ecdf × eps 0 / 0.2)
#   5–6. Posterior and trend plots → obsidian/figures/
#   7. Write all results into the Obsidian markdown in-place
#
# Output:
#   results/fit_01_neuro_{d0,d1,d2,d3,linear}.rds
#   obsidian/figures/260331_*.png
#   obsidian/02 Experiment Log/260331 neuro prior comparison.md
# ============================================================

source("src/_setup.R")

library(Lock5Data)
library(ggplot2)
library(tidyr)

data(FacebookFriends)
X <- matrix(sqrt(FacebookFriends$FBfriends), ncol = 1)
colnames(X) <- "sqrt_FBfriends"
y <- FacebookFriends$GMdensity


# ============================================================
# Step 1: Fit models
# ============================================================

fit_d0  <- fit_model(X, y, model = "full",   d_order = 0, chains = 2, iter = 5500, warmup = 500)
fit_d1  <- fit_model(X, y, model = "full",   d_order = 1, chains = 2, iter = 5500, warmup = 500)
fit_d2  <- fit_model(X, y, model = "full",   d_order = 2, chains = 2, iter = 5500, warmup = 500)
fit_d3  <- fit_model(X, y, model = "full",   d_order = 3, chains = 2, iter = 5500, warmup = 500)
fit_lin <- fit_model(X, y, model = "linear",              chains = 2, iter = 5500, warmup = 500)

saveRDS(fit_d0,  "results/fit_01_neuro_d0.rds")
saveRDS(fit_d1,  "results/fit_01_neuro_d1.rds")
saveRDS(fit_d2,  "results/fit_01_neuro_d2.rds")
saveRDS(fit_d3,  "results/fit_01_neuro_d3.rds")
saveRDS(fit_lin, "results/fit_01_neuro_linear.rds")

fits_full <- list(d0 = fit_d0, d1 = fit_d1, d2 = fit_d2, d3 = fit_d3)


# ============================================================
# Step 2: Capture estimation summaries
#
# summary_gp_fit() prints to console (useful progress info) and
# returns a list invisibly. We capture the structured return value
# to build cross-prior tables for the markdown.
# ============================================================

summaries <- lapply(fits_full, summary_gp_fit)

# Extract theta (linear slope) row for the single predictor
theta_summary_df <- do.call(rbind, lapply(names(summaries), function(nm) {
  lin <- summaries[[nm]]$linear
  row <- lin[lin$Predictor == "sqrt_FBfriends", ]
  data.frame(
    Prior          = nm,
    `theta mean`   = row$Estimate,
    `l-95% CI`     = row[["l-95% CI"]],
    `u-95% CI`     = row[["u-95% CI"]],
    `log BF10`     = row[["log BF10"]],   # character, already formatted
    ESS            = row$ESS,             # character
    Rhat           = row$Rhat,            # character
    check.names    = FALSE
  )
}))

# Extract delta (nonlinear magnitude) row for the single predictor
delta_summary_df <- do.call(rbind, lapply(names(summaries), function(nm) {
  dlt <- summaries[[nm]]$gp_delta
  row <- dlt[dlt$Predictor == "sqrt_FBfriends", ]
  data.frame(
    Prior          = nm,
    `delta mean`   = row$Estimate,
    `l-95% CI`     = row[["l-95% CI"]],
    `u-95% CI`     = row[["u-95% CI"]],
    `log BF10`     = row[["log BF10"]],   # character
    ESS            = row$ESS,             # character
    Rhat           = row$Rhat,            # character
    check.names    = FALSE
  )
}))


# ============================================================
# Step 3: Bridge sampling BF (full GP vs linear null)
# ============================================================

bf_bridge_df <- do.call(rbind, lapply(names(fits_full), function(nm) {
  logbf <- compute_bf_bridge(fits_full[[nm]], fit_lin, log = TRUE)
  data.frame(prior = nm, log_BF10_bridge = logbf)
}))


# ============================================================
# Step 4: Interval-null BFs
#
# compute_bf_table() uses logspline internally and has no method
# argument. For the ecdf method we call bf10_interval_null_delta()
# directly via a thin wrapper.
# ============================================================

compute_bf_table_method <- function(result, eps, method = "logspline") {
  post       <- rstan::extract(result$fit)
  X_s        <- result$stan_data$X
  P          <- result$stan_data$P
  d_order    <- result$d_order
  pred_names <- result$predictor_names

  # Linear BF: Savage-Dickey at theta = 0 (always logspline)
  prior0_theta <- prior0_theta_analytical(X_s)
  theta_draws  <- post$theta
  post0_theta  <- vapply(seq_len(P), function(j)
    post0_theta_logspline(theta_draws[, j]), numeric(1))
  log_bf10_linear <- log(prior0_theta / post0_theta)

  # Nonlinear BF: interval-null for delta with specified method
  delta_draws <- post$delta
  log_bf10_nonlinear <- vapply(seq_len(P), function(j) {
    draws_j <- if (is.matrix(delta_draws)) delta_draws[, j] else delta_draws
    bf <- bf10_interval_null_delta(draws_j, d_order,
                                   eps = eps, method = method,
                                   return_log = TRUE)
    bf$log_bf10[1]
  }, numeric(1))

  data.frame(
    predictor          = pred_names,
    log_bf10_linear    = round(log_bf10_linear, 4),
    log_bf10_nonlinear = round(log_bf10_nonlinear, 4)
  )
}

eps_grid <- list(
  list(eps = 0,   method = "logspline"),
  list(eps = 0.2, method = "logspline"),
  list(eps = 0.2, method = "ecdf")       # "ecdf" = "raw" method
)

bf_interval_df <- do.call(rbind, unlist(recursive = FALSE,
  lapply(names(fits_full), function(nm) {
    lapply(eps_grid, function(eg) {
      tbl        <- compute_bf_table_method(fits_full[[nm]], eg$eps, eg$method)
      tbl$prior  <- nm
      tbl$eps    <- eg$eps
      tbl$method <- eg$method
      tbl
    })
  })
))


# ============================================================
# Step 4b: BF comparison plot (bridge vs logspline vs ECDF, eps = 0.2)
# ============================================================

EPS <- 0.2

bf_plot_df <- do.call(rbind, lapply(names(fits_full), function(nm) {
  result  <- fits_full[[nm]]
  d_order <- result$d_order
  draws   <- rstan::extract(result$fit)$delta

  log_bf10_bridge    <- compute_bf_bridge(result, fit_lin, log = TRUE)
  log_bf10_logspline <- bf10_interval_null_delta(draws, d_order, eps = EPS,
                                                 method = "logspline",
                                                 return_log = TRUE)$log_bf10
  log_bf10_ecdf      <- bf10_interval_null_delta(draws, d_order, eps = EPS,
                                                 method = "ecdf",
                                                 return_log = TRUE)$log_bf10

  data.frame(
    order              = factor(nm, levels = c("d0","d1","d2","d3")),
    log_bf10_bridge    = log_bf10_bridge,
    log_bf10_logspline = log_bf10_logspline,
    log_bf10_ecdf      = log_bf10_ecdf
  )
}))

bf_plot_long <- tidyr::pivot_longer(
  bf_plot_df,
  cols      = c(log_bf10_bridge, log_bf10_logspline, log_bf10_ecdf),
  names_to  = "method",
  values_to = "log_bf10"
)
bf_plot_long$method <- factor(
  bf_plot_long$method,
  levels = c("log_bf10_bridge", "log_bf10_logspline", "log_bf10_ecdf"),
  labels = c("Bridge", "Logspline (eps=0.2)", "ECDF (eps=0.2)")
)

p_bf_comparison <- ggplot(
  bf_plot_long,
  aes(x = order, y = log_bf10,
      color = method, shape = method, linetype = method, group = method)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2.5) +
  scale_color_manual(
    values = c("Bridge"              = "#1b7837",
               "Logspline (eps=0.2)" = "#2166ac",
               "ECDF (eps=0.2)"      = "#d6604d")
  ) +
  scale_shape_manual(
    values = c("Bridge" = 16, "Logspline (eps=0.2)" = 17, "ECDF (eps=0.2)" = 15)
  ) +
  scale_linetype_manual(
    values = c("Bridge" = "solid", "Logspline (eps=0.2)" = "dashed", "ECDF (eps=0.2)" = "dotted")
  ) +
  labs(
    title  = "Neuroscience: log BF10 by method and prior order",
    x      = "Prior order",
    y      = "log BF10",
    color  = "Method", shape = "Method", linetype = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )




# ============================================================
# Steps 5 & 6: Save all figures to obsidian/figures/
# ============================================================

fig_dir <- "obsidian/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

for (nm in names(fits_full)) {
  p_delta <- plot_posterior_delta(fits_full[[nm]], bins = 50, plot_prior = TRUE, 
                                  xlim = c(0, 25), ylim = c(0, 0.6))
  p_theta <- plot_posterior_theta(fits_full[[nm]], bins = 50, plot_prior = TRUE)
  p_trend <- plot_gp_trends(fits_full[[nm]])

  ggsave(file.path(fig_dir, paste0("260331_posterior_delta_", nm, ".png")),
         p_delta, width = 6, height = 4, dpi = 150)
  ggsave(file.path(fig_dir, paste0("260331_posterior_theta_", nm, ".png")),
         p_theta, width = 6, height = 4, dpi = 150)
  ggsave(file.path(fig_dir, paste0("260331_gp_trend_", nm, ".png")),
         p_trend, width = 7, height = 4, dpi = 150)
}

ggsave(file.path(fig_dir, "260331_bf_comparison.png"), p_bf_comparison, width = 6, height = 4, dpi = 150)


# ============================================================
# Step 7: Write all results into the Obsidian markdown in-place
#
# replace_md_table_rows() is anchored on "^\\s*\\|" so it only
# ever matches table header rows — never body prose — making it
# robust even when the pattern word also appears in paragraphs.
# ============================================================

replace_md_table_rows <- function(md, header_pattern, new_data_rows) {
  # Match only lines that begin with a pipe (table rows)
  h <- grep(paste0("^\\s*\\|.*", header_pattern), md, perl = TRUE)[1]
  if (is.na(h)) {
    warning("Markdown table header not found for pattern: ", header_pattern)
    return(md)
  }
  sep <- h + 1   # |---|---| separator line
  end <- sep     # walk forward over all existing data rows
  while (end < length(md) && grepl("^\\s*\\|", md[end + 1])) end <- end + 1
  c(
    md[seq_len(sep)],
    new_data_rows,
    if (end < length(md)) md[seq(end + 1, length(md))] else character(0)
  )
}

md_path <- "obsidian/02 Experiment Log/260331 neuro prior comparison.md"
md      <- readLines(md_path, warn = FALSE)

# --- Theta summary ---
theta_rows <- sprintf("| %-5s | %-10s | %-8s | %-8s | %-9s | %-5s | %-5s |",
  theta_summary_df$Prior,
  sprintf("%.4f", as.numeric(theta_summary_df$`theta mean`)),
  sprintf("%.4f", as.numeric(theta_summary_df$`l-95% CI`)),
  sprintf("%.4f", as.numeric(theta_summary_df$`u-95% CI`)),
  theta_summary_df$`log BF10`,
  theta_summary_df$ESS,
  theta_summary_df$Rhat)
md <- replace_md_table_rows(md, "theta mean", theta_rows)

# --- Delta summary ---
delta_rows <- sprintf("| %-5s | %-10s | %-8s | %-8s | %-9s | %-5s | %-5s |",
  delta_summary_df$Prior,
  sprintf("%.4f", as.numeric(delta_summary_df$`delta mean`)),
  sprintf("%.4f", as.numeric(delta_summary_df$`l-95% CI`)),
  sprintf("%.4f", as.numeric(delta_summary_df$`u-95% CI`)),
  delta_summary_df$`log BF10`,
  delta_summary_df$ESS,
  delta_summary_df$Rhat)
md <- replace_md_table_rows(md, "delta mean", delta_rows)

# --- Bridge BF ---
bridge_rows <- sprintf("| %-5s | %-19s |",
  bf_bridge_df$prior,
  sprintf("%.4f", bf_bridge_df$log_BF10_bridge))
md <- replace_md_table_rows(md, "bridge", bridge_rows)

# --- Interval-null BF ---
fmt_eps <- function(x) ifelse(x == 0, "0  ", sprintf("%.1f", x))
interval_rows <- sprintf("| %-5s | %-3s | %-9s | %-18s | %-20s |",
  bf_interval_df$prior,
  fmt_eps(bf_interval_df$eps),
  bf_interval_df$method,
  sprintf("%.4f", bf_interval_df$log_bf10_linear),
  sprintf("%.4f", bf_interval_df$log_bf10_nonlinear))
md <- replace_md_table_rows(md, "method", interval_rows)

writeLines(md, md_path)
