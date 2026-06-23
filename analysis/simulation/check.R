# Quick diagnostic: print exact values behind Plot 6
# Run after the main plotting script has loaded and processed bf_bridge

df6_check <- bf_bridge |>
  dplyr::group_by(delta_true, N, order) |>
  dplyr::summarise(
    median_logBF = median(log_bf10_bridge, na.rm = TRUE),
    q25 = quantile(log_bf10_bridge, 0.25, na.rm = TRUE),
    q75 = quantile(log_bf10_bridge, 0.75, na.rm = TRUE),
    n_reps = sum(!is.na(log_bf10_bridge)),
    .groups = "drop"
  ) |>
  dplyr::arrange(delta_true, N, order)

# Print full table — no rounding, so we can see if values are identical or just close
options(width = 120, pillar.sigfig = 6)
print(df6_check, n = Inf)

# Also check: are they literally identical within each (delta, N) group?
cat("\n\n=== Max spread (max - min) of median_logBF within each (delta, N) ===\n")
spread <- df6_check |>
  dplyr::group_by(delta_true, N) |>
  dplyr::summarise(
    max_diff = max(median_logBF) - min(median_logBF),
    values = paste(sprintf("%.4f", median_logBF), collapse = " | "),
    .groups = "drop"
  )
print(spread, n = Inf)
