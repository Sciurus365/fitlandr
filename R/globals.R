# Silence R CMD check NSE notes for non-standard evaluation column names.
utils::globalVariables(c(
  "a_pred", "angle", "b_pred", "boot_index", "cluster", "cx", "cy", "d2",
  "eg", "is_noise", "lam", "mean_x", "mean_y", "n_runs", "s_xx", "s_xy",
  "s_yy", "Sigma", "V"
))
