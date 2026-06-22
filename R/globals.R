# Silence R CMD check NSE notes for non-standard evaluation column names.
utils::globalVariables(c(
  "a_pred", "angle", "b_pred", "boot_index", "cluster", "count", "cx", "cy", "d2",
  "eg", "hessian_ok", "is_noise", "lam", "mean_U", "mean_x", "mean_y", "n_mins",
  "n_runs", "s_xx", "s_xy", "s_yy", "sd_x", "Sigma", "V", "var_x",
  "x_conf_lower", "x_conf_upper", "x_pred_lower", "x_pred_upper", "U"
))
