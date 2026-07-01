##########
# Declare variables used in non-standard evaluation (dplyr/ggplot2 column
# references, the magrittr `.` pronoun, etc.) so that `R CMD check` does not
# report them as "no visible binding for global variable". These are not real
# objects in the package namespace; they are column names evaluated lazily
# inside data-masking verbs.
##########

utils::globalVariables(c(
  ".",
  ".data",
  "RSS_two",
  "algorithm",
  "algorithm_vt1",
  "algorithm_vt2",
  "bin",
  "bp",
  "determinant_bp",
  "dist_MSE_ratio",
  "dist_x_sq",
  "dist_y_sq",
  "est_ci",
  "f_stat",
  "hr",
  "idx",
  "inside_ci",
  "int_point_x",
  "inv_weight",
  "key",
  "method",
  "normalized_dist",
  "outlier",
  "p",
  "p_val_f",
  "pct_tot_weight",
  "pf_two",
  "pos_change",
  "pos_slope_after_bp",
  "sum_sq",
  "time",
  "vco2",
  "ve",
  "ve_vco2",
  "ve_vo2",
  "vo2",
  "x_var",
  "x_vt1",
  "x_vt2",
  "y_hat_left",
  "y_hat_right",
  "y_var",
  "y_vt1",
  "y_vt2"
))
