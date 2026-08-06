#' Second derivative inflection
#'
#' Find a breakpoint in ventilatory data based on the 2nd derivative inflection point after Santos et al. (2004). In their paper, they use pulmonary ventilation (VE) vs. VO2 to find the first ventilatory (VT1) / gas exchange threshold (GET) and VE vs. VCO2 to find the respiratory compensation point (RCP) / second ventilatory threshold (VT2). It is generally suitable to also use VE vs. time to find VT2. One can truncate the test and rerun VE vs. time to find VT1.
#'
#' @param .data Gas exchange data frame or tibble.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param bp Is this the first (\code{vt1}) or the second (\code{vt2}) breakpoint?
#' @param degree The polynomial degree fit to the curve. Default \code{5}, matching Santos and Giannella-Neto (2004). When \code{search_degree = TRUE}, this is the lowest degree searched. Must be at least 5.
#' @param search_degree Search for the best polynomial degree rather than using the fixed \code{degree}? Default \code{FALSE}.
#' @param max_degree Highest polynomial degree considered when \code{search_degree = TRUE}. Default \code{8}.
#' @param criterion Model-selection rule when \code{search_degree = TRUE}: \code{"aicc"} (default) or \code{"lrt"}.
#' @param vo2 Name of the \code{vo2} variable
#' @param vco2 Name of the \code{vco2} variable
#' @param ve Name of the \code{ve} variable
#' @param time Name of the \code{time} variable
#' @param alpha_linearity Significance value used by the \code{"lrt"} degree search.
#' @param ... Dot dot dot mostly allows this function to work properly if breakpoint() passes arguments that is not strictly needed by this function.
#' @param ordering Prior to fitting any functions, should the data be reordered by the x-axis variable or by time? Default is to use the current x-axis variable and use the time variable to break any ties.
#' @param front_trim_vt1 Seconds of data to remove from the beginning of the test before fitting (VT1). Default \code{60}.
#' @param front_trim_vt2 Seconds of data to remove from the beginning of the test before fitting (VT2). Default \code{60}.
#' @param ci Not yet implemented for this function; accepted for interface consistency with the other threshold algorithms but currently has no effect.
#' @param conf_level Not yet implemented for this function; see \code{ci}.
#' @param plots Should this function generate plots? Set to `FALSE` to save time.
#'
#' @returns A list containing breakpoint data, best-fit and related functions, and a plot.
#' @export
#'
#' @details To mimic the Santos et al. (2004) paper, this fits a 5th-order polynomial regression to the VE vs. VO2 or VCO2 data by default. Set \code{search_degree = TRUE} to instead select the degree from the data (between \code{degree} and \code{max_degree}) using the corrected AIC (\code{"aicc"}, default) or a sequential likelihood-ratio test (\code{"lrt"}). After establishing the polynomial, we find the real roots of the 2nd derivative within the given range of the `.x` variable to find the inflection points. We then filter by concavity changes that go from up to down because this matches the expected shape of the pulmonary ventilation curves.
#'
#' We are unsure how to handle multiple inflection points at this time, but this is unlikely to be an issue if the polynomial degree is kept equal to 5.
#'
#' @references
#' Santos, E. L., & Giannella-Neto, A. (2004). Comparison of computerized methods for detecting the ventilatory thresholds. European journal of applied physiology, 93, 315-324.
#'
#' @examples
#'
#' cpet_bbb <- utils::read.csv(
#'   system.file("extdata", "anton_vo2max_clean.csv", package = "gasExchangeR")
#' )
#'
#' cpet_15b_roll <- avg_exercise_test(
#'   cpet_bbb,
#'   method = "breath", calc_type = "rolling", time_col = "time", window = 15
#' )
#'
#' # VT1/GET from the inflection point of VE vs. VO2
#' d2_inflection(
#'   cpet_15b_roll,
#'   .x = "vo2",
#'   .y = "ve",
#'   bp = "vt1",
#'   vo2 = "vo2",
#'   vco2 = "vco2",
#'   ve = "ve",
#'   time = "time",
#'   plots = FALSE
#' )
#'
d2_inflection <- function(.data,
                          .x,
                          .y,
                          bp,
                          ...,
                          degree = 5,
                          search_degree = FALSE,
                          max_degree = 8,
                          criterion = c("aicc", "lrt"),
                          vo2 = "vo2",
                          vco2 = "vco2",
                          ve = "ve",
                          time = "time",
                          alpha_linearity = 0.05,
                          ordering = c("by_x", "time"),
                          front_trim_vt1 = 60,
                          front_trim_vt2 = 60,
                          ci = FALSE,
                          conf_level = 0.95,
                          plots = TRUE) {
  # check for crucial missing data
  rlang::check_required(.data)
  rlang::check_required(.x)
  rlang::check_required(.y)
  rlang::check_required(bp)

  bp <- rlang::arg_match(bp, values = c("vt1", "vt2"))
  ordering <- rlang::arg_match(ordering)
  criterion <- rlang::arg_match(criterion)

  .data <- order_cpet_df(.data, .x = .x, time = time, ordering = ordering)

  # remove the noisy start of the test before fitting (as the other threshold
  # functions do), to avoid spurious early inflection points.
  front_trim <- set_front_trim(
    bp = bp,
    front_trim_vt1 = front_trim_vt1,
    front_trim_vt2 = front_trim_vt2
  )
  .data <- .data %>%
    dplyr::filter(.data[[time]] >= min(.data[[time]] + front_trim))

  lm_poly <- fit_poly_reg(
    .data = .data, .x = .x, .y = .y,
    degree = degree, search_degree = search_degree,
    max_degree = max_degree, criterion = criterion,
    alpha_linearity = alpha_linearity
  )

  # return quick summary if generating the model fails
  if (is.null(lm_poly) || anyNA(stats::coef(lm_poly))) {
    bp_dat <- return_indeterminant_findings(
      .data = .data, bp = bp,
      algorithm = "d2_inflection",
      .x = .x, .y = .y, est_ci = "estimate"
    )
    return(list(breakpoint_data = bp_dat))
  }

  poly_expr <- expr_from_coefs(lm_poly$coefficients)
  deriv1 <- Deriv::Deriv(poly_expr, x = "x", nderiv = 1) # slope
  deriv2 <- Deriv::Deriv(poly_expr, x = "x", nderiv = 2) # acceleration

  # find x-values of roots of second derivative (inflection points)
  inflection_points_x <- deriv2 %>%
    Ryacas::y_fn("Simplify") %>% # simplify derivative for ease of extraction coefficients
    Ryacas::yac_str() %>%
    stringr::str_extract_all("(?<!\\^)-?\\d+\\.?\\d*e?-?\\d*") %>% # extract coefficients
    purrr::map(as.numeric) %>%
    unlist() %>%
    rev() %>% # reverse order for polyroot()
    polyroot() %>%
    find_real()

  # filter by roots within range of x values
  inflection_points_x <- inflection_points_x[inflection_points_x >=
    min(.data[[.x]]) &
    inflection_points_x <=
      max(.data[[.x]])]

  concavity_chgs <- concavity_changes(poly_expr,
    inflection_points = inflection_points_x
  )
  # filter by concave up to concave down
  inflection_points_x <- inflection_points_x[concavity_chgs == "up to down"]
  # if there's more than one up to down inflection point, choose whichever
  # point has the higher slope
  inflection_points_x <- inflection_points_x[which.max(
    eval(deriv1, envir = list(x = inflection_points_x))
  )]

  inflection_points_y <- eval(poly_expr, envir = list(x = inflection_points_x))

  if (length(inflection_points_x) > 0 & length(inflection_points_y) > 0) {
    bp_dat <- find_threshold_vals(
      .data = .data, thr_x = inflection_points_x,
      thr_y = inflection_points_y, .x = .x,
      .y = .y, ...
    )
  } else {
    bp_dat <- tibble::tibble()
  }

  # get values at threshold. Make this a deriv helpers function to prep bp_dat?

  if (nrow(bp_dat) == 0) { # no breakpoint found
    bp_dat <- bp_dat %>%
      dplyr::add_row() %>%
      dplyr::mutate(determinant_bp = FALSE)
    # add character or factor columns that are all the same value (e.g. ids)
    non_numeric_df <- .data %>%
      dplyr::select(tidyselect::where(function(x) {
        is.character(x) | is.factor(x) &
          all(x == x[1])
      })) %>%
      dplyr::slice(1)
    bp_dat <- dplyr::bind_cols(bp_dat, non_numeric_df)
  } else { # breakpoint found
    bp_dat <- bp_dat %>%
      dplyr::mutate(determinant_bp = TRUE)
  }
  bp_dat <- bp_dat %>%
    dplyr::mutate(
      bp = bp,
      algorithm = "d2_inflection",
      x_var = .x,
      y_var = .y
    ) %>%
    dplyr::relocate(bp, algorithm, x_var, y_var, determinant_bp)

  # make plot for output
  if (plots) {
    bp_plot <- ggplot2::ggplot(
      data = .data,
      ggplot2::aes(x = .data[[.x]], y = .data[[.y]])
    ) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::geom_line(ggplot2::aes(y = eval(poly_expr,
        envir = list(x = .data[[.x]])
      ))) +
      ggplot2::geom_vline(xintercept = bp_dat[[.x]]) +
      ggplot2::theme_minimal()
  } else {
    bp_plot <- NULL
  }

  return(list(
    breakpoint_data = bp_dat,
    lm_poly_reg = lm_poly,
    deriv1_expr = deriv1,
    deriv2_expr = deriv2,
    bp_plot = bp_plot
  ))
}

#' Find when a function changes from concave up to concave down or vice versa
#' @keywords internal
#' @noRd
#'
concavity_changes <- function(f, inflection_points, tol = 0.1) {
  # find second derivative
  f_dd <- Deriv::Deriv(f, x = "x", nderiv = 2)
  if (is.expression(f_dd)) { # change expression to function if needed
    f_dd <- expr_to_func(f_dd)
  }
  left_sign <- sign(f_dd(inflection_points - tol))
  right_sign <- sign(f_dd(inflection_points + tol))
  conc_changes <- dplyr::if_else(left_sign > 0 & right_sign < 0,
    "up to down", "down to up"
  )
  conc_changes
}
