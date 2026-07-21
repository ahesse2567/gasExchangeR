#' Determine if a graded exercise test demonstrated a VO2 plateau. Should there be an "all" test?
#'
#' @param .data Gas exchange .data
#' @param test Test(s) to determine the presence of a plateau. Choose from slope at end of test (\code{slope_eot}), slope not significantly different from zero (\code{zero_slope}), and \code{vo2max_neighbor}. Supply a character vector (e.g. \code{c("slope_eot", "zero_slope")}) to run multiple tests at once and compare whether they agree; each test contributes one row to the output.
#' @param method Choose between \code{time} (default) or \code{breath}. Only used by the \code{slope_eot} and \code{zero_slope} tests, to determine how this function finds the end of the test: does \code{last_x} refer to the last \code{last_x} seconds, or the last \code{last_x} rows (breaths)? Most published methods (e.g. Yoon et al., 2007) use time; Myers et al. (1989, 1990) used breaths. \code{method = "breath"} assumes each row of \code{.data} represents one breath and issues a warning to that effect.
#' @param vo2_col The name of the VO2 column to use in .data. Can be relative or absolute. Make sure units align with \code{units}. The default is 150 mL/min, but other papers use 50, 100, 200, and even 280 mL/min (Midgley et al., 2007).
#' @param time_col The name of the time column in .data. This function assumes your time is in seconds.
#' @param units Specify if the O2 part of the VO2 units are in mL or in L
#' @param last_x How many seconds (\code{method = "time"}) or breaths (\code{method = "breath"}) at the end of test will be included for the slope_eot and zero_slope tests.
#' @param delta_vo2 What is the acceptable variation in VO2 at the end of the test to consider if someone demonstrated a VO2 plateau.
#' @param alpha = The level of significance to use when applying the \code{zero_slope} method.
#' TODO add a plot function? This could highlight the VO2max point and locate the closest
#' # neighboring data point. For slope_eot it could add a regression line to the last_x
#' # of the graph.
#'
#' @returns A tibble with a `test` column, a `plateau` column denoting if there was a plateau by the chosen test, and additional columns specific to each test. One row is returned per test requested; when tests with different output columns are combined, unused cells are filled with `NA`.
#' @export
#'
#' @details
#'
#' See Yoon et al. (2007) for a description of the \code{slope_eot} test. Briefly, if VO2-time slope in the last 30-60 seconds of exercise is less than a specified value, the participant demonstrated a plateau. Yoon et al. (2007) used a conservative value of 50 mL/min. Most papers appear to use 150 mL/min per Taylor et al. (1955). This test is likely more suitable when using rolling averages or digital filtering.
#'
#' See Myers et al. (1989) and Myers et al. (1990) for a description of the \code{zero_slope} test. In short, this test tests if the slope is not significantly different from zero. If it was not significantly different from zero, the participant demonstrated a plateau. In this code, we also consider a VO2 plateau to be present if the slope was significantly different from zero while the sign is negative. In our personal experience, some subjects will demonstrate a plateau or a peak followed by a decrease in VO2 within the last 1-2 minutes of the test.
#' Technically this implementation is different from the Myers et al. studies. They used an 8-breath rolling average of the entire test and then performed a rolling 30-point regression of those 8-breath rolling average of VO2 vs. time. Once they had the slopes, they regressed that vs. time and determined which points were statistically significantly different from 0. To reiterate this point, they performed TWO regressions. First: VO2 vs. time using the 30 eight-breath, right-aligned rolling averages. This gave them the slope of VO2. They are repeatedly calculating the slope of the 8-breath right-aligned rolling averages across 30 of such data points. Then, they regressed this slope of VO2 against time (again). One question I have is: how much test time does a single one of these data points include? Recommendations center on calculating VO2max with no more than 30 seconds of data. If these 30-point slopes near the end of the test "stretch back" too far into the test, they might include more than 30 seconds of data. As such, that may may be too much time and you would be more likely to observe a non-zero slope of VO2.
#'
#' See Astorino (2009), Astorino et al. (2000), and Robergs (2001) for a description of the vo2max_neighbor test. This test finds VO2max and then looks at the closest neighboring data point. If the difference between those two points is less than the specified \code{delta_VO2}, the participant demonstrated a plateau. This test is likely better suited when using bin averaging or in historical work using Douglas Bags when there were far fewer data points.
#'
#' @references
#' Astorino, T. A. (2009). Alterations in VO2 max and the VO2 plateau with manipulation of sampling interval. Clinical Physiology and Functional Imaging, 29(1), 60–67. https://doi.org/10.1111/j.1475-097X.2008.00835.x
#'
#' Astorino, T. A., Robergs, R. A., Ghiasvand, F., Marks, D., & Burns, S. (2000). Incidence Of The Oxygen Plateau at VO2max During Exercise Testing To Volitional Fatigue. An International Electronic Journal, 3(4), 1–12. http://eprints.qut.edu.au/96933/1/96933.pdf
#'
#' Midgley, A. W., McNaughton, L. R., Polman, R., & Marchant, D. (2007). Criteria for Determination of Maximal Oxygen Uptake: A Brief Critique and Recommendations for Future Research. Sports Medicine, 37(12), 1019–1028. https://doi.org/10.2165/00007256-200737120-00002
#'
#' Midgley, A. W., Carroll, S., Marchant, D., McNaughton, L. R., & Siegler, J. (2009). Evaluation of true maximal oxygen uptake based on a novel set of standardized criteria. Applied Physiology, Nutrition, and Metabolism, 34(2), 115–123. https://doi.org/10.1139/H08-146
#'
#' Myers, J., Walsh, D., Buchanan, N., & Froelicher, V. F. (1989). Can maximal cardiopulmonary capacity be recognized by a plateau in oxygen uptake? Chest, 96(6), 1312–1316. https://doi.org/10.1378/chest.96.6.1312
#'
#' Myers, J., Walsh, D., Sullivan, M., & Froelicher, V. (1990). Effect of sampling on variability and plateau in oxygen uptake. Journal of Applied Physiology, 68(1), 404–410. https://doi.org/10.1152/jappl.1990.68.1.404
#'
#' Robergs, R. A. (2001). An exercise physiologist’s “contemporary” interpretations of the “ugly and creaking edifices” of the VO2max concept. Journal of Exercise Physiology Online, 4(1), 1–44.
#'
#' Taylor, H. L., Buskirk, E., & Henschel, A. (1955). Maximal Oxygen Intake as an Objective Measure of Cardio-Respiratory Performance. Journal of Applied Physiology, 8(1), 73–80. https://doi.org/10.1152/jappl.1955.8.1.73
#'
#' Yoon, B. K., Kravitz, L., & Robergs, R. (2007). VO2max, protocol duration, and the VO2 plateau. Medicine and Science in Sports and Exercise, 39(7), 1186–1192. https://doi.org/10.1249/mss.0b13e318054e304
#'
#' @examples
#'
#' # Load breath-by-breath graded exercise testing file
#' cpet_bbb <- utils::read.csv(
#'   system.file("extdata", "anton_vo2max_clean.csv", package = "gasExchangeR")
#' )
#'
#' # VO2max neighbor test
#'
#' # A 10-second bin average to showcase the VO2max neighbor test
#' cpet_10s_bin <- avg_exercise_test(
#'   cpet_bbb,
#'   test = "time",
#'   calc_type = "bin",
#'   time_col = "time",
#'   window = 10
#' )
#'
#' vo2_plat_10s_bin <- vo2_plateau(cpet_10s_bin,
#'   test = "vo2max_neighbor",
#'   vo2_col = "vo2",
#'   time_col = "time",
#'   units = "mL",
#'   last_x = 60,
#'   delta_vo2 = 150,
#'   alpha = 0.05
#' )
#'
#' vo2_plat_10s_bin
#'
#' # A 30-second bin average to showcase the VO2max neighbor test
#' cpet_30s_bin <- avg_exercise_test(
#'   cpet_bbb,
#'   test = "time",
#'   calc_type = "bin",
#'   time_col = "time",
#'   window = 30
#' )
#'
#' vo2_plat_30s_bin <- vo2_plateau(cpet_30s_bin,
#'   test = "vo2max_neighbor",
#'   vo2_col = "vo2",
#'   time_col = "time",
#'   units = "mL",
#'   last_x = 60,
#'   delta_vo2 = 150,
#'   alpha = 0.05
#' )
#'
#' vo2_plat_30s_bin
#'
#' # Notice how the presence of a VO2 plateau differs by a 10 (no plateau) vs.
#' # 30-second (yes plateau) bin average. This is due in part to the change in
#' # VO2 of 50 ml/min. Though opinions differ on the best bin average duration
#' # and VO2 differences for for determining VO2max and the presence of a VO2
#' # plateau, just be sure to clearly document your tests!
#'
#' # Slope at end-of-test test
#'
#' # We generally recommend a rolling breath average or digital filter for this
#' # approach as bin averages may result in too few data points.
#'
#' cpet_15b_roll <- avg_exercise_test(
#'   cpet_bbb,
#'   test = "breath",
#'   calc_type = "rolling",
#'   time_col = "time",
#'   window = 15
#' )
#'
#' vo2_plat_15b_roll_slope_eot <- vo2_plateau(cpet_15b_roll,
#'   test = "slope_eot",
#'   vo2_col = "vo2",
#'   time_col = "time",
#'   units = "mL",
#'   last_x = 60,
#'   delta_vo2 = 150,
#'   alpha = 0.05
#' )
#'
#' vo2_plat_15b_roll_slope_eot
#'
#' vo2_plat_15b_roll_zero_slope <- vo2_plateau(cpet_15b_roll,
#'  test = "zero_slope",
#'  vo2_col = "vo2",
#'  time_col = "time",
#'  units = "mL",
#'  last_x = 60,
#'  alpha = 0.05
#')
#'
#'vo2_plat_15b_roll_zero_slope
#'
vo2_plateau <- function(.data,
                        test,
                        method = "time",
                        vo2_col = "vo2",
                        time_col = "time",
                        units = "mL",
                        last_x = 30,
                        delta_vo2 = 150,
                        alpha = 0.05) {
  stopifnot(
    !missing(.data),
    # !missing(vo2_col),
    last_x > 0,
    delta_vo2 > 0
  )

  test <- rlang::arg_match(test,
    values = c(
      "slope_eot",
      "vo2max_neighbor",
      "zero_slope"
    ),
    multiple = TRUE
  )
  method <- rlang::arg_match(method, values = c("time", "breath"))
  units <- rlang::arg_match(units, values = c("mL", "L"))

  # Run each requested test and stack the results. When more than one test
  # is requested (e.g. test = c("slope_eot", "zero_slope")), each contributes
  # one row so the user can compare whether the tests agree.
  results <- purrr::map(test, function(m) {
    switch(m,
      slope_eot = vo2_plateau_slope_eot(
        .data, m, method, vo2_col, time_col, units, last_x, delta_vo2, alpha
      ),
      zero_slope = vo2_plateau_zero_slope(
        .data, m, method, vo2_col, time_col, units, last_x, delta_vo2, alpha
      ),
      vo2max_neighbor = vo2_plateau_vo2max_neighbor(
        .data, m, method, vo2_col, time_col, units, last_x, delta_vo2, alpha
      )
    )
  })
  dplyr::bind_rows(results)
}

#' Convert VO2 from L to mL when needed.
#' @keywords internal
#' @noRd
convert_vo2_units <- function(.data, vo2_col, units) {
  if (units == "L") {
    .data[[vo2_col]] <- .data[[vo2_col]] / 1000
  }
  .data
}

#' Subset to the last `last_x` seconds (\code{method = "time"}) or the last
#' `last_x` rows/breaths (\code{method = "breath"}) of the test.
#'
#' `method = "breath"` assumes each row of `.data` represents one breath
#' (e.g. breath-by-breath data, or a breath-based average from
#' `avg_exercise_test(test = "breath", ...)`); a one-time warning is issued
#' since this assumption can't be verified from the data itself.
#' @keywords internal
#' @noRd
end_of_test <- function(.data, time_col, last_x, method) {
  if (method == "breath") {
    rlang::warn(
      paste(
        "`method = \"breath\"` assumes each row of `.data` represents one",
        "breath (e.g. breath-by-breath data, or a breath-based average from",
        "avg_exercise_test(test = \"breath\", ...)). If your data is a",
        "time-based average, use `method = \"time\"` instead."
      ),
      .frequency = "once",
      .frequency_id = "vo2_plateau_eot_breath_method"
    )
    return(utils::tail(.data, last_x))
  }
  eot_idx <- which((.data[[time_col]] - max(.data[[time_col]]) + last_x) >= 0)
  .data[eot_idx, ]
}

#' @keywords internal
#' @noRd
vo2_plateau_slope_eot <- function(.data, test, method, vo2_col, time_col,
                                  units, last_x, delta_vo2, alpha) {
  .data <- convert_vo2_units(.data, vo2_col, units)
  eot <- end_of_test(.data, time_col, last_x, method)

  form <- stats::as.formula(paste(vo2_col, "~ 1 +", time_col))
  lm_vo2_time <- lm(form, data = eot)
  vo2_time_slope_min <- stats::coef(lm_vo2_time)[2] * 60
  plateau <- vo2_time_slope_min < delta_vo2
  tibble::tibble(
    test = test,
    plateau = plateau,
    vo2_time_slope_min = vo2_time_slope_min,
    delta_vo2 = delta_vo2,
    p.value = broom::tidy(lm_vo2_time)[["p.value"]][2]
  )
}

#' @keywords internal
#' @noRd
vo2_plateau_zero_slope <- function(.data, test, method, vo2_col, time_col,
                                   units, last_x, delta_vo2, alpha) {
  .data <- convert_vo2_units(.data, vo2_col, units)
  eot <- end_of_test(.data, time_col, last_x, method)

  form <- stats::as.formula(paste(vo2_col, "~ 1 +", time_col))
  lm <- lm(form, data = eot)

  if (broom::tidy(lm)[["p.value"]][2] > alpha) {
    plateau <- TRUE
  } else if (broom::tidy(lm)[["p.value"]][2] < alpha & stats::coef(lm)[2] * 60 < 0) {
    plateau <- TRUE
  } else {
    plateau <- FALSE
  }
  tibble::tibble(
    test = test,
    plateau = plateau,
    vo2_time_slope_min = stats::coef(lm)[2] * 60,
    p.value = broom::tidy(lm)[["p.value"]][2]
  )
}

#' @keywords internal
#' @noRd
vo2_plateau_vo2max_neighbor <- function(.data, test, method, vo2_col, time_col,
                                        units, last_x, delta_vo2, alpha) {
  .data <- convert_vo2_units(.data, vo2_col, units)

  vo2max_idx <- which.max(.data[[vo2_col]])
  vo2max <- max(.data[[vo2_col]])

  if (vo2max_idx == nrow(.data)) {
    # this is a special case when the VO2max index is the last data point
    vo2max_neighbor_diff <- vo2max -
      .data[[vo2_col]][vo2max_idx - 1]
  } else {
    # find the neighbors, i.e. data points on either side of VO2max
    t_diffs <- abs(c(
      .data[[time_col]][vo2max_idx - 1],
      .data[[time_col]][vo2max_idx + 1]
    ) - .data[[time_col]][vo2max_idx])

    # don't use which.min() for this b/c that only returns 1 value, and with
    # time averages there will often be ties
    min_t_diff_idx <- which(t_diffs == min(t_diffs))

    if (length(min_t_diff_idx) > 1) {
      # there's a tie for the nearest data point with respect to time,
      # so find the VO2 (y-value) at each index and use the one closest
      # to VO2max as the comparison
      closer_neighbor_idx <- which.max(c(
        .data[[vo2_col]][vo2max_idx - 1],
        .data[[vo2_col]][vo2max_idx + 1]
      ))
      # higher VO2 value is to the left
      if (closer_neighbor_idx == 1) {
        vo2max_neighbor_diff <- vo2max -
          .data[[vo2_col]][vo2max_idx - 1]
      } else { # higher VO2 value is to the right
        vo2max_neighbor_diff <- vo2max -
          .data[[vo2_col]][vo2max_idx + 1]
      }
    } else {
      if (min_t_diff_idx == 1) {
        # closer point is to the left
        vo2max_neighbor_diff <- vo2max -
          .data[[vo2_col]][vo2max_idx - 1]
      } else {
        # closer point is to the right
        vo2max_neighbor_diff <- vo2max -
          .data[[vo2_col]][vo2max_idx + 1]
      }
    }
  }

  if (vo2max_neighbor_diff < delta_vo2) {
    plateau <- TRUE
  } else {
    plateau <- FALSE
  }

  out <- tibble::tibble(
    test = test,
    plateau = plateau,
    delta_vo2 = delta_vo2,
    vo2max_neighbor_diff = vo2max_neighbor_diff
  )
  out
}

#' Determine if a graded exercise test demonstrated VO2max according to "secondary" criteria. Importantly, there is considerable criticism and controversy over using secondary criteria to reach VO2max (Midgley et al., 2009).
#'
#' @param .data Gas exchange .data
#' @param age Participant age in years.
#' @param sex Participant sex, \code{"male"} or \code{"female"}. Only required
#'   if \code{aphrm_eq} includes a sex-specific equation (\code{gulati},
#'   \code{fairbarn-m}, or \code{fairbarn-f}). Leave as \code{NULL} (default) if you only need the
#'   sex-agnostic equations (\code{fox}, \code{gellish}, \code{tanaka},
#'   \code{arena}, \code{astrand}, \code{nes}).
#' @param aphrm_eq Age predicted heart rate (HR) maximum (APHRM) equation(s)
#'   to use. Choose one or more of \code{fox}, \code{gellish}, \code{tanaka},
#'   \code{arena}, \code{astrand}, \code{nes}, \code{gulati} (female only),
#'   \code{fairbarn-m} (male only), or \code{fairbarn-f} (female only). Set to
#'   \code{"all"} to compute every equation available given \code{sex}: if
#'   \code{sex} is \code{NULL}, this runs only the sex-agnostic equations; if
#'   \code{sex} is supplied, sex-specific equations for that sex are included
#'   too.
#' @param hr_cutoff How close to APHRM should the participant get to consider
#'   reaching a VO2max effort?
#' @param hr_unit Will the HR cutoff be determined using a percentage
#'   (\code{pct}) of APHRM or a specific number of \code{beats}? For
#'   \code{pct}, the criterion is met when measured HRmax is at least
#'   \code{hr_cutoff} times APHRM. For \code{beats}, the criterion is met when
#'   measured HRmax is within \code{hr_cutoff} beats of APHRM in either
#'   direction.
#' @param hr_col The name of the heart rate (HR) column to use in .data.
#' @param rer_col The name of the respiratory exchange ratio (RER) column to use in .data.
#' @param rer_cutoff What is the required RER to consider reaching VO2max?
# #' @param RPE_col The name of the rating of perceived exertion (RPE) column to use in .data. Unused at this time.
# #' @param BLa_col The name of the blood lactate column to use in .data. Unused at this time.
#'
#' @returns A tibble with one row per requested APHRM equation, containing
#'   \code{aphrm_eq}, \code{equation_sex} (the sex the equation's formula
#'   applies to, or \code{NA} if sex-agnostic), \code{age}, \code{hr_max}
#'   (measured), \code{aphrm} (predicted), \code{hr_cutoff}, \code{hr_unit},
#'   \code{hr_criteria_met}, \code{rer_max} (measured), \code{rer_cutoff},
#'   \code{rer_criteria_met}, and \code{criteria_met} (\code{TRUE} only if
#'   both the HR and RER criteria were met).
#' @export
#'
#' @references
#' Midgley, A. W., Carroll, S., Marchant, D., McNaughton, L. R., & Siegler, J. (2009). Evaluation of true maximal oxygen uptake based on a novel set of standardized criteria. Applied Physiology, Nutrition, and Metabolism, 34(2), 115–123. https://doi.org/10.1139/H08-146
#'
#' Shookster, D., Lindsey, B., Cortes, N., & Martin, J. R. (2020). Accuracy of commonly used age-predicted maximal heart rate equations. International Journal of Exercise Science, 13(7), 1242–1250.
#'
#' @details
#' The nine APHRM equations come from Shookster et al. (2020), who compared
#' them against measured HRmax from graded exercise tests:
#' \itemize{
#'   \item fox: HRmax = 220 - Age
#'   \item gellish: HRmax = 207 - 0.7 x Age
#'   \item tanaka: HRmax = 208 - 0.7 x Age
#'   \item arena: HRmax = 209.3 - 0.72 x Age
#'   \item astrand: HRmax = 216.6 - 0.84 x Age
#'   \item nes: HRmax = 211 - 0.64 x Age
#'   \item gulati (female only): HRmax = 206 - 0.88 x Age
#'   \item fairbarn-m (male): HRmax = 208 - 0.8 x Age
#'   \item fairbarn-f (female): HRmax = 201 - 0.63 x Age
#' }
#'
#' @examples
#'
#' cpet_bbb <- utils::read.csv(
#'   system.file("extdata", "anton_vo2max_clean.csv", package = "gasExchangeR")
#' )
#'
#' # Selecting two common heart APMHR equations.
#' vo2_max_secondary(
#'   cpet_bbb,
#'   age = 28,
#'   aphrm_eq = c("fox", "tanaka"),
#'   hr_cutoff = 0.9,
#'   hr_unit = "pct",
#'   hr_col = "hr",
#'   rer_col = "rer",
#'   rer_cutoff = 1.1
#' )
#'
#' # All sex-agnostic equations only (sex left as NULL)
#' vo2_max_secondary(
#'   cpet_bbb,
#'   age = 28,
#'   aphrm_eq = "all",
#'   hr_cutoff = 0.9,
#'   hr_unit = "pct",
#'   hr_col = "hr",
#'   rer_col = "rer",
#'   rer_cutoff = 1.1
#' )
#'
#' # Include sex-specific equations by supplying sex
#' vo2_max_secondary(
#'   cpet_bbb,
#'   age = 28,
#'   sex = "male",
#'   aphrm_eq = "all",
#'   hr_cutoff = 10,
#'   hr_unit = "beats",
#'   hr_col = "hr",
#'   rer_col = "rer",
#'   rer_cutoff = 1.1
#' )
#'
vo2_max_secondary <- function(.data,
                              age,
                              sex = NULL,
                              aphrm_eq = "fox",
                              hr_cutoff = 0.9,
                              hr_unit = "pct",
                              hr_col = "hr",
                              rer_col = "rer",
                              rer_cutoff = 1.1) {
  rlang::check_required(.data)
  rlang::check_required(age)
  if (!is.null(sex)) {
    sex <- rlang::arg_match(sex, values = c("male", "female"))
  }
  hr_unit <- rlang::arg_match(hr_unit, values = c("pct", "beats"))

  eqs <- resolve_aphrm_eq(aphrm_eq, sex)

  hr_max <- max(.data[[hr_col]], na.rm = TRUE)
  rer_max <- max(.data[[rer_col]], na.rm = TRUE)
  rer_criteria_met <- rer_max >= rer_cutoff

  eqs %>%
    dplyr::mutate(
      age = age,
      aphrm = intercept - slope * age,
      hr_max = hr_max,
      hr_criteria_met = switch(hr_unit,
        pct = hr_max >= hr_cutoff * aphrm,
        beats = abs(hr_max - aphrm) <= hr_cutoff
      ),
      hr_cutoff = hr_cutoff,
      hr_unit = hr_unit,
      rer_max = rer_max,
      rer_cutoff = rer_cutoff,
      rer_criteria_met = rer_criteria_met,
      criteria_met = hr_criteria_met & rer_criteria_met
    ) %>%
    dplyr::select(
      aphrm_eq, equation_sex = sex_required, age, hr_max, aphrm,
      hr_cutoff, hr_unit, hr_criteria_met, rer_max, rer_cutoff,
      rer_criteria_met, criteria_met
    )
}

#' Table of published age-predicted maximal heart rate (APHRM) equations.
#'
#' See Shookster et al. (2020) for the source of each equation.
#' @keywords internal
#' @noRd
aphrm_equations <- function() {
  tibble::tribble(
    ~aphrm_eq, ~sex_required, ~intercept, ~slope,
    "fox", NA_character_, 220, 1,
    "gellish", NA_character_, 207, 0.7,
    "tanaka", NA_character_, 208, 0.7,
    "arena", NA_character_, 209.3, 0.72,
    "astrand", NA_character_, 216.6, 0.84,
    "nes", NA_character_, 211, 0.64,
    "gulati", "female", 206, 0.88,
    "fairbarn-m", "male", 208, 0.8,
    "fairbarn-f", "female", 201, 0.63
  )
}

#' Resolve which APHRM equations to compute, given the requested `aphrm_eq`
#' and (optionally) the participant's `sex`.
#'
#' If `aphrm_eq == "all"`, returns every sex-agnostic equation, plus
#' sex-specific equations for `sex` if `sex` is supplied. Otherwise, resolves
#' each requested equation name, erroring if a sex-specific equation is
#' requested without `sex`, or with a `sex` it isn't defined for.
#'
#' @keywords internal
#' @noRd
resolve_aphrm_eq <- function(aphrm_eq, sex) {
  all_eqs <- aphrm_equations()
  sex_agnostic <- all_eqs[is.na(all_eqs$sex_required), ]

  if (identical(aphrm_eq, "all")) {
    if (is.null(sex)) {
      return(sex_agnostic)
    }
    sex_specific <- all_eqs[
      !is.na(all_eqs$sex_required) & all_eqs$sex_required == sex,
    ]
    return(dplyr::bind_rows(sex_agnostic, sex_specific))
  }

  requested <- unique(aphrm_eq)
  valid_names <- unique(all_eqs$aphrm_eq)
  unknown <- setdiff(requested, valid_names)
  if (length(unknown) > 0) {
    rlang::abort(paste0(
      "Unknown `aphrm_eq` value(s): ", paste(unknown, collapse = ", "),
      ". Choose from: ", paste(valid_names, collapse = ", "), ", or \"all\"."
    ))
  }

  matched_rows <- lapply(requested, function(eq) {
    rows <- all_eqs[all_eqs$aphrm_eq == eq, ]
    if (all(is.na(rows$sex_required))) {
      return(rows)
    }
    if (is.null(sex)) {
      rlang::abort(paste0(
        "`sex` must be supplied to use the \"", eq, "\" equation."
      ))
    }
    matched <- rows[!is.na(rows$sex_required) & rows$sex_required == sex, ]
    if (nrow(matched) == 0) {
      defined_for <- paste(rows$sex_required, collapse = " or ")
      rlang::abort(paste0(
        "The \"", eq, "\" equation is only defined for ", defined_for,
        " participants; `sex` was \"", sex, "\"."
      ))
    }
    matched
  })

  dplyr::bind_rows(matched_rows)
}
