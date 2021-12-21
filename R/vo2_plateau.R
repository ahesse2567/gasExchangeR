# VO2 plateau methods
# slope last x seconds is <= 50 or 150 mL O2/min (also 2.1 mL O2/kg/min)
# slope not significantly different from 0
# compare neighboring .data point to VO2 max .data point


# how to involve the "despite a change in workload" part of definition

#' Determine if a graded exercise test demonstrated a VO2 plateau
#'
#' @param .data Gas exchange .data
#' @param method Method to determine the presence of a plateau. Choose from slope at end of test (\code{slope_eot}), slope not significantly different from zero (\code{zero_sleop}), and \code{vo2max_neighbor}.
#' @param last_x_s How many seconds at the end of test will be included for slope_eot and zero_slope methods
#' @param delta_vo2 What is the acceptable variation in VO2 at the end of the test to conside if someone demonstrated a VO2 plateau.
#'
#' @return
#' @export
#'
#' @details
#'
#' See Yoon et al. (2007) for a description of the slope_eot method. Briefly, if VO2-time slope in the last 30 seconds of exercise is less than 50 mL O2/min, the participant demonstrated a plateau. This method is likely more suitable when using rolling averages or digital filtering.
#'
#' See Myers et al. (1989) and Myers et al. (1990) for a description of the zero_slope method. In short, this method tests if the slope is not significantly different from zero. If it was not significantly different from zero, the participant demonstrated a plateau. In this code, we also consider a VO2 plateau to be present if the slope was significantly different from zero while the sign is negative. Some subjects will demonstrate a plateau or a peak followed by a decrease in VO2 within the last 1-2 minutes of the test.
#'
#' See Astorino (2009), Astorino et al. (2000), and Robergs (2001) for a description of the vo2max_neighbor method. This method finds VO2max and then looks at the closest neighboring .data point. If the difference between those two points is less than the specified delta_VO2, the participant demonstrated a plateau. This method is likely better suited when using bin averaging or in historical work when there were far fewer .data points.
#'
#' @references
#' Astorino, T. A. (2009). Alterations in VO2 max and the VO2 plateau with manipulation of sampling interval. Clinical Physiology and Functional Imaging, 29(1), 60–67. https://doi.org/10.1111/j.1475-097X.2008.00835.x
#'
#' Astorino, T. A., Robergs, R. A., Ghiasvand, F., Marks, D., & Burns, S. (2000). Incidence Of The Oxygen Plateau at VO2max During Exercise Testing To Volitional Fatigue. An International Electronic Journal, 3(4), 1–12. http://eprints.qut.edu.au/96933/1/96933.pdf
#'
#' Myers, J., Walsh, D., Buchanan, N., & Froelicher, V. F. (1989). Can maximal cardiopulmonary capacity be recognized by a plateau in oxygen uptake? Chest, 96(6), 1312–1316. https://doi.org/10.1378/chest.96.6.1312
#'
#' Myers, J., Walsh, D., Sullivan, M., & Froelicher, V. (1990). Effect of sampling on variability and plateau in oxygen uptake. Journal of Applied Physiology, 68(1), 404–410. https://doi.org/10.1152/jappl.1990.68.1.404
#'
#' Robergs, R. A. (2001). An exercise physiologist’s “contemporary” interpretations of the “ugly and creaking edifices” of the VO2max concept. Journal of Exercise Physiology Online, 4(1), 1–44.
#'
#' Yoon, B. K., Kravitz, L., & Robergs, R. (2007). V̇O2max, protocol duration, and the V̇O2 plateau. Medicine and Science in Sports and Exercise, 39(7), 1186–1192. https://doi.org/10.1249/mss.0b13e318054e304
#'
#' @examples
#' # write example later
vo2_plateau <- function(.data,
                        method = c("slope_eot",
                                   "zero_slope",
                                   "vo2max_neighbor"),
                        vo2_col,
                        time_col = "time",
                        last_x_s = 30,
                        delta_vo2 = 50) {
    stopifnot(!missing(.data),
              !missing(vo2_col),
              !missing(method),
              last_x_s > 0,
              delta_vo2 > 0)

    method = match.arg(
        method,
        choices = c("slope_eot", "zero_slope", "vo2max_neighbor"))

    class(.data) <- append(class(.data), method)
    UseMethod("vo2_plateau", .data)
}

#' @export
vo2_plateau.slope_eot <- function(.data,
                        method = c("slope_eot",
                                   "zero_slope",
                                   "vo2max_neighbor"),
                        vo2_col,
                        time_col = "time",
                        last_x_s = 30,
                        delta_vo2 = 50) {
    # browser()
    eot_idx <- which((.data[[time_col]] - max(.data[[time_col]]) + last_x_s) >= 0)
    eot <- .data[eot_idx,]

    form <- as.formula(paste(vo2_col, "~ 1 +", time_col))
    lm <- lm(form, data = eot)
    if(coef(lm)[2]*60 < delta_vo2) {
        plateau <- TRUE
    } else {
        plateau <- FALSE
    }
    plateau
}

