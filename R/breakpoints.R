#' Find breakpoints in gas exchange data
#'
#' @param .data Gas exchange data
#' @param .x Variable to use on the x axis, usually time, VO2, or VCO2.
#' @param .y Variable to use on the x axis
#' @param method
#'
#' @return
#' @export
#'
#' @details
#' There are multiple ways to find a breakpoint. The first is by specifying a commonly used \code{method} and the appropriate \code{split} and \code{.x}, and \code{y}, variables will be selected for you. For example, entering 'v-slope' for \code{method} will pre-select 'v-slope' for \code{split}, 'vo2' for .x, and 'vco2' for .y. If you enter a specific \code{method}, you can override the default \code{split}, \code{.x}, and \code{.y}.
#'
#' Since there aren't common names for every way to find a breakpoint, you can also leave \code{method} as \code{NULL} and specify the \code{split}, \code{.x}, and \code{.y}.
#'
#' @examples
#' TODO write an an example
#'
#' @references
#' Beaver, W. L., Wasserman, K., & Whipp, B. J. (1986). A new method for detecting anaerobic threshold by gas exchange. Journal of Applied Physiology, 121(6), 2020–2027.
#'
breakpoint <- function(.data,
                       method = NULL,
                       split = NULL, # does split need to be helper function?
                       .x = NULL,
                       .y = NULL,
                       vo2_col = "vo2",
                       cco2_col = "vco2",
                       ve_col = "ve",
                       time_col = "time") { # need to add time_col?
    stopifnot(!missing(.data),
              !missing(method) & !missing(split) & !missing(.x) & !missing(.y))

    if(!is.null(method)) {
        method <- match.arg(method,
                           choices = c("excess_co2", "v-slope", "v_slope_simple"))
        # class(.data) <- append(class(.data), method) # could be useful to
        # use class to pre-select .x and .y variable?
    }

    # if a specific method was chosen, should a function pre-fill split, .x, and .y
    # based on the method so long as .x and .y are not already specified?

    split <- match.arg(split,
                  choices = c("dmax", "dmax_mod", "jm",
                              "orr", "v-slope", "v_slope_simple",
                              "splines"))
    # there's some lactate breakpoints that may be worth adding

}

# breakpoint.jm <- function(.data, .x, .y, method) {
#
# }


#' @keywords internal
split <- function(.data, # should loop be a helper function within split?
                 split,
                 .x = NULL,
                 .y = NULL,
                 .vo2 = "vo2",
                 .vco2 = "vco2",
                 .ve = "ve",
                 .time = "time") { # need to add time_col?

    stopifnot(!missing(.data), !missing(split), !missing(.x), !missing(.y))

    split = match.arg(split,
                      choices = c("dmax", "dmax_mod", "jm",
                                  "orr", "v-slope", "v_slope_simple",
                                  "splines"))
    class(.data) <- append(class(.data), split)
    UseMethod("split", .data)
}

#' @keywords internal
split.jm <- function(.data, # should loop be a helper function within split?
                     split,
                     .x = NULL,
                     .y = NULL,
                     .vo2 = "vo2",
                     .vco2 = "vco2",
                     .ve = "ve",
                     .time = "time") {

    min_ss_idx <- which.min(loop(.data = .data, split = split,
                                 .x = .x, .y = .y))
}

#' @keywords internal
loop <- function(.data, split, .x, .y) {
    stopifnot(!missing(.data), !missing(split), !missing(.x), !missing(.y))
    loop <- match.arg(split, choices = c("jm", "orr", "v-slope",
                                         "v_slope_simple", "splines"))
    class(.data) <- append(class(.data), loop)
    UseMethod("loop")
}

#' @keywords internal
loop.jm <- function(.data, .x, .y) {
    # browser()
    ss_both <- numeric(length = nrow(.data)-2)
    for(i in 2:(nrow(.data)-1)) {
        # split data into left and right halves. x0 = i+1
        .data_left <- .data[1:i,] # x < x0
        .data_right <- .data[(i+1):nrow(.data),] # x >= x0

        # the JM formula is y = b0 + b1*x0 + b3(x-x0). Therefore, we will create a new column in .data_right that is equal to x - x0.
        .data_right <- .data_right %>%
            mutate(s1 = .data_right[[.x]] - .data[[.x]][i+1])

        .data_right$s1 = .data_right[[.x]] - .data[[.x]][i+1]

        lm_left <- lm(.data_left[[.y]] ~ 1 + .data_left[[.x]], data = .data_left)

        # according to the JM method, the right regression line will have a constant equal to b0 + b1*x0
        b0_plus_b1x0 <- lm_left$coefficients[1] + lm_left$coefficients[2] * .data[[.x]][i+1]

        # for some reason, using I() function to try to force the regression to use b0_plus_b1x0 as an intercept doesn't work, while the offset argument does. The slope coefficient (b3) is the same, but the predicted values are not.
        lm_right <- lm(.data_right[[.y]] ~ 0 + s1, data = .data_right,
                       offset = rep(b0_plus_b1x0, nrow(.data_right)))

        # This should be the same as using the offset argument, but it's not
        # lm_right <- lm(I(.data[[.y]] - b0_plus_b1x0) ~ 0 + .data[[.x]], data = .data_right)

        ss_left <- sum((lm_left$residuals)^2)
        ss_right <- sum((lm_right$residuals)^2)
        ss_both[i-1] <- (ss_left + ss_right)

        #ss_both <- c(ss_both, ss_left + ss_right)
    }

    ss_both
}
