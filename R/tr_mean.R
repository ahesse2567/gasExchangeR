#' Trimmed mean helper function to emulate the 'mid 5 of 7' and similar averaging methods offered by MGC and Breeze
#'
#' @details Function for internal use only. It will not be exported.
#'
#' @param .v A vector of numeric or logical values
#' @param x The number of values to keep
#' @param y The total number of values to consider
#'
#' @return
#' @keywords internal
#'
#' @details This function trims the mean according to the difference between the x and y parameters. For example, if \code{x = 5} and \code{y = 7}, the difference is\code{2}. Half of the difference value is removed from the highest side of the vector and half the difference is removed from the lowest side of the vector. Therefore, the higest and lowest values would be removed before taking the mean. If instead \code{x = 5} and \code{y = 7}, the two lowest and two highest values would be removed before taking the mean.
#'
#' @examples
#'
#' set.seed(2357234)
#' v <- 1:7
#' v <- round(v + rnorm(n = length(v), sd = 3), 2)
#' v
#' tr_mean(v, x = 5, y = 7)
#'
tr_mean <- function(.v, x, y) {
    # idx_col/decision_col = "ve_vo2_prod" should be implemented in the future
    # because Breeze uses the ve * vco2 to determine which rows to remove
    stopifnot(is.numeric(.v) | is.logical(.v),
              x >= 1 & x %% 1 ==0 & x %% 2 == 1,
              y >= 1 & y %% 1 ==0 & y %% 2 == 1 & y - x >= 2)
    # browser()
    num2remove <- y - x
    for (i in 1:(num2remove / 2)) {
        min_val <- which(.v == min(.v, na.rm = TRUE))
        if (length(min_val) > 1) { # if two values are tied, remove one at random
            min_val <- sample(min_val, 1)
        }
        .v <- .v[-min_val]
        max_val <- which(.v == max(.v, na.rm = TRUE))
        if (length(max_val) > 1) { # if two values are tied, remove one at random
            max_val <- sample(max_val, 1)
        }
        .v <- .v[-max_val]
        out <- mean(.v)
    }
    out
}
