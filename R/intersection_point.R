#' Intersection point between two lines
#'
#' This function takes two \code{lm} objects in the form of \code{y = mx + b} and returns a vector of their \code{(x, y)} intersection point.
#' @param lm1 An \code{lm} ojbect.
#' @param lm2 A second \code{lm} ojbect.
#'
#' @return Returns a named vector of point \code{(x, y)}
#' @export
#'
intersection_point <- function(lm1, lm2) {
    stopifnot(class(lm1) == "lm",
              class(lm2) == "lm",
              length(lm1$coefficients) == 2,
              length(lm2$coefficients) == 2)
    # y = mx + b
    # mx1 + b1 = mx2 + b2 at the intersection point
    # x = (b2 - b1) / (m1 - m2)
    x <- (lm2$coefficients[1] - lm1$coefficients[1]) /
        (lm1$coefficients[2] - lm2$coefficients[2])
    y <- lm1$coefficients[1] + lm1$coefficients[2] * x

    point <- c(x, y)
    names(point) <- c('x', 'y')
    point
}
