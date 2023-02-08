library(SplinesUtils)

reparam_spline <- RegSplineAsPiecePoly(
    lm_spline,
    SplineTerm = "bs(time, df = 6)", # needs to be dynamic
    shift = FALSE)

num_deriv <- function(x, y, n = 1L) {
    # browser()
    # this finds the numerical derivative with respect to x
    dy <- diff(zoo(y), na.pad = TRUE)
    dx <- diff(zoo(x), na.pad = TRUE)
    dydx <- dy/dx
    if(n > 1) {
        dydx <- num_deriv(x = x, y = dydx, n = n - 1)
    }
    dydx
}



x <- seq(-7, 5, by = 0.1)
y <- (x - 4)*(x + 2)*(x - 1)*(x + 6)
plot(x, y)

my_expr <- expression((x - 4)*(x + 2)*(x - 1)*(x + 6))
my_func <- function(x) (x - 4)*(x + 2)*(x - 1)*(x + 6)

change = "maxima"
threshold = 1e-6
# our function will be diff(y) / diff(x)? or will it be
# a splineFun?
find_extrema <- function(x, fx, change = "both", threshold = 1e-6, ...) {
    # find extreme within a certain resolution
    # x is a vector of x coordinates
    # fx is a function that gives y values based on x as an input
    browser()
    change <- match.arg(change,
                        choices = c("both", "maxima", "minima"),
                        several.ok = FALSE)

    if(is.function(fx)) {
        y <- fx(x)
    } else if (is.expression(fx)) {
        y <- eval(fx, envir = list(x = x))
    } else {
        stop()
    }
    # find the indices closest to maxima
    maxima_idx <- sign_changes(y, change = change)
    # TODO will need to update to for loop at some point

    # get x points on either side of index
    x_either_side <- x[c(maxima_idx-1, maxima_idx+1)]
    # increase resolution
    new_x <- increase_resolution(x_either_side, by_order = 1)
    new_y <- fx(new_x)

    new_maxima_idx <- sign_changes(new_y, change = change)
    diffs <- c(new_y[new_maxima_idx] - new_y[new_maxima_idx-1],
               new_y[new_maxima_idx] - new_y[new_maxima_idx+1])

    if(all(diffs < threshold)) {
        return(c("x" = new_x[new_maxima_idx],
                 "y" = new_y[new_maxima_idx]))
    } else { # call functions recursively until until threshold is met
        find_extrema(x = new_x, fx = fx,
                     change = change, threshold = threshold)
    }
}

x <- seq(-7, 5, by = 0.1)
my_func <- function(x) (x - 4)*(x + 2)*(x - 1)*(x + 6)
y <- my_func(x)
plot(x, y)
find_extrema(x = x, fx = my_func, change = "maxima")

current_log10_order <- function(x) {
    if(x < 0) {
        x <- abs(x)
    }
    if(x == 0) {
        return(0)
    }
    floor(log10(x))
}

increase_resolution <- function(x, by_order = 1) {
    stopifnot(length(x) == 2) # x must be a vector
    new_resolution <- 10^(current_log10_order(diff(x)) - by_order)
    out <- seq(x[1], x[2], by = new_resolution)
    out
}
