library(gasExchangeR)
library(tidyverse)
library(janitor)
library(broom)
library(splines)
library(zoo)
# library(Deriv)
# library(Ryacas)
# library(SplinesUtils)

# big thanks to this stack overflow post
# https://stackoverflow.com/questions/44739192/export-fitted-regression-splines-constructed-by-bs-or-ns-as-piecewise-poly

# maybe try this one?
# https://stackoverflow.com/questions/29499686/how-to-extract-the-underlying-coefficients-from-fitting-a-linear-b-spline-regres

sign_changes <- function(y, change = "both") {
    # return the index of sign changes in a vector
    change <- match.arg(change,
                        choices = c("both", "maxima", "minima"),
                        several.ok = FALSE)
    if(change == "both") {
        return(which(diff(sign(diff(y))) != 0) + 1)
    } else if (change == "maxima") {
        return(which(diff(sign(diff(y))) < 0) + 1)
    } else {
        return(which(diff(sign(diff(y))) > 0) + 1)
    }
}

loop_reg_spline <- function(.data, .x, .y, df = NULL, knots = NULL,
                            degree = 3, alpha_linearity = 0.05) {
    # TODO allow users to specify b-spline or natural-spline basis
    # would that use do.call()?

    # if statement for if users specify the knots

    spline_mod_list = vector(mode = "list", length = 0)
    cont <- TRUE
    i <- 0
    # reference model without any interior knots
    lm_spline <- paste0(.y, " ~ ", "1 + ",
                        "bs(", .x, ", df = ", i + degree, ")") %>%
        as.formula() %>%
        lm(data = .data)
    spline_mod_list <- append(spline_mod_list, list(lm_spline))
    # while loop beings with 1 knot (3 df already used assuming a 3rd spline)
    while(cont == TRUE) {
        # browser()
        i <- i + 1
        # TODO add options for user-defined knots and df
        lm_spline <- paste0(.y, " ~ ", "1 + ",
                            "bs(", .x, ", df = ", i + degree, ")") %>%
            as.formula() %>%
            lm(data = .data)
        spline_mod_list <- append(spline_mod_list, list(lm_spline))
        lrt <- anova(spline_mod_list[[i]], spline_mod_list[[i+1]])
        if (is.na(lrt$`Pr(>F)`[2]) | lrt$`Pr(>F)`[2] >= alpha_linearity) {
            cont = FALSE
            lm_spline <- spline_mod_list[[i]] # take the previous model
        }
    }
    lm_spline
}

df_raw <- read_csv("inst/extdata/mar22_140_pre.csv",
                   show_col_types = FALSE)

df_unavg <- df_raw %>%
    clean_names() %>%
    select(-c(time, clock_time, ramp_time, speed, grade, ramp_speed)) %>%
    rename(speed = fixed_speeds,
           grade = fixed_grades,
           vo2_rel = vo2,
           vt = vt_btps,
           ve = ve_btps,
           time = ex_time) %>%
    mutate(ve_vo2 = ve*1000 / vo2_abs,
           ve_vco2 = ve*1000 / vco2,
           excess_co2 = (vco2^2 / vo2_abs) - vco2) %>%
    filter(!is.na(time)) %>%
    filter(stage >=6) %>%
    relocate(time, speed, grade)
df_unavg

# plot raw VO2 vs. time data
ggplot(data = df_unavg, aes(x = time, y = vo2_rel)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_line(alpha = 0.5) +
    theme_minimal() +
    ggtitle("Raw Data")

df_avg <- avg_exercise_test(df_unavg, type = "breath", subtype = "rolling",
                            time_col = "time", roll_window = 9, roll_trim = 4)

ggplot(data = df_avg, aes(x = time, y = vo2_rel)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_line(alpha = 0.5) +
    theme_minimal() +
    ggtitle("Averaged Data")

.data <- df_avg
.x <- "time"
.y <- "ve_vco2"
time <- "time"
degree <- 3
alpha_linearity = 0.05

spline_reg <- function(.data,
                       .x,
                       .y,
                       degree = NULL,
                       vo2 = "vo2",
                       vco2 = "vco2",
                       ve = "ve",
                       time = "time",
                       alpha_linearity = 0.05, # change to just alpha?
                       bp) {

    # check if there is crucial missing data
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))

    .data <- .data %>% # rearrange by x variable. Use time var to break ties.
        dplyr::arrange(.data[[.x]], .data[[time]])

    # find best number of knots
    lm_spline <- loop_reg_spline(.data = .data, .x = .x, .y = .y,
                                 degree = degree, knots = NULL)
    # find maxima in 2nd derivative. This is accomplished by making
    # a spline interpolation function. Use this with optimize() function.
    spline_func <- splinefun(x = .data[[.x]],
                            y = lm_spline$fitted.values)
    # find number of maxima
    sign_change_idx <- sign_changes(y = spline_func(x = .data[[.x]], deriv = 2),
                                   change = "maxima")
    lm_spline$fitted.values[sign_change_idx]
    # build on this for if there is more than one sign change (for loop?)
    if(length(sign_change_idx) < 2) {
        interval <- c(.data[[.x]][sign_change_idx-1],
                      .data[[.x]][sign_change_idx+1])

        extrema_deriv2 <- optimize(f = spline_func,
                                interval = interval,
                                deriv = 2, maximum = TRUE)
    }

    ggplot(data = .data, aes(x = time, y = ve_vco2)) +
        geom_line(color = "green") +
        geom_line(aes(y = lm_spline$fitted.values)) +
        # geom_vline(xintercept = extrema_deriv2$maximum,
        #            linetype = "dotted") +
        geom_line(aes(y = spline_func(x = .data[[.x]], deriv = 1)),
                  linetype = "dashed") +
        geom_line(aes(y = spline_func(x = .data[[.x]], deriv = 2)),
                  linetype = "dotted") +
        ylim(c(-2, 38)) +
        # geom_vline(xintercept = 23.5) +
        geom_vline(xintercept = .data[[.x]][sign_change_idx]) +
        theme_minimal()

    threshold_idx <- which.min(abs(.data[[.x]] - extrema_deriv2$maximum))
    .data[threshold_idx,]

    # for use with grad from numDeriv package?
    fx <- function(mod, x, col_name) {
        tib <- tibble("{col_name}" := x)
        predict(mod, newdata = tib)
    }

    fx(lm_spline, seq(20, 26, by = 0.01), .x)

    high_res_x <- seq(min(.data[[.x]]),
                      max(.data[[.x]]),
                      length.out = 1000)
    high_res_y <- predict(lm_spline, newdata = tibble("{.x}" := high_res_x))

    spline_interp_func <- splinefun(x = high_res_x, y = high_res_y)

    plot_data <- tibble("{.x}" := high_res_x) %>%
        mutate(y_hat = high_res_y,
               deriv1 = spline_interp_func(x = high_res_x, deriv = 1),
               deriv2 = spline_interp_func(x = high_res_x, deriv = 2),
        )
    plot_data

    # plot to see how the model fits the data
    ggplot(data = .data, aes(x = time)) +
        geom_line(aes(y = ve_vco2), color = "green", alpha = 0.5) +
        ylab("mm Hg") +
        xlab("Time (min)") +
        ylim(c(-1,2.5)) +
        geom_line(data = plot_data, aes(x = time, y = y_hat)) +
        geom_line(data = plot_data, aes(x = time, y = deriv1),
                  linetype = "dashed") +
        geom_line(data = plot_data, aes(x = time, y = deriv2),
                  linetype = "dotted") +
        # geom_vline(xintercept = solve(reparam_spline, b = 0, deriv = 1)) +
        # geom_vline(xintercept = solve(reparam_spline, b = 0, deriv = 2),
        #            linetype = "dashed") +
        theme_minimal()

    spline_interp_func(high_res_x, deriv=0)
    spline_interp_func(high_res_x, deriv=1)
    spline_interp_func(high_res_x, deriv=2)

    # use this? http://blog.quantitations.com/tutorial/2013/02/12/numerical-derivatives-in-r


    # given that most papers use 3rd order polynomial regression or
    # smoothing splines, I don't think I can technically take more derivatives
    # than three in order to find local maxima.
    # You probably could with a 5th order spline, but people don't do that.
    # Given that, I may need to find numerical derivatives anyway

    # x^3
    # 3x^2
    # 6x accel
    # 6 getting to the 4th derivative doesn't let you solve for the maxima
    # of the 3rd derivative.



    # I think I need to use a splineFun so I can reavaluate the points



    # is the second derivative just a bunch of linear equations
    # that meet at specific points?


    tidy(lm_spline)
    tidy(lm_spline)$estimate[-1] * .data[[.x]][1]

    reparam_spline <- RegSplineAsPiecePoly(
        lm_spline,
        SplineTerm = "bs(time, df = 6)", # needs to be dynamic
        shift = FALSE)

    reparam_spline
    reparam_spline_coefs <- reparam_spline$PiecePoly$coef
    reparam_spline_coefs[1, ] <-
        reparam_spline_coefs[1, ] + lm_spline$coefficients[1]
    poly1 <- expr_from_coefs(reparam_spline_coefs[,1])
    poly2 <- expr_from_coefs(reparam_spline_coefs[,2])
    poly3 <- expr_from_coefs(reparam_spline_coefs[,3])
    poly4 <- expr_from_coefs(reparam_spline_coefs[,4])

    library(Deriv)

    spline_dat <- .data %>%
        select(time) %>%
        mutate(y_hat = lm_spline$fitted.values,
               grad = grad(fx, x = .data[[.x]], mod = lm_spline, col_name = .x),
               d1 = spline_interp_func(time, deriv=1),
               d2 = spline_interp_func(time, deriv=2),
               d2p1 = eval(Deriv(poly1, "x", nderiv=2), list(x = time)),
               d2p2 = eval(Deriv(poly2, "x", nderiv=2), list(x = time)),
               d2p3 = eval(Deriv(poly3, "x", nderiv=2), list(x = time)),
               d2p4 = eval(Deriv(poly4, "x", nderiv=2), list(x = time)),
               )
    spline_dat

    ggplot(data = .data, aes(x = time)) +
        # geom_line(aes(y = ve_vo2), color = "purple", alpha = 0.5) +
        geom_line(aes(y = ve_vco2), color = "green", alpha = 0.5) +
        geom_line(data = spline_dat, aes(y = y_hat)) +
        ylab("mm Hg") +
        xlab("Time (min)") +
        ylim(c(-1, 3)) +
        ylim(c(-5, 40)) +
        geom_line(data = spline_dat, aes(y = grad),
                  linetype = "dashed") +
        geom_line(data = spline_dat, aes(y = d2),
                  linetype = "dotted") +
        geom_line(data = spline_dat, aes(y = d2p1),
                  color = "red", alpha = 0.5) +
        geom_line(data = spline_dat, aes(y = d2p2),
                  color ="blue", alpha = 0.5) +
        geom_line(data = spline_dat, aes(y = d2p3),
                  color = "purple", alpha = 0.5) +
        geom_line(data = spline_dat, aes(y = d2p4),
                  color = "orange", alpha = 0.5) +
        geom_vline(xintercept = solve(reparam_spline, b = 0, deriv = 1)) +
        geom_vline(xintercept = solve(reparam_spline, b = 0, deriv = 2),
                   linetype = "dashed") +
        geom_vline(xintercept = peak_d2$maximum, color = "cyan") +
        theme_minimal()


    # note to myself that I can take the crossing point of
    # some of the the 2nd derivatives
    # of the polynomial transformation of the b-spline to find the peak
    # in the 2nd derivative.

    # this seems like a good page on derivatives in R

    x <- c(1, 2, 3, 4, 5, 6, 7, 8)
    y <- c(1, 4, 9, 16, 25, 16, 9, 3)
    spline_fun <- splinefun(x, y)
    plot(x, y)

    # Use optimize() to find the maximum of the spline function over the interval [1,5]
    result <- optimize(spline_fun, c(1, 6), maximum = TRUE)
    result

    sp_func <- splinefun(x = .data[[.x]], y = lm_spline$fitted.values)

    peak_d2 <- optimize(f = spline_interp_func,
             interval = c(min(.data[[.x]]), max(.data[[.x]])),
             deriv = 2, maximum = TRUE)


    # this fucking works. Not sure how I didn't get this before
    # however, I'm pretty sure you need to switch to the next
    # polynomial when you go past a knot
    eval(poly1, envir = list(x = 18))

    knot_locations <- attr(bs(.data[[.x]],
            df = i + degree), "knots")
    knot_locations

    which(.data[[.x]] <= knot_locations[1])
    eval(poly1, envir = list(x = .data[[.x]][1:68])) -
        lm_spline$fitted.values[1:68] # essentially the same


    idx <- 70

    .data[[.x]][idx]
    eval(poly1, envir = list(x = .data[[.x]][idx]))
    eval(poly2, envir = list(x = .data[[.x]][idx]))
    lm_spline$fitted.values[idx]

    which(.data[[.x]] > knot_locations[1] & .data[[.x]] <= knot_locations[2])

    eval(poly2, envir = list(x = .data[[.x]][68:70]))

    idx = length(.data[[.x]])
    eval(poly4, envir = list(x = .data[[.x]][idx]))
    lm_spline$fitted.values[idx]


    f <- rms::ols(.data[[.y]] ~ 1 + bs(.data[[.x]],
                            df = i + degree), data = .data)
    Function(f)
    latex(f)


    # i'm stuck on getting the local maxima of the second derivative
    # how am I supposed to do that when the original order is only 3?!
    # do I need to do so numerically?


    # prepare bp output data
}

spline_func <- splinefun(x = .data[[.x]], y = lm_spline$fitted.values)
spline_func(x = .data[[.x]], deriv = 1) %>% head

num_deriv(x = .data[[.x]], y = lm_spline$fitted.values, n = 1) %>% head()

knot_locations <- attr(bs(.data[[.x]],
                          df = i + degree), "knots")
knot_locations

spline_preds <- .data %>%
    select(time) %>%
    mutate(poly1 = eval(poly1, envir = list(x = time)),
           poly2 = eval(poly2, envir = list(x = time)),
           poly3 = eval(poly3, envir = list(x = time)),
           poly4 = eval(poly4, envir = list(x = time)))
spline_preds

plot_data <- .data %>%
    select(time) %>%
    mutate(y_hat = lm_spline$fitted.values,
           deriv1 = num_deriv(x = time, y = y_hat, n = 1),
           deriv2 = num_deriv(x = time, y = y_hat, n = 2),
           deriv3 = num_deriv(x = time, y = y_hat, n = 3),
           grad = grad(fx, x = .data[[.x]], mod = lm_spline, col_name = .x),
           # deriv1_sp = spline_func(time, 1),
           # deriv2_sp = spline_func(time, 2),
           # deriv3_sp = spline_func(time, 3),
           )
plot_data

# plot to see how the model fits the data
ggplot(data = .data, aes(x = time)) +
    # geom_line(aes(y = ve_vo2), color = "purple", alpha = 0.5) +
    geom_line(aes(y = ve_vco2), color = "green", alpha = 0.5) +
    geom_vline(xintercept = knot_locations) +
    ylab("mm Hg") +
    xlab("Time (min)") +
    # ylim(c(28,38)) +
    geom_line(data = plot_data, aes(x = time, y = y_hat)) +
    # geom_line(data = spline_preds, aes(x = time, y = poly1),
    #           linetype = "dashed", color = "red") +
    # geom_line(data = spline_preds, aes(x = time, y = poly2),
    #           linetype = "dashed", color = "blue") +
    # geom_line(data = spline_preds, aes(x = time, y = poly3),
    #           linetype = "dashed", color = "purple") +
    # geom_line(data = spline_preds, aes(x = time, y = poly4),
    #           linetype = "dashed", color = "orange") +
    geom_line(data = plot_data, aes(x = time, y = deriv1),
              linetype = "dashed") +
    geom_line(data = plot_data, aes(x = time, y = deriv2),
              linetype = "dotted") +
    # geom_line(data = plot_data, aes(x = time, y = deriv3),
    #           linetype = "dotdash") +
    geom_line(data = plot_data, aes(x = time, y = grad),
              linetype = "longdash", color = "magenta") +
    # geom_line(data = plot_data, aes(x = time, y = deriv1_sp),
    #           linetype = "dashed", color = "red") +
    # geom_line(data = plot_data, aes(x = time, y = deriv2_sp),
    #           linetype = "dotted", color = "red") +
    # geom_line(data = plot_data, aes(x = time, y = deriv3_sp),
    #           linetype = "dotdash", color = "red") +
    geom_vline(xintercept = solve(reparam_spline, b = 0, deriv = 1)) +
    geom_vline(xintercept = solve(reparam_spline, b = 0, deriv = 2),
               linetype = "dashed") +
    theme_minimal()



expr_from_coefs <- function(poly_coefs, expr = TRUE) {
    string_expr <- paste("x", seq_along(poly_coefs) - 1, sep = "^")
    string_expr <- paste(string_expr, poly_coefs, sep = " * ")
    string_expr <- paste(string_expr, collapse = " + ")
    if (expr) {
        return(parse(text = string_expr))
    } else {
        return(string_expr)
    }
}

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
