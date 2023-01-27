library(gasExchangeR)
library(devtools)
library(tidyverse)
# library(devtools)
library(janitor)
library(broom)
library(Deriv)
library(Ryacas)
library(readxl)

# Different authors who use derivatives define the threshold differently
# Santos et al. (2004) found the 2nd derivative of a 5th degree polynomial and found where that 2nd derivative function crossed 0. Those were inflection point (concavity changes) in the 5th order polynomial function. They defined those points as the threshold.

# The study by Leo et al. (2017) was somewhat different from Santos by my reading because they found a local maxima in the 2nd derivative. That is, they essentially found the zero-crossing point (root) of the third derivative to find the maxima in the second derivative. # 10.1249/MSS.0000000000001226

# Cross et al. (2012) also used extrema in the 2nd derivative, in this case a 6th degree polynomial spline fit. Unfortunately, this paper did NOT mention anything about knots, so I'm curious if they actually just used a plain polynomial regression.

# Wisén and Wohlfart (2004) found the first derivative of the vo2 and vco2 vs. time curves using a 6th degree polynomial and found the last time where the vco2 vs. time derivative surpassed the vo2 vs. time derivative.


# big thanks to this stack overflow thread: https://stackoverflow.com/questions/40438195/function-for-polynomials-of-arbitrary-order-symbolic-method-preferred

df_colnames <- read_xlsx("inst/extdata/Foreman_Ramp_10_18_2022.xlsx",
                         n_max = 0) %>%
    colnames()
df_raw <- read_xlsx("inst/extdata/Foreman_Ramp_10_18_2022.xlsx",
                    col_names = FALSE, skip = 3) %>%
    set_names(df_colnames)

df_unavg <- df_raw %>%
    clean_names() %>%
    rename(time = t) %>%
    mutate(time = lubridate::minute(time) * 60 + lubridate::second(time)) %>%
    filter(phase == "EXERCISE") %>%
    relocate(time, speed, grade)

df_avg <- avg_exercise_test(df_unavg, type = "breath", subtype = "rolling",
                            time_col = "time", roll_window = 9, roll_trim = 4)

# plot raw VO2 vs. time data
ggplot(data = df_avg, aes(x = time, y = vo2)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_line(alpha = 0.5) +
    geom_point(aes(y = vco2), alpha = 0.5, color = "blue") +
    geom_line(aes(y = vco2), alpha = 0.5) +
    theme_minimal() +
    ggtitle("Raw Data")


ggplot(data = df_avg, aes(x = vo2, y = vco2)) +
    geom_point(color = "blue", alpha = 0.5) +
    scale_x_continuous(limits = c(0, max(df_avg$vco2))) +
    scale_y_continuous(limits = c(0, max(df_avg$vco2))) +
    ggtitle("V-Slope") +
    theme_minimal()

bp_dat <- breakpoint(df_avg,
           algorithm_vt1 = "v-slope",
           x_vt1 = "vo2",
           y_vt1 = "vco2",
           algorithm_vt2 = "jm",
           x_vt2 = "vco2",
           y_vt2 = "ve",
           vo2 = "vo2")
bp_dat$bp_dat

ggplot(data = df_avg, aes(x = vco2, y = ve)) +
    geom_point(color = "orange", alpha = 0.5) +
    theme_minimal() +
    geom_vline(xintercept = bp_dat$bp_dat %>%
                   filter(bp == "vt2") %>%
                   select(vco2) %>%
                   pull())

ggplot(data = df_avg, aes(x = time)) +
    geom_point(aes(y = ve_vo2), color = "purple", alpha = 0.5) +
    geom_line(aes(y = ve_vo2), color = "purple", alpha = 0.5) +
    geom_point(aes(y = ve_vco2), color = "green", alpha = 0.5) +
    geom_line(aes(y = ve_vco2), color = "green", alpha = 0.5) +
    # geom_vline(xintercept = bp_dat$bp_dat$time) +
    theme_minimal() +
    ggtitle("Ventilatory Equivalents")

loop_poly_reg <- function(.data, .x, .y,
                          degree = NULL, alpha_linearity = 0.05) {
    # browser()
    # if the user specifies a degree, find that and be done with it
    if (!is.null(degree)) {
        lm_poly <- lm(.data[[.y]] ~ 1 + poly(.data[[.x]], degree = degree,
                                             raw = TRUE), data = .data)
        # if the user does NOT specify a degree, find the best degree using
        # likelihood ratio test
    } else {
        degree = 5 # start at degree = 5 so you can take 4 derivatives
        lm_list = vector(mode = "list", length = 0) # hold lm data
        # keep raw = TRUE for now b/c it's way easier for me to test
        # if the derivative expressions are working
        # starting with degree = 5 b/c that's the minimum I've seen in
        # other papers that seemed to use polynomial fits
        lm_poly <- lm(.data[[.y]] ~ 1 + poly(.data[[.x]], degree = degree,
                                             raw = TRUE), data = .data)
        lm_list <- append(lm_list, list(lm_poly))

        cont <- TRUE
        i <- 1 # start at 2 b/c we already made linear (degree = 1) model
        while(cont == TRUE) {
            lm_poly <- lm(.data[[.y]] ~ 1 + poly(.data[[.x]],
                                                 degree = degree + i,
                                                 raw = TRUE),
                          data = .data)
            lm_list <- append(lm_list, list(lm_poly))
            lrt <- anova(lm_list[[i]], lm_list[[i+1]])
            if (is.na(lrt$`Pr(>F)`[2]) | lrt$`Pr(>F)`[2] >= alpha_linearity) {
                cont = FALSE
                lm_poly <- lm_list[[i]] # take the previous model
            }
            i <- i + 1
        }
    }

    lm_poly
}

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

find_real_roots <- function(v, threshold = 1e-6) {
    # find real roots by fixing rounding errors
    Re(v)[abs(Im(v)) < threshold]
}

.x <- "time"
.y <- "ve_vco2"
time <- "time"
.data <- df_avg
# .data <- df_avg[1:threshold_idx,]
degree <- NULL
alpha_linearity = 0.05
bp = "vt2"

poly_regression <- function(.data,
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

    lm_poly <- loop_poly_reg(.data = .data, .x = .x, .y = .y,
                             degree = degree,
                             alpha_linearity = alpha_linearity)

    # TODO what about if the best model is linear? Then this method probably
    # won't work well and you should get a warning or an error about that
    # also, I think you need a minimum of a 4th order equation if you want to
    # find the maxima in the acceleration (2nd derivative) of the relationship
    # b/c you essentially need to take a third derivative and still have an x term

    # lm_poly <- lm_list[[3]] # keep for testing, delete later
    # this definitely breaks if the best-fit polynomial order is too low
    poly_expr <- expr_from_coefs(lm_poly$coefficients)
    deriv1 <- Deriv(poly_expr, x = "x", nderiv = 1) # slope
    deriv2 <- Deriv(poly_expr, x = "x", nderiv = 2) # acceleration
    deriv3 <- Deriv(poly_expr, x = "x", nderiv = 3) # jerk
    deriv4 <- Deriv(poly_expr, x = "x", nderiv = 4) # snap. This kinda feels hacky and dumb at this point

    # find x-values of roots
    roots_deriv3 <- deriv3 %>%
        y_fn("Simplify") %>% # simplify derivative for ease of extraction coefficients
        yac_str() %>%
        str_extract_all("(?<!\\^)-?\\d+\\.?\\d*e?-?\\d*") %>% # extract coefficients
        map(as.numeric) %>%
        unlist() %>%
        rev() %>% # reverse order for polyroot()
        polyroot() %>%
        find_real_roots()

    # filter by roots within range of x values
    roots_deriv3 <- roots_deriv3[roots_deriv3 >= min(.data[[.x]]) &
                     roots_deriv3 <= max(.data[[.x]])]

    # use 4th derivative to find which are maxima (y @ roots < 0)
    local_maxima_x <- roots_deriv3[eval(deriv4,
                                      envir = list(x = roots_deriv3)) < 0]
    local_maxima_y <- eval(deriv2, envir = list(x = local_maxima_x))

    # you would generally expect the threshold to be the highest of any possible
    # local maxima (at least for VT2)

    # select higher of the two local maxima b/c we expect just one large change
    # I feel like this goes against my visual detection
    # according to Cross et al. (2012), one should select the lower of any two
    # maxima for VT1 and the higher of any two maxima for VT2
    if(length(local_maxima_x) > 1) {
        if(bp == "vt1") {
            min_local_maxima_x <- local_maxima_x[which.min(local_maxima_y)]
            threshold_idx <-
                which.min(abs(.data[[.x]] - min_local_maxima_x))
        }
        if(bp == "vt2") {
            max_local_maxima_x <- local_maxima_x[which.max(local_maxima_y)]
            threshold_idx <- which.min(abs(.data[[.x]] - max_local_maxima_x))
        }
    } else {
        threshold_idx <- which.min(abs(.data[[.x]] - local_maxima_x))
    }

    # get values at threshold
    bp_dat <- .data[threshold_idx,] %>%
        mutate(bp = bp,
               algorithm = "poly_reg",
               x_var = .x,
               y_var = .y,
               # determinant_bp = determinant_bp,
               # pct_slope_change = pct_slope_change,
               # f_stat = f_stat,
               # p_val_f = pf_two,
               ) %>%
        relocate(bp, algorithm, x_var, y_var,
                 # determinant_bp,
                 # pct_slope_change, f_stat, p_val_f
                 )
    bp_dat
}

plot_data <- .data %>%
    select(all_of(.x)) %>%
    mutate(y_hat = lm_poly$fitted.values,
           y_hat_deriv1 = eval(deriv1,
                                      envir = list(x = .data[[.x]])),
           y_hat_deriv2 = eval(deriv2,
                                      envir = list(x = .data[[.x]])),
           # y_hat_deriv3 = eval(deriv3,
           #                     envir = list(x = .data[[.x]])),
           )
plot_data

ggplot(data = .data, aes(x = .data[[.x]])) +
    geom_point(aes(y = .data[[.y]]), alpha = 0.5) +
    geom_line(aes(y = lm_poly$fitted.values)) +
    theme_minimal()

# plot to see how the model fits the data
ggplot(data = .data, aes(x = .data[[.x]])) +
    geom_point(aes(y = .data[[.y]]), alpha = 0.5) +
    geom_line(data = plot_data, aes(y = y_hat)) +
    # ylim(c(-2.5,5)) +
    # geom_line(data = plot_data,
    #           aes(x = time, y = y_hat_deriv1),
    #           linetype = "dashed") +
    # geom_line(data = plot_data,
    #           aes(x = time, y = y_hat_deriv2), linetype = "dotted") +
    geom_vline(xintercept = local_maxima_x, linetype = "dotted") +
    # geom_vline(xintercept = bp_dat$bp_dat %>%
    #                filter(bp == "vt2") %>%
    #                select(time) %>%
    #                pull()) +
    # geom_line(aes(y = vo2_kg), color = "red") +
    # geom_line(aes(y = ve_vo2), color = "purple") +
    # ylim(c(-3, 50)) +
    theme_minimal()

# zoomed in on derivatives
ggplot(data = .data, aes(x = .data[[.x]])) +
    # geom_line(data = plot_data,
    #           aes(x = time, y = y_hat_deriv1), linetype = "dashed") +
    geom_line(data = plot_data,
              aes(x = time, y = y_hat_deriv2), linetype = "dotted") +
    # geom_line(data = plot_data, aes(y = y_hat_deriv3)) +
    geom_vline(xintercept = local_maxima_x, linetype = "dotted") +
    geom_hline(yintercept = local_maxima_y, linetype = "dotted") +
    theme_minimal()

vo2max <- max(.data$vo2_kg)

.data[threshold_idx,] %>%
    select(time, speed, grade, rq, vo2_kg) %>%
    mutate(vo2max = vo2max,
           pct_vo2_max = vo2_kg / vo2max)

bp_dat$bp_dat %>%
    select(time, speed, grade, rer, vo2_rel) %>%
    mutate(pct_vo2max = vo2_rel / vo2max)

