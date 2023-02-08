library(gasExchangeR)
library(tidyverse)
library(janitor)
library(broom)
library(Deriv)
library(Ryacas)
library(readxl)
library(devtools)

# This method is closest to Santos et al. (2004)
# does this inflection point method only apply to VE vs. VO2 or VE vs. VCO2?
# in my anecdotal experience analyzing Nick's data, VE vs. Time does not work well

# will we need some sort of curve type argument for derivative methods?
# e.g., poly, reg spline, smoothing spline?

# note to Anton: inflection points seem like (almost?) always different from
# maxima. Would that introduce a systematic bias between then?

# Does this (and maybe other derivative functions) need to assess if there aren't any
# inflection points within the interval, and if there aren't, try a different polynomial
# regression order?

# TBH this doesn't always seem to work if the polynomial order is different than 5

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
    relocate(time, speed, grade) %>%
    ventilatory_outliers()

df_avg <- avg_exercise_test(df_unavg, type = "breath", subtype = "rolling",
                            time_col = "time", roll_window = 15)

# plot raw VO2 vs. time data
ggplot(data = df_avg, aes(x = time, y = vo2)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_line(alpha = 0.5) +
    geom_point(aes(y = vco2), alpha = 0.5, color = "blue") +
    geom_line(aes(y = vco2), alpha = 0.5) +
    theme_minimal()

ggplot(data = df_avg, aes(x = time, y = vo2_kg)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_line(alpha = 0.5) +
    geom_point(aes(y = ve), alpha = 0.5, color = "orange") +
    theme_minimal() +
    ggtitle("Raw Data")

bp_ve_vs_vco2 <- d2_inflection(.data = df_avg, .x = "vco2", .y = "ve", bp = "vt2",
                               degree = 5)
bp_ve_vs_time <- d2_inflection(.data = df_avg, .x = "time", .y = "ve", bp = "vt2")

ggplot(data = df_avg, aes(x = vco2, y = ve)) +
    geom_point(color = "orange", alpha = 0.5) +
    geom_vline(xintercept = bp_ve_vs_vco2$vco2)

ggplot(data = df_avg, aes(x = time, y = ve_vco2)) +
    geom_point(color = "green", alpha = 0.5) +
    geom_vline(xintercept = bp_ve_vs_vco2$time)



# assuming we're using polynomial functions to make the model for the time being
d2_inflection <- function(.data,
                          .x,
                          .y,
                          bp,
                          degree = NULL,
                          vo2 = "vo2",
                          vco2 = "vco2",
                          ve = "ve",
                          time = "time",
                          alpha_linearity = 0.05) # change to just alpha?
    {

    # check if there is crucial missing data
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))

    .data <- .data %>% # rearrange by x variable. Use time var to break ties.
        dplyr::arrange(.data[[.x]], .data[[time]])

    lm_poly <- loop_d2_inflection(.data = .data, .x = .x, .y = .y,
                                  degree = degree,
                                  alpha_linearity = alpha_linearity)

    poly_expr <- expr_from_coefs(lm_poly$coefficients)
    deriv1 <- Deriv::Deriv(poly_expr, x = "x", nderiv = 1) # slope
    deriv2 <- Deriv::Deriv(poly_expr, x = "x", nderiv = 2) # acceleration

    optimize(f = expr_to_func(deriv1),
             interval = c(min(.data[[.x]]), max(.data[[.x]])),
             maximum = TRUE)

    # find x-values of roots of second derivative (inflection points)
    inflection_points <- deriv2 %>%
        Ryacas::y_fn("Simplify") %>% # simplify derivative for ease of extraction coefficients
        Ryacas::yac_str() %>%
        stringr::str_extract_all("(?<!\\^)-?\\d+\\.?\\d*e?-?\\d*") %>% # extract coefficients
        purrr::map(as.numeric) %>%
        unlist() %>%
        rev() %>% # reverse order for polyroot()
        polyroot() %>%
        find_real()

    # browser()

    # filter by roots within range of x values
    inflection_points <- inflection_points[inflection_points >= min(.data[[.x]]) &
                                               inflection_points <= max(.data[[.x]])]

    concavity_chgs <- concavity_changes(poly_expr,
                                        inflection_points = inflection_points,
                                        tol = 0.1)
    # filter by concave up to concave down
    inflection_points <- inflection_points[concavity_chgs == "up to down"]

    # what do we do if there's more than one inflection point?
    # do we take the one with the highest 1st derivative slope? But then,
    # how different is that from just finding the highest slope?
    threshold_idx <- which.min(abs(.data[[.x]] - inflection_points))

    # get values at threshold
    bp_dat <- .data[threshold_idx,] %>%
        dplyr::mutate(bp = bp,
                      algorithm = "d2_inflection",
                      x_var = .x,
                      y_var = .y,
                      # determinant_bp = determinant_bp,
                      # pct_slope_change = pct_slope_change,
                      # f_stat = f_stat,
                      # p_val_f = pf_two,
        ) %>%
        dplyr::relocate(bp, algorithm, x_var, y_var,
                        # determinant_bp,
                        # pct_slope_change, f_stat, p_val_f
        )
    bp_dat

}

loop_d2_inflection <- function(.data, .x, .y,
                               degree = NULL, alpha_linearity = 0.05) {
    # if the user specifies a degree, find that and be done with it
    if (!is.null(degree)) {
        lm_poly <- paste0(.y, " ~ ", "1 + ",
                          "poly(", .x, ", degree = ", degree, ", raw = TRUE)") %>%
            as.formula() %>%
            stats::lm(data = .data)
        # if the user does NOT specify a degree, find the best degree using
        # likelihood ratio test
    } else {
        degree = 5 # from testing and previous papers it seems like you need to
        # force a higher derivative if you want a local maxima within the
        # range of x values. That doesn't feel wonderful.
        lm_list = vector(mode = "list", length = 0) # hold lm data
        # keep raw = TRUE for now b/c it's way easier for me to test
        # if the derivative expressions are working
        # starting with degree = 5 b/c that's the minimum I've seen in
        # other papers that seemed to use polynomial fits

        lm_poly <- paste0(.y, " ~ ", "1 + ",
                          "poly(", .x, ", degree = ", degree, ", raw = TRUE)") %>%
            as.formula() %>%
            stats::lm(data = .data)

        lm_list <- append(lm_list, list(lm_poly))

        cont <- TRUE
        i <- 1 # start at 2 b/c we already made linear (degree = 1) model
        while(cont == TRUE) {
            lm_poly <- paste0(.y, " ~ ", "1 + ",
                              "poly(", .x, ", degree = ", degree + i, ", raw = TRUE)") %>%
                as.formula() %>%
                stats::lm(data = .data)
            lm_list <- append(lm_list, list(lm_poly))
            lrt <- stats::anova(lm_list[[i]], lm_list[[i+1]])
            if (is.na(lrt$`Pr(>F)`[2]) | lrt$`Pr(>F)`[2] >= alpha_linearity) {
                cont = FALSE
                lm_poly <- lm_list[[i]] # take the previous model
            }
            i <- i + 1
        }
    }

    lm_poly
}

concavity_changes <- function(f, inflection_points, tol = 0.1) {
    # find second derivative
    f_dd <- Deriv::Deriv(f, x = "x", nderiv = 2)
    if(is.expression(f_dd)) { # expression to function if needed
        f_dd <- expr_to_func(f_dd)
    }
    left_sign <- sign(f_dd(inflection_points - tol))
    right_sign <- sign(f_dd(inflection_points + tol))
    conc_changes <- if_else(left_sign > 0 & right_sign < 0,
                            "up to down", "down to up")
}


.x <- "time"
.y <- "ve"
time <- "time"
.data <- df_avg
# .data <- df_avg[1:threshold_idx,]
degree <- 4
alpha_linearity = 0.05
bp = "vt2"


bp_dat <- d2_inflection(.data = df_avg, .x = "time", .y = "ve", bp = "vt2")
bp_dat
bp_dat$vo2_kg / max(df_avg$vo2_kg)


plot_data <- .data %>%
    select(all_of(.x)) %>%
    mutate(y_hat = lm_poly$fitted.values,
           y_hat_deriv1 = eval(deriv1,
                               envir = list(x = .data[[.x]])),
           y_hat_deriv2 = eval(deriv2,
                               envir = list(x = .data[[.x]])),
    )
plot_data

.data %>% # rearrange by x variable. Use time var to break ties.
    dplyr::arrange(.data[[.x]], .data[[time]]) %>%
    ggplot(aes(x = .data[[.x]], y = .data[[.y]])) +
    geom_point(alpha = 0.5, color = "orange") +
    geom_line(aes(y = lm_poly$fitted.values)) +
    geom_vline(xintercept = inflection_points, linetype = "dotted") +
    # geom_line(alpha = 0.5) +
    # ylim(c(0, 175)) +
    theme_minimal()

ggplot(data = .data, aes(x = .data[[.x]])) +
    geom_line(data = plot_data, aes(y = y_hat_deriv1),
              linetype = "dashed") +
    geom_line(data = plot_data, aes(y = y_hat_deriv2),
              linetype = "dotted") +
    geom_vline(xintercept = inflection_points, linetype = "dotted") +
    geom_hline(yintercept = 0) +
    # ylim(c(-0.005, 0.0051)) +
    theme_minimal()

ggplot(data = df_avg, aes(x = time)) +
    geom_line(aes(y = ve_vo2), color = "purple") +
    geom_line(aes(y = ve_vco2), color = "green") +
    geom_vline(xintercept = inflection_points, linetype = "dotted") +
    theme_minimal()

ggplot(data = df_avg, aes(x = vco2, y = ve)) +
    geom_point(color = "orange", alpha = 0.5) +
    geom_vline(xintercept = bp_dat$vco2) +
    theme_minimal()
