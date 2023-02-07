library(gasExchangeR)
library(tidyverse)
library(janitor)
library(broom)
library(Deriv)
library(Ryacas)
library(readxl)
library(devtools)

# for derivative methods, do users need to choose how they want the function
# created, and then what derivative method to use? e.g. they say they want
# to use a polynomial regression, regression splines, or smoothing splines.
# then, they say if they want to use 2nd derivative maxima, 1st derivative crossing
# or 2nd derivative inflection

# d1 crossing is most similar to Wisén and Wohlfart (2004)
# this is a VT1 method ONLY
# but, could this be applied to VE vs. VCO2? probably not b/c they aren't on the same scale
# but what if we scaled them?
# so I ran through this and it almost seemed like it worked, but the very end
# of the VCO2 vs. time curve saw a rapid drop in it's slope after a sudden rise
# towards the end. The previous final crossing seemed to come a little to early, however.
# Also, just getting them to the same units (mL/min), doesn't fix that either
# because their scales are so different

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

# plot VO2 vs. time data
ggplot(data = df_avg, aes(x = time, y = vo2)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_line(alpha = 0.5) +
    geom_point(aes(y = vco2), alpha = 0.5, color = "blue") +
    geom_line(aes(y = vco2), alpha = 0.5) +
    theme_minimal()

ggplot(data = df_avg, aes(x = time)) +
    geom_point(aes(y = ve_vo2), alpha = 0.5, color = "purple") +
    geom_point(aes(y = ve_vco2), alpha = 0.5, color = "green") +
    theme_minimal()

vt1_dat <- d1_crossing(df_avg, .x = "time", .y = "vo2", bp = "vt1")
# vt1_dat

ggplot(data = df_avg, aes(x = time)) +
    geom_point(aes(y = ve_vo2), alpha = 0.5, color = "purple") +
    geom_point(aes(y = ve_vco2), alpha = 0.5, color = "green") +
    geom_vline(xintercept = vt1_dat$time, color = "purple") +
    theme_minimal()

est_speed_vt1 <- df_avg %>%
    slice(which.min(abs(df_avg$time - (vt1_dat$time - 30)))) %>%
    select(speed) %>%
    pull()

pace_per_mi <- est_speed_vt1^-1*60
# ~ 6:40 min/mi

# find 1st derivative of vo2 vs. time
# find 1st derivative of vco2 vs. time
# find where they cross

# assuming we're using polynomial functions to make the model for the time being
d1_crossing <- function(.data,
                        .x,
                        .y,
                        degree = NULL,
                        vo2 = "vo2",
                        vco2 = "vco2",
                        ve = "ve",
                        time = "time",
                        alpha_linearity = 0.05, # change to just alpha?
                        pos_change = TRUE,
                        bp = "vt1") {

    # check if there is crucial missing data
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))

    .data <- .data %>% # rearrange by x variable. Use time var to break ties.
        dplyr::arrange(.data[[.x]], .data[[time]])

    # best-fit polynomial for vo2
    lm_poly_vo2 <- loop_poly_d1_crossing(.data = .data, .x = time, .y = vo2,
                                         degree = degree,
                                         alpha_linearity = alpha_linearity)

    # best-fit polynomial for vco2
    lm_poly_vco2 <- loop_poly_d1_crossing(.data = .data, .x = time, .y = vco2,
                                          degree = degree,
                                          alpha_linearity = alpha_linearity)
    # 1st derivative for vo2
    poly_expr_vo2 <- expr_from_coefs(lm_poly_vo2$coefficients)
    deriv1_vo2 <- Deriv::Deriv(poly_expr_vo2, x = "x", nderiv = 1) # slope

    # 1st derivative for vco2
    poly_expr_vco2 <- expr_from_coefs(lm_poly_vco2$coefficients)
    deriv1_vco2 <- Deriv::Deriv(poly_expr_vco2, x = "x", nderiv = 1) # slope

    # turn derivative expressions into functions to use with uniroot.all()
    deriv1_vo2_func <- expr_to_func(deriv1_vo2)
    deriv1_vco2_func <- expr_to_func(deriv1_vco2)

    # calculate derivative of difference in 1st derivative functions
    # this lets us know if vco2 was surpassing vo2 or the other way around
    deriv_diff_deriv1 <- paste0(deriv1_vo2 %>%
                                    Ryacas::y_fn("Simplify") %>%
                                    Ryacas::yac_str(),
                                " - (",
                                deriv1_vco2 %>%
                                    Ryacas::y_fn("Simplify") %>%
                                    Ryacas::yac_str(),
                                ")") %>%
        parse(text = .) %>%
        Deriv(x = "x", nderiv = 1)

    roots <- rootSolve::uniroot.all(function(x) deriv1_vo2_func(x) - deriv1_vco2_func(x),
                                    interval = c(min(.data[[.x]]), max(.data[[.x]])))

    # filter by roots within range of x values
    roots <- roots[roots >= min(.data[[.x]]) &
                       roots <= max(.data[[.x]])]

    # filter by which roots have a negative derivative. This indicates that
    # co2 is rising above o2. Finding max finds the last time this occurs.
    final_crossing <- roots[eval(deriv_diff_deriv1, envir = list(x = roots)) < 0] %>%
        max()

    threshold_idx <- which.min(abs(.data[[.x]] - final_crossing))

    # get values at threshold
    bp_dat <- .data[threshold_idx,] %>%
        dplyr::mutate(bp = bp,
                      algorithm = "d1_crossing",
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
    bp_dat # return breakpoint data
}

#' @keywords internal
loop_poly_d1_crossing <- function(.data, .x, .y,
                                  degree = NULL, alpha_linearity = 0.05) {
    # browser()
    # if the user specifies a degree, find that and be done with it
    if (!is.null(degree)) {
        lm_poly <- paste0(.y, " ~ ", "1 + ",
                          "poly(", .x, ", degree = ", degree, ", raw = TRUE)") %>%
            as.formula() %>%
            lm(data = .data)
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
            lm(data = .data)

        lm_list <- append(lm_list, list(lm_poly))

        cont <- TRUE
        i <- 1 # start at 2 b/c we already made linear (degree = 1) model
        while(cont == TRUE) {
            lm_poly <- paste0(.y, " ~ ", "1 + ",
                              "poly(", .x, ", degree = ", degree + i, ", raw = TRUE)") %>%
                as.formula() %>%
                lm(data = .data)
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

.x <- "time"
.y <- "ve_ml"
vo2 <- "vco2"
vco2 <- "ve_ml"
time <- "time"
degree <- NULL
alpha_linearity = 0.05
bp = "vt2"

.data <- df_avg %>%
    mutate(ve_ml = ve * 1000)

bp_dat <- d1_crossing(.data = df_avg, .x = "time", .y = "ve",
                      vo2 = "vco2", vco2 = "ve_ml", time = "time", bp = "vt1")

bp_dat

which.min(abs(.data[["time"]] - bp_dat$time))

df_avg %>%
    mutate(vo2max = max(vo2_kg),
           hr_max = max(hr)) %>%
    slice(which.min(abs(.data[["time"]] - (bp_dat$time - 30)))) %>%
    select(speed, grade, vo2_kg, rq, hr, vo2max, hr_max) %>%
    mutate(pct_vo2max = round(vo2_kg / vo2max * 100,1),
           pct_hr_max = round(hr / hr_max * 100,1)) %>%
    clipr::write_clip()

plot_data <- .data %>%
    select(all_of(.x)) %>%
    mutate(y_hat_vo2 = lm_poly_vo2$fitted.values,
           y_hat_vco2 = lm_poly_vco2$fitted.values,
           y_hat_deriv1_vo2 = eval(deriv1_vo2,
                               envir = list(x = .data[[.x]])),
           y_hat_deriv1_vco2 = eval(deriv1_vco2,
                               envir = list(x = .data[[.x]])),
    )
plot_data

ggplot(data = .data, aes(x = .data[[time]])) +
    geom_point(aes(y = .data[[vo2]]), color = "blue") +
    geom_point(aes(y = .data[[vco2]]), color = "orange")

ggplot(data = .data, aes(x = .data[[.x]])) +
    # geom_point(aes(y = .data[[vo2]]), alpha = 0.5, color = "red")+
    # geom_point(aes(y = .data[[vco2]]), alpha = 0.5, color = "blue") +
    # geom_line(data = plot_data, aes(y = y_hat_vo2)) +
    # geom_line(data = plot_data, aes(y = y_hat_vco2)) +
    geom_line(data = plot_data, aes(y = y_hat_deriv1_vo2), color = "red") +
    geom_line(data = plot_data, aes(y = y_hat_deriv1_vco2), color = "blue") +
    geom_vline(xintercept = bp_dat$time) +
    ggtitle("Derivative Crossing") +
    theme_minimal()

ggplot(data = .data, aes(x = .data[[.x]])) +
    geom_point(aes(y = .data[[vo2]]), alpha = 0.5, color = "red")+
    geom_point(aes(y = .data[[vco2]]), alpha = 0.5, color = "blue") +
    geom_line(data = plot_data, aes(y = y_hat_vo2)) +
    geom_line(data = plot_data, aes(y = y_hat_vco2)) +
    # geom_line(data = plot_data, aes(y = y_hat_deriv1_vo2), color = "red") +
    # geom_line(data = plot_data, aes(y = y_hat_deriv1_vco2), color = "blue") +
    geom_vline(xintercept = bp_dat$time) +
    ggtitle("Derivative Crossing") +
    theme_minimal()

ggplot(data = .data, aes(x = .data[[vo2]])) +
    geom_point(aes(y = .data[[vco2]]), alpha = 0.5, color = "blue") +
    geom_vline(xintercept = bp_dat$vo2) +
    ggtitle("V-slope") +
    theme_minimal()

ggplot(data = .data, aes(x = .data[[time]])) +
    geom_point(aes(y = .data[["ve_vo2"]]), alpha = 0.5, color = "purple") +
    geom_vline(xintercept = bp_dat$time) +
    ggtitle("Ventiliatory Equivalents") +
    theme_minimal()

