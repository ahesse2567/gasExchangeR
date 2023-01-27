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

# plot VO2 vs. time data
ggplot(data = df_avg, aes(x = time, y = vo2)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_line(alpha = 0.5) +
    geom_point(aes(y = vco2), alpha = 0.5, color = "blue") +
    geom_line(aes(y = vco2), alpha = 0.5) +
    theme_minimal()

# find 1st derivative of vo2 vs. time
# find 1st derivative of vco2 vs. time
# find where they cross

# assuming we're using polynomial functions to make the model for the time being
loop_poly_reg <- function(.data, .x, .y,
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
        degree = 1 # with this method, I don't think we have to start at a higher degree
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
        i <- 1 # start at 2 (degree + i) b/c we already made linear (degree = 1) model
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

expr_to_func <- function(expr) {
    if(is.character(expr)) {
        expr <- parse(text = expr) # parse() turns a character string into an expression
    }
    function(x) {eval(expr, envir = list(x = x))} # return a function
}


d1_crossing <- function(.data,
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

    # best-fit polynomial for vo2
    lm_poly_vo2 <- loop_poly_reg(.data = .data, .x = .x, .y = vo2,
                             degree = degree,
                             alpha_linearity = alpha_linearity)

    # best-fit polynomial for vco2
    lm_poly_vco2 <- loop_poly_reg(.data = .data, .x = .x, .y = vco2,
                                 degree = degree,
                                 alpha_linearity = alpha_linearity)
    # 1st derivative for vo2
    poly_expr_vo2 <- expr_from_coefs(lm_poly_vo2$coefficients)
    deriv1_vo2 <- Deriv(poly_expr_vo2, x = "x", nderiv = 1) # slope

    # 1st derivative for vco2
    poly_expr_vco2 <- expr_from_coefs(lm_poly_vco2$coefficients)
    deriv1_vco2 <- Deriv(poly_expr_vco2, x = "x", nderiv = 1) # slope

    # turn derivative expressions into functions to use with uniroot.all()
    deriv1_vo2_func <- expr_to_func(deriv1_vo2)
    deriv1_vco2_func <- expr_to_func(deriv1_vco2)

    # calculate derivative of difference in 1st derivative functions
    # this lets us know if vco2 was surpassing vo2 or the other way around
    deriv_diff_deriv1 <- paste0(deriv1_vo2 %>%
                                 y_fn("Simplify") %>%
                                 yac_str(),
                             " - (",
                             deriv1_vco2 %>%
                                 y_fn("Simplify") %>%
                                 yac_str(),
                             ")") %>%
        parse(text = .) %>%
        Deriv(x = "x", nderiv = 1)

    roots <- rootSolve::uniroot.all(function(x) deriv1_vo2_func(x) - deriv1_vco2_func(x),
                           interval = c(min(.data[[.x]]), max(.data[[.x]])))

    # filter by which roots have a negative derivative. This indicates that
    # co2 is rising above o2. Finding max finds the last time this occurs.
    final_crossing <- roots[eval(deriv_diff_deriv1, envir = list(x = roots)) < 0] %>%
        max()

    threshold_idx <- which.min(abs(.data[[.x]] - final_crossing))

    # get values at threshold
    bp_dat <- .data[threshold_idx,] %>%
        mutate(bp = bp,
               algorithm = "d1_crossing",
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
    bp_dat # return breakpoint data
}

.data <- df_avg
.x <- "time"
.y <- "vo2"
vo2 <- "vo2"
vco2 <- "vco2"
time <- "time"
degree <- NULL
alpha_linearity = 0.05
bp = "vt1"

bp_dat <- d1_crossing(.data = df_avg, .x = "time", .y = "vo2", vo2 = "vo2",
            vco2 = "vco2", time = "time", bp = "vt1")

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

