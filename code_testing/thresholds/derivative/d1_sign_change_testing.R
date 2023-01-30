library(gasExchangeR)
library(tidyverse)
library(janitor)
library(broom)
library(Deriv)
library(Ryacas)
library(readxl)
library(devtools)

# d1 sign change is most similar to Wisén and Wohlfart (2004)
# This function finds the first derivative of the VE/VCO2 vs. time
# and notes VT2 as the LAST time the first derivative crosses above zero
# The same procedure can be used to find VT1, but use VE/VO2 vs. time instead
# We're also adding a criteria that the derivative needs to be INCREASING.
# The reason is that the at the very end of the test, there can be a large
# decrease in the fitted (polynomial) equation

# my take on this method is that it has trouble when there are slow and steady
# rises in VE/VO2 or VE/VCO2 over time. They rely on those variables decreasing
# or staying steady over time before a marked rise. However, there is sometimes
# an initial dip in ventilatory equivalents vs. time at the onset of exercise prior
# to a period of relative stability. However, that period of relative stability
# may actually be a very slow, but consistent rise. Given that, the only time
# the first derivative crosses above zero is after the initial onset of exercise
# dip, and not the slow but steady rise prior to the actual threshold. Nick's
# test is a good example of this.

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

ggplot(data = df_avg, aes(x = time)) +
    geom_point(aes(y = ve_vo2), alpha = 0.5, color = "purple") +
    geom_point(aes(y = ve_vco2), alpha = 0.5, color = "green") +
    ggtitle("Ventiliatory Equivalents") +
    ylab("mm Hg") +
    xlab("Time (s)") +
    theme_minimal()


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
        degree = 5 # start at degree = 5 to mimic Wisén and Wolhfart, who found degree = 6 was best
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

d1_sign_change <- function(.data,
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
    lm_poly <- loop_poly_reg(.data = .data, .x = .x, .y = .y,
                                 degree = degree,
                                 alpha_linearity = alpha_linearity)
    # 1st derivative
    poly_expr <- expr_from_coefs(lm_poly$coefficients)
    deriv1 <- Deriv(poly_expr, x = "x", nderiv = 1) # slope
    deriv2 <- Deriv(poly_expr, x = "x", nderiv = 2) # acceleration
    # find x-values of roots
    roots_deriv1 <- deriv1 %>%
        y_fn("Simplify") %>% # simplify derivative for ease of extraction coefficients
        yac_str() %>%
        str_extract_all("(?<!\\^)-?\\d+\\.?\\d*e?-?\\d*") %>% # extract coefficients
        map(as.numeric) %>%
        unlist() %>%
        rev() %>% # reverse order for polyroot()
        polyroot() %>%
        find_real_roots()

    # filter by roots within range of x values
    roots_deriv1 <- roots_deriv1[roots_deriv1 >= min(.data[[.x]]) &
                     roots_deriv1 <= max(.data[[.x]])]

    # filter by positive slopes in the 1st derivative (going from below to above 0)
    final_crossing <- roots_deriv1[eval(deriv2, envir = list(x = roots_deriv1)) > 0] %>%
        max()

    threshold_idx <- which.min(abs(.data[[.x]] - final_crossing))

    # get values at threshold
    bp_dat <- .data[threshold_idx,] %>%
        mutate(bp = bp,
               algorithm = "d1_sign_change",
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
.y <- "ve_vo2"
vo2 <- "vo2"
vco2 <- "vco2"
time <- "time"
degree <- NULL
alpha_linearity = 0.05
bp = "vt1"

bp_dat <- d1_sign_change(.data = df_avg, .x = "time", .y = "ve_vo2",
                         vo2 = "vo2", vco2 = "vco2", time = "time",
                         degree = NULL, bp = "vt1")

bp_dat

plot_data <- .data %>%
    select(all_of(.x)) %>%
    mutate(y_hat = lm_poly$fitted.values,
           y_hat_deriv1 = eval(deriv1,
                                   envir = list(x = .data[[.x]])),
           y_hat_deriv2 = eval(deriv2,
                               envir = list(x = .data[[.x]])),
    )
plot_data

ggplot(data = .data, aes(x = .data[[.x]])) +
    geom_point(aes(y = .data[[.y]])) +
    geom_line(data = plot_data, aes(y = y_hat)) +
    geom_vline(xintercept = bp_dat$time)

ggplot(data = .data, aes(x = .data[[.x]])) +
    geom_line(data = plot_data, aes(y = y_hat_deriv1),
              linetype = "dashed") +
    geom_vline(xintercept = bp_dat$time)

ggplot(data = .data, aes(x = .data[[.x]])) +
    geom_line(data = plot_data, aes(y = y_hat_deriv2),
              linetype = "dotted")
