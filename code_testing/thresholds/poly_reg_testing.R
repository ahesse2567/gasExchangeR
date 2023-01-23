library(gasExchangeR)
library(tidyverse)
# library(devtools)
library(janitor)
library(broom)
library(Deriv)
library(Ryacas)

# 10.1249/MSS.0000000000001226 this articles finds the local maxima
# of the second deriative. They want to know where the relationship is
# accelerating suddenly, so they want the second derivative of the
# original curve. However, since they want the maxima, we need to take
# another (third) derivative.

# big thanks to this stack overflow thread: https://stackoverflow.com/questions/40438195/function-for-polynomials-of-arbitrary-order-symbolic-method-preferred

# df_raw <- read_delim("inst/extdata/Aaron_VO2max.txt",
#                      delim = "\t",
#                      na = c("", "NA", " "),
#                      show_col_types = FALSE)

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
           ve_vco2 = ve*1000 / vco2) %>%
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

# plot averaged VO2 vs. time data
ggplot(data = df_avg, aes(x = time, y = vo2_rel)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_line(alpha = 0.5) +
    theme_minimal() +
    ggtitle("Averaged Data")

ggplot(data = df_avg, aes(x = vo2_abs, y = vco2)) +
    geom_point(color = "blue", alpha = 0.5) +
    theme_minimal()

ggplot(data = df_avg, aes(x = time)) +
    geom_point(aes(y = ve_vo2), color = "purple", alpha = 0.5) +
    geom_line(aes(y = ve_vo2), color = "purple", alpha = 0.5) +
    geom_point(aes(y = ve_vco2), color = "green", alpha = 0.5) +
    geom_line(aes(y = ve_vco2), color = "green", alpha = 0.5) +
    theme_minimal() +
    ggtitle("Ventilatory Equivalents")


.x <- "time"
.y <- "ve_vco2"
time <- "time"
.data <- df_avg
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


    lm_poly <- loop_poly_reg(.data = df_avg, .x = .x, .y = .y,
                             degree = NULL,
                             alpha_linearity = alpha_linearity)

    # TODO what about if the best model is linear? Then this method probably
    # won't work well and you should get a warning or an error about that
    # also, I think you need a minimum of a 4th order equation if you want to
    # find the maxima in the acceleration (2nd derivative) of the relationship
    # b/c you essentially need to take a third derivative and still have an x term

    # lm_poly <- lm_list[[3]] # keep for testing, delete later
    poly_expr <- expr_from_coefs(lm_poly$coefficients)
    deriv1 <- Deriv(poly_expr, x = "x", nderiv = 1) # slope
    deriv2 <- Deriv(poly_expr, x = "x", nderiv = 2) # acceleration
    deriv3 <- Deriv(poly_expr, x = "x", nderiv = 3) # jerk
    deriv4 <- Deriv(poly_expr, x = "x", nderiv = 4) # snap. This kinda feels hacky and dumb at this point

    roots_deriv3 <- deriv3 %>%
        y_fn("Simplify") %>% # simplify derivative for ease of extraction coefficients
        yac_str() %>%
        str_extract_all("(?<!\\^)-?\\d+\\.?\\d*e?\\d?") %>% # extract coefficients
        map(as.numeric) %>%
        unlist() %>%
        rev() %>%
        polyroot() %>%
        find_real_roots()

    local_maxima <- roots_deriv3[eval(deriv4, envir = list(x = roots_deriv3)) < 0]
    # you would generally expect the threshold to be the highest of any possible
    # local maxima
    # TODO local maxima need to be within the range of valid x values, too
    threshold_idx <- which.min(abs(.data[[.x]] - max(local_maxima)))

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
    select(time) %>%
    mutate(y_hat = predict(lm_poly),
           y_hat_deriv1 = eval(deriv1,
                                      envir = list(x = df_avg$time)),
           y_hat_deriv2 = eval(deriv2,
                                      envir = list(x = df_avg$time)),
           # y_hat_deriv3 = eval(deriv3,
           #                            envir = list(x = df_avg$time)),
           )
plot_data

# plot to see how the model fits the data
ggplot(data = .data, aes(x = time)) +
    # geom_line(aes(y = ve_vo2), color = "purple", alpha = 0.5) +
    geom_line(aes(y = ve_vco2), color = "green", alpha = 0.5) +
    # scale_color_manual(name = "Ventilatory\nEquivalents",
    #                    values = c("ve_vo2" = "purple", "ve_vco2" = "green"),
    #                    labels = c("VE/VO2", "VE/VCO2")) +
    ylab("mm Hg") +
    xlab("Time (min)") +
    # ylim(c(-2.5,5)) +
    geom_line(data = plot_data, aes(x = time, y = y_hat/max(y_hat))) +
    geom_line(data = plot_data,
              aes(x = time, y = y_hat_deriv1/y_hat_deriv1), linetype = "dashed") +
    geom_line(data = plot_data,
              aes(x = time, y = y_hat_deriv2/y_hat_deriv2), linetype = "dotted") +
    # geom_line(data = plot_data,
    #           aes(x = time, y = y_hat_deriv3/max(abs(y_hat_deriv3))), linetype = "dotdash") +
    geom_vline(xintercept = max(local_maxima)) +
    theme_minimal()

# zoomed in on derivatives
ggplot(data = .data, aes(x = time)) +
    theme_minimal() +
    ylab("mm Hg") +
    xlab("Time (min)") +
    geom_line(data = plot_data,
              aes(x = time, y = y_had_ve_vo2_deriv1), linetype = "dashed") +
    geom_line(data = plot_data,
              aes(x = time, y = y_had_ve_vo2_deriv2), linetype = "dotted")


loop_poly_reg <- function(.data, .x, .y,
                          degree = NULL, alpha_linearity = 0.05) {
    # browser()

    # if the user specifies a degree, find that and be done with it
    if (!is.null(degree)) {
        lm_poly <- lm(.data[[.y]] ~ 1 + poly(.data[[.x]], degree = degree),
                      data = .data)
    # if the user does NOT specify a degree, find the best degree using
    # likelihood ratio test
    } else {
        lm_list = vector(mode = "list", length = 0) # hold lm data
        # keep raw = TRUE for now b/c it's way easier for me to test
        # if the derivative expressions are working
        lm_poly <- lm(.data[[.y]] ~ 1 + poly(.data[[.x]],
                                             degree = 1, raw = TRUE),
                      data = .data)
        lm_list <- append(lm_list, list(lm_poly))

        cont <- TRUE
        i <- 2 # start at 2 b/c we already made linear (degree = 1) model
        while(cont == TRUE) {
            lm_poly <- lm(.data[[.y]] ~ 1 + poly(.data[[.x]],
                                                 degree = i,
                                                 raw = TRUE),
                          data = .data)
            lm_list <- append(lm_list, list(lm_poly))
            lrt <- anova(lm_list[[i-1]], lm_list[[i]])
            if (is.na(lrt$`Pr(>F)`[2]) | lrt$`Pr(>F)`[2] >= alpha_linearity) {
                cont = FALSE
                lm_poly <- lm_list[[i-1]] # take the previous model
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
