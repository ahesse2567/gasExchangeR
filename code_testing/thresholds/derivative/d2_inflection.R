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

# will we need some sort of curve type argument for derivative methods?
# e.g., poly, reg spline, smoothing spline?

# note to Anton: inflection points seem like (almost?) always different from
# maxima. Would that introduce a systematic bias between then?

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
    theme_minimal()

ggplot(data = df_avg, aes(x = time, y = vo2_kg)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_line(alpha = 0.5) +
    theme_minimal() +
    ggtitle("Raw Data")

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
                              "poly(", .x, ", degree = ", degree, ", raw = TRUE)") %>%
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

.x <- "vo2"
.y <- "ve"
time <- "time"
.data <- df_avg
# .data <- df_avg[1:threshold_idx,]
degree <- NULL
alpha_linearity = 0.05
bp = "vt2"

d2_inflection <- function(.data,
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

    poly_expr <- expr_from_coefs(lm_poly$coefficients)
    deriv1 <- Deriv(poly_expr, x = "x", nderiv = 1) # slope
    deriv2 <- Deriv(poly_expr, x = "x", nderiv = 2) # acceleration

    # find x-values of roots
    inflection_points <- deriv2 %>%
        y_fn("Simplify") %>% # simplify derivative for ease of extraction coefficients
        yac_str() %>%
        str_extract_all("(?<!\\^)-?\\d+\\.?\\d*e?-?\\d*") %>% # extract coefficients
        map(as.numeric) %>%
        unlist() %>%
        rev() %>% # reverse order for polyroot()
        polyroot() %>%
        find_real_roots()

    # filter by roots within range of x values
    inflection_points <- inflection_points[inflection_points >= min(.data[[.x]]) &
                                               inflection_points <= max(.data[[.x]])]

    # what do we do if there's more than one inflection point?
    threshold_idx <- which.min(abs(.data[[.x]] - inflection_points))

    # get values at threshold
    bp_dat <- .data[threshold_idx,] %>%
        mutate(bp = bp,
               algorithm = "d2_inflection",
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

bp_dat <- d2_inflection(.data = df_avg, .x = "vo2", .y = "ve", bp = "vt2")
bp_dat$vo2_kg / max(df_avg$vo2_kg)


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

ggplot(data = .data, aes(x = vo2, y = ve)) +
    geom_point(alpha = 0.5, color = "orange") +
    geom_line(aes(y = lm_poly$fitted.values)) +
    geom_vline(xintercept = inflection_points) +
    # geom_line(alpha = 0.5) +
    # ylim(c(0, 175)) +
    theme_minimal()

ggplot(data = .data, aes(x = .data[[.x]])) +
    # geom_line(data = plot_data, aes(y = y_hat_deriv1),
    #           linetype = "dashed") +
    geom_line(data = plot_data, aes(y = y_hat_deriv2),
              linetype = "dotted") +
    geom_vline(xintercept = inflection_points) +
    theme_minimal()

ggplot(data = df_avg, aes(x = time)) +
    geom_line(aes(y = ve_vo2), color = "purple") +
    geom_line(aes(y = ve_vco2), color = "green") +
    geom_vline(xintercept = bp_dat$time) +
    theme_minimal()

ggplot(data = df_avg, aes(x = vco2, y = ve)) +
    geom_point(color = "orange", alpha = 0.5) +
    geom_vline(xintercept = bp_dat$vco2) +
    theme_minimal()
