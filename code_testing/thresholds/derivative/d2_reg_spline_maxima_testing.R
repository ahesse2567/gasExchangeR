library(gasExchangeR)
library(tidyverse)
library(janitor)
library(broom)
library(splines)
library(zoo)
library(devtools)
library(readxl)

# big thanks to this stack overflow post
# https://stackoverflow.com/questions/44739192/export-fitted-regression-splines-constructed-by-bs-or-ns-as-piecewise-poly

# maybe try this one?
# https://stackoverflow.com/questions/29499686/how-to-extract-the-underlying-coefficients-from-fitting-a-linear-b-spline-regres

# details
# This forces a minimum of one knot, uses a b-spline basis, uses splinefun
# to generate derivatives of best-fit spline

df_unavg <- read_xlsx("../gasExchangeR_validation/data/processed/rand_15_cpet_exercisethresholds/mar22_139_pre_gxt.xlsx") %>%
    rename(vo2 = vo2_abs,
           vo2_kg = vo2)


# file_lines <- readLines("inst/extdata/Anton_vo2max.txt")
# df_raw <- read.table(textConnection(file_lines[-2]), header = TRUE, sep="\t")
#
# df_unavg <- df_raw %>%
#     as_tibble() %>%
#     clean_names() %>%
#     separate(`time`, into = c("m1", "s1"), sep = ":") %>%
#     separate(ex_time, into = c("m2", "s2"), sep = ":") %>%
#     separate(time_clock,
#              into = c("h3", "m3", "s3"),
#              sep = ":") %>%
#     mutate(across(where(is.character), as.numeric)) %>%
#     mutate(time = (m1*60 + s1), .keep = "unused") %>%
#     mutate(ex_time = (m2*60 + s2 ), .keep = "unused") %>%
#     mutate(clock_time = hms::hms(s3, m3, h3), .keep = "unused") %>%
#     relocate(contains("time")) %>%
#     filter(!is.na(ex_time)) %>%
#     filter(speed >= 4.5 & ex_time >= 750) %>%
#     select(-time) %>%
#     rename(time = ex_time,
#            vo2_kg = vo2,
#            vo2 = vo2_1,
#            ve = ve_btps) %>%
#     mutate(ve_vo2 = ve / vo2 * 1000,
#            ve_vco2 = ve/vco2*1000,
#            excess_co2 = vco2^2 / vo2 - vco2) %>%
#     ventilatory_outliers(plot_outliers = FALSE)

# df_colnames <- read_xlsx("inst/extdata/Foreman_Ramp_10_18_2022.xlsx",
#                          n_max = 0) %>%
#     colnames()
# df_raw <- read_xlsx("inst/extdata/Foreman_Ramp_10_18_2022.xlsx",
#                     col_names = FALSE, skip = 3) %>%
#     set_names(df_colnames)
#
# df_unavg <- df_raw %>%
#     clean_names() %>%
#     rename(time = t) %>%
#     mutate(time = lubridate::minute(time) * 60 + lubridate::second(time)) %>%
#     filter(phase == "EXERCISE") %>%
#     relocate(time, speed, grade) %>%
#     ventilatory_outliers()
#
df_avg <- avg_exercise_test(df_unavg, type = "time", subtype = "bin",
                            time_col = "time", bin_w = 15)

ggplot(data = df_avg, aes(x = time, y = vo2)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_line(alpha = 0.5) +
    geom_point(aes(y = vco2), alpha = 0.5, color = "blue") +
    geom_line(aes(y = vco2), alpha = 0.5) +
    theme_minimal() +
    ggtitle("Averaged Data")

bp_dat <- d2_reg_spline_maxima(.data = df_avg,
                               .x = "vo2",
                               .y = "ve_vco2",
                               bp = "vt2",
                               degree = 5)

bp_dat$bp_plot

.data <- df_avg
.x <- "vo2" # Leo et al and most others specifically used VO2, not time
.y <- "ve_vco2"
bp <- "vt2"
time <- "time"
degree <- 5
# df <- 16 + 3
alpha_linearity = 0.05
df <- NULL
pos_change <- TRUE

.data %>%
    arrange(.data[[.x]], .data[[time]]) %>%
    ggplot(aes(x = .data[[.x]])) +
    geom_point(aes(y = .data[[.y]]), alpha = 0.5) +
    theme_minimal() +
    ggtitle("Averaged Data")

d2_reg_spline_maxima <- function(.data,
                                 .x,
                                 .y,
                                 bp,
                                 degree = 3,
                                 df = NULL,
                                 vo2 = "vo2",
                                 vco2 = "vco2",
                                 ve = "ve",
                                 time = "time",
                                 alpha_linearity = 0.05, # change to just alpha?
                                 pos_change = TRUE,
                                 ...) {

    # check if there is crucial missing data
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))

    .data <- .data %>% # rearrange by x variable. Use time var to break ties.
        dplyr::arrange(.data[[.x]], .data[[time]])

    # find best number of knots
    lm_spline <- loop_d2_reg_spline(.data = .data, .x = .x, .y = .y,
                                    df = df, degree = degree)

    # It feels like there's a better way than splinefun to find a more exact
    # 2nd derivative, but I don't know it

    # find maxima in 2nd derivative. This is accomplished by making
    # a spline interpolation function. Use this with optimize() function.

    # get new values at equal spacing for a smoother splinefun result
    equi_spaced_x <- seq(from = min(.data[[.x]]), to = max(.data[[.x]]),
                         length.out = nrow(.data))
    pred <- predict(lm_spline, newdata = tibble("{.x}" := equi_spaced_x))

    spline_func <- splinefun(x = equi_spaced_x, y = pred)
    # find number of maxima
    sign_change_idx <- slope_sign_changes(y = spline_func(x = equi_spaced_x,
                                                          deriv = 2),
                                          change = "pos_to_neg")
    # filter by expected slope (usually positive)
    # there can be a spike in accel. during the initial drop in vent eqs
    # pos_change &
    if(pos_change) {
        sign_change_idx <- sign_change_idx[spline_func(
            x = .data[[.x]][sign_change_idx],
            deriv = 1) > 0]
    } else {
        sign_change_idx <- sign_change_idx[spline_func(
            x = .data[[.x]][sign_change_idx],
            deriv = 1) < 0]
    }

    extrema_deriv2 <- vector(mode = "list", length = length(sign_change_idx))
    for(i in seq_along(sign_change_idx)) {
        interval <- c(.data[[.x]][sign_change_idx[i]-1],
                      .data[[.x]][sign_change_idx[i]+1])

        extrema_deriv2[[i]] <- optimize(f = spline_func,
                                  interval = interval,
                                  deriv = 2, maximum = TRUE)
    }
    # find highest maxima
    threshold <- purrr::map(extrema_deriv2, function(x) x["objective"]) %>%
        unlist() %>%
        which.max() %>%
        {extrema_deriv2[[.]]}
    threshold_x <- threshold$maximum
    threshold_y <- predict(lm_spline, newdata = tibble("{.x}" := threshold_x))

    # use highest maxima as threshold
    threshold_idx <- sign_change_idx[which.max(y_val_sign_changes)]

    # x_threshold <- .data[[.x]][threshold_idx]
    # y_hat_threshold <- predict(lm_spline,
    #                                tibble("{.x}" := x_threshold))

    # get values at threshold
    bp_dat <- .data[threshold_idx,] %>%
        dplyr::mutate(bp = bp,
                      algorithm = "d2_reg_spline_maxima",
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

    bp_plot <- ggplot(data = .data, aes(x = .data[[.x]], y = .data[[.y]])) +
        geom_point(alpha = 0.5) +
        geom_line(aes(x = equi_spaced_x, y = pred)) +
        geom_vline(xintercept = bp_dat[[.x]]) +
        theme_minimal()

    return(list(bp_dat = bp_dat,
                lm_reg_spline = lm_spline,
                spline_deriv_func = spline_func,
                bp_plot = bp_plot))
}

loop_d2_reg_spline <- function(.data, .x, .y, df = NULL,
                               degree = 3, alpha_linearity = 0.05) {
    # TODO allow users to specify b-spline or natural-spline basis
    # would that use do.call()?

    # if statement for if users specify the knots or df
    if(!is.null(df)) {
        lm_spline <- paste0(.y, " ~ ", "1 + ",
                            "bs(", .x,
                            ", df = ", df,
                            ", degree = ", degree, ")") %>%
            as.formula() %>%
            lm(data = .data)
    } else {
        spline_mod_list = vector(mode = "list", length = 0)
        cont <- TRUE
        i <- 1
        # reference model with one interior knot
        lm_spline <- paste0(.y, " ~ ", "1 + ",
                            "bs(", .x,
                            ", df = ", i + degree,
                            ", degree = ", degree, ")") %>%
            as.formula() %>%
            lm(data = .data)
        spline_mod_list <- append(spline_mod_list, list(lm_spline))
        # while loop beings with 1 knot (3 df already used assuming a 3rd spline)
        while(cont == TRUE) {
            # browser()
            i <- i + 1
            # TODO add options for user-defined knots and df
            lm_spline <- paste0(.y, " ~ ", "1 + ",
                                "bs(", .x,
                                ", df = ", i + degree,
                                ", degree = ", degree, ")") %>%
                as.formula() %>%
                lm(data = .data)
            spline_mod_list <- append(spline_mod_list, list(lm_spline))
            lrt <- anova(spline_mod_list[[i-1]], spline_mod_list[[i]])
            if (is.na(lrt$`Pr(>F)`[2]) | lrt$`Pr(>F)`[2] >= alpha_linearity) {
                cont = FALSE
                lm_spline <- spline_mod_list[[i-1]] # take the previous model
            }
        }
    }
    lm_spline
}

bp_dat <- d2_reg_spline_maxima(.data = df_avg, .x = "vo2", .y = "ve_vco2",
                               bp = "vt2", degree = 5)
bp_dat$bp_plot

.data %>%
    arrange(.data[[.x]], .data[[.y]]) %>%
    ggplot(aes(x = .data[[.x]], y = .data[[.y]])) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = bp_dat$lm_reg_spline$fitted.values)) +
    geom_vline(xintercept = bp_dat$bp_dat[[.x]]) +
    geom_text(aes(label = time), vjust = 0) +
    # geom_vline(xintercept = .data[[.x]][threshold_idx]) +
    # geom_vline(xintercept = .data[[.x]][sign_change_idx]) +
    ggtitle("Fitted Data and threshold") +
    theme_minimal()

.data %>%
    arrange(.data[[time]], .data[[.y]]) %>%
    ggplot(aes(x = .data[[time]], y = .data[[.y]])) +
    geom_point(alpha = 0.5) +
    # geom_line(aes(x = .data[[time]],
    #               y = bp_dat$lm_reg_spline$fitted.values)) +
    geom_vline(xintercept = bp_dat$bp_dat[[time]]) +
    # geom_vline(xintercept = .data[[.x]][threshold_idx]) +
    # geom_vline(xintercept = .data[[.x]][sign_change_idx]) +
    ggtitle("Fitted Data and threshold") +
    theme_minimal()

.data %>%
    arrange(.data[[.x]], .data[[.y]]) %>%
    ggplot(aes(x = .data[[.x]], y = .data[[.y]])) +
    geom_line(aes(x = equi_spaced_x,
                  y = bp_dat$spline_deriv_func(x = equi_spaced_x, deriv = 2)),
              linetype = "dotted") +
    # geom_vline(xintercept = .data[[.x]][threshold_idx]) +
    geom_vline(xintercept = bp_dat$bp_dat[[.x]]) +
    ggtitle("2nd derivative maxima") +
    theme_minimal()

.data %>%
    arrange(.data[[.x]], .data[[.y]]) %>%
    ggplot(aes(x = .data[[.x]], y = .data[[.y]])) +
    geom_line(aes(y = bp_dat$spline_deriv_func(x = .data[[.x]], deriv = 1)),
              linetype = "dotted") +
    # geom_vline(xintercept = .data[[.x]][threshold_idx]) +
    geom_vline(xintercept = bp_dat$bp_dat[[.x]]) +
    ggtitle("1st derivative") +
    # ylim(c(-0.001, 0.001)) +
    theme_minimal()


# VIM or mice packages looks like they could be good for imputation / interpolation
# technically



x_threshold_idx <- .data[[.x]][threshold_idx]
y_hat_threshold_idx <- predict(lm_spline,
                               tibble("{.x}" := x_threshold_idx))
threshold_point <- tibble("{.x}" := x_threshold_idx,
                          "{.y}" := y_hat_threshold_idx)
# calculate Euclidean distance to threshold point in case the associated
# y value is very high or very low

# do I need to normalize everything?

d <- threshold_point %>%
    rbind(.data %>%
              select(all_of(c(.x, .y)))) %>%
    # normalize %>%
    dist()
n <- attr(d, "Size")
threshold_idx <- which.min(d[dist_idx_conv(2:n, 1, d)])

# future Anton can implement a better way to extract distance
# values than use as.matrix(dist_obj)
# https://stackoverflow.com/questions/39005958/r-how-to-get-row-column-subscripts-of-matched-elements-from-a-distance-matri
