library(gasExchangeR)
library(tidyverse)
library(devtools)

file_lines <- readLines("inst/extdata/Anton_vo2max.txt")
df_raw <- read.table(textConnection(file_lines[-2]), header = TRUE, sep="\t")

df_unavg <- df_raw %>%
    as_tibble() %>%
    janitor::clean_names() %>%
    separate(`time`, into = c("m1", "s1"), sep = ":") %>%
    separate(ex_time, into = c("m2", "s2"), sep = ":") %>%
    separate(time_clock,
             into = c("h3", "m3", "s3"),
             sep = ":") %>%
    mutate(across(where(is.character), as.numeric)) %>%
    mutate(time = (m1*60 + s1), .keep = "unused") %>%
    mutate(ex_time = (m2*60 + s2 ), .keep = "unused") %>%
    mutate(clock_time = hms::hms(s3, m3, h3), .keep = "unused") %>%
    relocate(contains("time")) %>%
    filter(!is.na(ex_time)) %>%
    filter(speed >= 4.5 & ex_time >= 750) %>%
    select(-time) %>%
    rename(time = ex_time,
           vo2_kg = vo2,
           vo2 = vo2_1,
           ve = ve_btps) %>%
    mutate(ve_vo2 = ve / vo2 * 1000,
           ve_vco2 = ve/vco2*1000,
           excess_co2 = vco2^2 / vo2 - vco2) %>%
    ventilatory_outliers(plot_outliers = FALSE)

df_avg <- avg_exercise_test(df_unavg, type = "time", subtype = "bin",
                            time_col = "time", bin_w = 10)

wb_rc <- function(.data,
         .x,
         .y,
         bp,
         vo2 = "vo2",
         vco2 = "vco2",
         ve = "ve",
         time = "time",
         alpha_linearity = 0.05,
         front_trim_vt1 = 60,
         front_trim_vt2 = 60,
         pos_change = TRUE,
         min_pct_change = 0.15,
         ...) {

    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))
    bp = match.arg(bp, choices = c("vt1", "vt2"), several.ok = FALSE)

    .data <- .data %>% # rearrange by x variable. Use time var to break ties.
        dplyr::arrange(.data[[.x]], .data[[time]])
    plot_df <- .data

    front_trim <- set_front_trim(bp = bp,
                                 front_trim_vt1 = front_trim_vt1,
                                 front_trim_vt2 = front_trim_vt2)

    .data <- .data %>%
        dplyr::filter(.data[[time]] >= min(.data[[time]] + front_trim))

    loop_idx <- loop_wb_rc(.data = .data, .x = .x, .y = .y,
                      min_pct_change = min_pct_change)

    lm_simple <- stats::lm(paste0(.y, " ~ ", .x), data = .data)

    if(is.na(loop_idx)) {
        bp_dat <- .data[0,] %>%
            dplyr::mutate(algorithm = "wb_rc",
                          x_var = .x,
                          y_var = .y,
                          determinant_bp = FALSE,
                          bp = bp,
                          pct_slope_change = NA,
                          f_stat = NA,
                          p_val_f = NA) %>%
            dplyr::relocate(bp, algorithm, determinant_bp,
                            pct_slope_change, f_stat, p_val_f)
        return(list(breakpoint_data = bp_dat,
                    fitted_vals = pred,
                    lm_left = NA,
                    lm_right = NA,
                    lm_simple = lm_simple,
                    bp_plot = NA))
    } else {
        df_left <- .data[1:loop_idx,]
        df_right <- .data[loop_idx:nrow(.data),]

        lm_left <- stats::lm(paste0(.y, " ~ ", .x), data = df_left)
        lm_right <- stats::lm(paste0(.y, " ~ ", .x), data = df_right)

        pw_stats <- piecewise_stats(lm_left, lm_right, lm_simple)
        list2env(pw_stats, envir = environment())

        pct_slope_change <- 100*(lm_right$coefficients[2] - lm_left$coefficients[2]) /
            abs(lm_left$coefficients[2])

        y_hat_left <- tibble::tibble("{.x}" := df_left[[.x]],
                                     "{.y}" := lm_left$fitted.values,
                                     algorithm = "wb_rc")
        y_hat_right <- tibble::tibble("{.x}" := df_right[[.x]],
                                      "{.y}" := lm_right$fitted.values,
                                      algorithm = "wb_rc")
        pred <- dplyr::bind_rows(y_hat_left, y_hat_right)

        determinant_bp <- dplyr::if_else(pf_two < alpha_linearity &
                                             (pos_change == (pct_slope_change > 0)),
                                         TRUE, FALSE)

        int_point <- intersection_point(lm_left, lm_right)

        # for some reason I'm not able to use get_threshold_vals()
        # bp_dat <- find_threshold_vals(.data = .data, thr_x = int_point["x"],
        #                               int_point["y"], .x = .x,
        #                               .y = .y, ...)
        threshold_idx <- which.min(abs(.data[[.x]] - int_point["x"]))
        bp_dat <- .data[threshold_idx,]

        bp_dat <- bp_dat %>%
            dplyr::mutate(algorithm = "wb_rc",
                          x_var = .x,
                          y_var = .y,
                          determinant_bp = determinant_bp,
                          bp = bp,
                          pct_slope_change = pct_slope_change,
                          f_stat = f_stat,
                          p_val_f = pf_two) %>%
            dplyr::relocate(bp, algorithm, determinant_bp,
                            pct_slope_change, f_stat, p_val_f)

        bp_plot <- make_piecewise_bp_plot(plot_df, .x, .y, lm_left, lm_right, bp_dat)

        return(list(breakpoint_data = bp_dat,
                    fitted_vals = pred,
                    lm_left = lm_left,
                    lm_right = lm_right,
                    lm_simple = lm_simple,
                    bp_plot = bp_plot))
    }
}

loop_wb_rc <- function(.data, .x, .y, min_pct_change = 0.15) {

    # fit a simple linear model with all the data for a comparison
    lm_simple <- lm(paste0(.y, " ~ ", .x), data = .data)
    RSS_simple <- sum(stats::resid(lm_simple)^2)
    df2 <- nrow(lm_simple$model) - 4 # -4 b/c estimating 4 parameters

    # initialize results storage
    f_stats <- numeric(length = nrow(.data))
    p_values <- numeric(length = nrow(.data))
    pct_changes <- numeric(length = nrow(.data))

    # if at the first or last data points, set values to NA. Doing this so
    # the length of the vectors match the number of rows in the data frame
    # this way the index matches up
    f_stats[c(1, length(f_stats))] <- NA
    p_values[c(1, length(p_values))] <- NA
    pct_changes[c(1, length(pct_changes))] <- NA

    # loop through divisions in data and calculate metrics
    for (i in 2:(nrow(.data) - 1)) {
        # split data into left and right portions
        df_left <- .data[1:i, ]
        df_right <- .data[i:nrow(.data), ]

        # fit piecewise linear model to left and right portions using the .x and .y variables
        lm_left <- stats::lm(paste0(.y, " ~ ", .x), data = df_left)
        lm_right <- stats::lm(paste0(.y, " ~ ", .x), data = df_right)

        # calculate and store the pooled sums of squares in the left and right portions
        ss_left <- sum(stats::resid(lm_left)^2)
        ss_right <- sum(stats::resid(lm_right)^2)
        RSS_two <- sum(c(ss_left, ss_right))
        MSE_two <- RSS_two / df2
        # compare the sums of squares in the simple linear model against
        # the pooled sums of squares in the piecewise model using an F test
        f_stat <- (RSS_simple - RSS_two) / (2 * MSE_two)
        p_value <- stats::pf(f_stat, df1 = 2, df2 = df2,
                                lower.tail = FALSE)
        # calculate and store the percent change from the left to the right regression
        pct_change <- 100 * (stats::coefficients(lm_right)[2] - stats::coefficients(lm_left)[2]) /
            abs(stats::coefficients(lm_left)[2])

        # store values
        f_stats[i] <- f_stat
        p_values[i] <- p_value
        pct_changes[i] <- pct_change
    }
    browser()
    # return the index of the division of the data where there was a significant departure from linearity
    # and also a percent change equal to or greater than min_pct_change
    if (any(p_values < 0.05) & any(pct_changes >= min_pct_change)) {
        return(which(p_values < 0.05 & pct_changes == max(pct_changes,
                                                          na.rm = TRUE)))
    } else {
        return(NA)
    }
}

wb_rc(.data = df_avg, .y = "ve", .x = "vco2", bp = "vt2")
