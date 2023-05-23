library(gasExchangeR)
library(tidyverse)
library(devtools)
library(readxl)
library(janitor)

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

df_avg <- avg_exercise_test(df_unavg, method = "time", calc_type = "bin",
                            time_col = "time", bin_w = 10)

ggplot(data = df_avg, aes(x = time, y = vo2)) +
    geom_point(alpha = 0.5) +
    theme_bw()

debug(orr)
breakpoint(.data = df_avg,
           x_vt1 = "vo2",
           y_vt1 = "vco2",
           algorithm_vt1 = "orr",
           x_vt2 = "vco2",
           y_vt2 = "ve",
           algorithm_vt2 = "orr",
           bp = "both",
           vo2 = "vo2",
           vco2 = "vco2",
           ve = "ve",
           time = "time",
           k = 5)

orr_dat <- orr(.data = df_avg, .x = "vo2", .y = "vco2",
               bp = "vt1", vo2 = "vo2", vco2 = "vco2", ve = "ve",
               time = "time", alpha_linearity = 0.05,
               thr_calc_method = "inv_dist", k = 5)
orr_dat$bp_plot


plot_data <- tibble(vo2_abs = seq(min(df_avg$vo2_abs), max(df_avg$vo2_abs), by = 1),
                    vco2_left = seq(min(df_avg$vo2_abs),
                                    max(df_avg$vo2_abs),
                                    by = 1) * orr_dat$lm_left$coefficients[2] +
                        orr_dat$lm_left$coefficients[1],
                    vco2_right = seq(min(df_avg$vo2_abs),
                                    max(df_avg$vo2_abs),
                                    by = 1) * orr_dat$lm_right$coefficients[2] +
                        orr_dat$lm_right$coefficients[1])

min_ss_idx <- orr_dat$lm_left$model %>% nrow()


df_ordered <- df_avg %>%
    arrange(vo2_abs, time)
df_ordered[c("vco2", "vo2_abs")] %>% head()
orr_dat$lm_left$model %>% head()

point <- df_ordered[min_ss_idx,]
fitted_point <- orr_dat$fitted_vals[min_ss_idx,]


ggplot(data = df_ordered, aes(x = vo2_abs, y = vco2)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    geom_line(data = orr_dat$fitted_vals[1:min_ss_idx,], aes(x = vo2_abs, y = vco2),
              linetype = "dashed", color = "red") +
    geom_line(data = orr_dat$fitted_vals[(min_ss_idx+1):nrow(orr_dat$fitted_vals),],
              aes(x = vo2_abs, y = vco2),
              linetype = "dashed") +
    geom_point(data = point, color = "green", size = 8, alpha = 0.5) +
    geom_point(data = fitted_point, color = "blue", size = 3, alpha = 0.5) +
    xlim(2600, 3200) +
    ylim(2400, 3200)


ggplot(data = df_avg, aes(x = vo2_abs, y = vco2)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    geom_line(data = plot_data,
              aes(x = vo2_abs, y = vco2_left),
              color = "red") +
    geom_line(data = plot_data,
              aes(x = vo2_abs, y = vco2_right),
              color = "green") +
    geom_line(data = orr_dat$fitted_vals[1:min_ss_idx,], aes(x = vo2_abs, y = vco2),
              linetype = "dashed", size = 2) +
    geom_line(data = orr_dat$fitted_vals[(min_ss_idx+1):nrow(orr_dat$fitted_vals),],
          aes(x = vo2_abs, y = vco2),
          linetype = "dashed", size = 2) +
    geom_vline(xintercept = orr_dat$breakpoint_data$vo2_abs) +
    geom_point(data = point, color = "green", size = 8, alpha = 0.5) +
    geom_point(data = fitted_point, color = "blue", size = 3, alpha = 0.5)








##### code to vectorize loop_orr()

library(profvis)

calc_ss <- function(df, x, y) {
    splits <- vector(mode = "list", length = nrow(df) - 2)
    for(i in 2:(nrow(df)-1)) {
        splits[[i]] <- list(left = data[1:i,], right = data[(i+1):nrow(data),])
    }
}

splits <- map(2:(nrow(df)-1), function(i) {
    list(left = df[1:i, ], right = df[(i+1):nrow(df), ])
})

profvis({
    df <- tibble(x = 1:10000, y = 2 * x)

    splits <- vector(mode = "list", length = nrow(df))
    for(i in seq_along(splits)) {
        splits[[i]] <- list(df_left = df[1:i,], df_right = df[i:nrow(df),])
    }

    splits2 <- map2(1:nrow(df), df, function(i, df) {
        list(left = df[1:i, ], right = df[i:nrow(df), ])
    })
})

