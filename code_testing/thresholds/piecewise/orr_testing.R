library(gasExchangeR)
library(tidyverse)
library(devtools)

df_raw <- read_csv("inst/extdata/mar16_101_pre.csv", show_col_types = FALSE)

df_unavg <- df_raw %>%
    rename_all(.funs = tolower) %>%
    rename(vo2_rel = vo2,
           vo2_abs = vo2.1,
           ve = ve.btps) %>%
    trim_pre_post(intensity_col = "speed") %>%
    trim_pre_post(intensity_col = "speed", pre_ex_intensity = 3) %>%
    trim_pre_post(intensity_col = "speed", pre_ex_intensity = 5.6) %>%
    select(time, speed, grade, vo2_rel, vo2_abs, vco2, ve, peto2, petco2) %>%
    mutate(ve_vo2 = ve*1000 / vo2_abs,
           ve_vco2 = ve*1000 / vco2)

df_avg <- avg_exercise_test(.data = df_unavg,
                            type = "breath",
                            subtype = "rolling",
                            time_col = "time",
                            align = "center",
                            mos = "mean",
                            roll_window = 9,
                            roll_trim = 4)

breakpoint(.data = df_avg,
           x_vt1 = "vo2_abs",
           y_vt1 = "vco2",
           algorithm_vt1 = "orr",
           x_vt2 = "vco2",
           y_vt2 = "ve",
           algorithm_vt2 = "orr",
           vo2 = "vo2_abs",
           vco2 = "vco2",
           ve = "ve",
           time = "time",
           bps = "both")

orr_dat <- orr(.data = df_avg, .x = "vo2_abs", .y = "vco2",
    vo2 = "vo2_abs", vco2 = "vco2", ve = "ve",
    time = "time", alpha_linearity = 0.05, bp = "vt1")
orr_dat


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
