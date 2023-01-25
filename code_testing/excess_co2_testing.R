library(gasExchangeR)
library(tidyverse)
library(devtools)

df_raw <- read_csv("inst/extdata/37818_vo2max_unavg.csv", show_col_types = FALSE)

df_unavg <- df_raw %>%
    rename_with(.fn = tolower) %>%
    rename(vo2_rel = vo2...4,
           vo2_abs = vo2...5,
           ve = `ve btps`,
           vt = `vt btps`) %>%
    select(time, speed, grade, vo2_rel, vo2_abs, vco2, ve) %>%
    trim_pre_post(intensity_col = "grade",
                  pre_ex_intensity = 300.1,
                  post_ex_intensity = 300.1) %>%
    mutate(ve_vo2 = ve / (vo2_abs/1000),
           ve_vco2 = ve / (vco2/1000))

df_avg <- avg_exercise_test(df_unavg, type = "breath", subtype = "bin",
                            time_col = "time", bin_w = 12)

df_avg <- df_avg %>%
    mutate(excess_co2 = (vco2^2 / vo2_abs) - vco2)

ggplot(data = df_avg, aes(x = time, y = excess_co2)) +
    geom_point() +
    theme_bw()

bp_info <- breakpoint(.data = df_avg,
           method = "excess_co2",
           algorithm_vt1 = "jm",
           x_vt2 = "vco2",
           y_vt2 = "ve",
           algorithm_vt2 = "jm",
           vo2 = "vo2_abs",
           bps = "both")

bp_info$bp_dat

ggplot(data = df_avg, aes(x = time, y = excess_co2)) +
    geom_point() +
    theme_bw() +
    geom_vline(xintercept = bp_info$vt1_dat$breakpoint_data$time) +
    geom_vline(xintercept = bp_info$vt2_dat$breakpoint_data$time)


bp_info$bp_dat %>% View
