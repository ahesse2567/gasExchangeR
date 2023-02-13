library(gasExchangeR)
library(tidyverse)
library(devtools)

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
    mutate(time = time * 60,
           ve_vo2 = ve*1000 / vo2_abs,
           ve_vco2 = ve*1000 / vco2,
           excess_co2 = (vco2^2 / vo2_abs) - vco2) %>%
    filter(!is.na(time)) %>%
    filter(stage >=6) %>%
    relocate(time, speed, grade) %>%
    ventilatory_outliers(outlier_cols = "vo2_abs", plot_outliers = TRUE)

df_avg <- avg_exercise_test(.data = df_avg, type = "time", subtype = "bin",
                            bin_w = 10)

ggplot(data = df_avg, aes(x = vco2, y = ve)) +
    geom_point(alpha = 0.5)

bp_dat_dmax <- dmax(.data = df_avg, .x = "vco2", .y = "ve", vo2 = "vo2_abs", bp = "vt2")
bp_dat_dmax$bp_plot
