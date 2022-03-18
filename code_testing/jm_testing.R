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

df_avg <- avg_exercise_test(df_unavg, type = "breath", subtype = "rolling",
                            time_col = "time", roll_window = 9, roll_trim = 4)

.data <- df_avg
.x <- "vo2_abs"
.y <- "vco2"
x_vt2 <- "vco2"
y_vt2 <- "ve"

jm(.data, x_vt2, y_vt2, vo2 = "vo2_abs")

breakpoint(.data = df_avg,
           x_vt1 = "vo2_abs",
           y_vt1 = "vco2",
           algorithm_vt1 = "jm",
           x_vt2 = "vco2",
           y_vt2 = "ve",
           algorithm_vt2 = "jm",
           vo2 = "vo2_abs",
           bps = "both")




####### test notes from breakpoints.R

ss <- loop_jm(.data = .data, .x = .x, .y = .y)
min_ss_idx <- which.min(ss)

jm_dat <- jm(.data = .data, .x = .x, .y = .y, vo2 = "vo2_abs")

jm_dat

ggplot(data = .data, aes(x = .data[[.x]], y = .data[[.y]])) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    geom_point(data = jm_dat$lm_left$model,
               aes(x = df_left[[.x]], y = df_left[[.y]]),
               color = "orange", alpha = 0.5) +
    geom_smooth(data = jm_dat$lm_left$model,
                aes(x = df_left[[.x]], y = df_left[[.y]]),
                method = "lm",
                color = "orange", alpha = 0.5) +
    geom_point(data = jm_dat$lm_right$model,
               aes(x = df_right[[.x]], y = df_right[[.y]]),
               color = "red", alpha = 0.5) +
    geom_smooth(data = jm_dat$lm_right$model,
                aes(x = df_right[[.x]], y = df_right[[.y]]),
                method = "lm", se = FALSE,
                color = "red", alpha = 0.5) +
    geom_vline(xintercept = jm_dat$breakpoint_data$vco2)
