library(gasExchangeR)
library(tidyverse)
library(devtools)
library(lubridate)

df_raw <- read_csv("inst/extdata/37818_vo2max_unavg.csv",
                   show_col_types = FALSE)
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
           ve_vco2 = ve / (vco2/1000),
           time = as.numeric(ms(str_remove(as.character(time), ":00"))))

df_avg <- avg_exercise_test(df_unavg,
                            type = "breath",
                            subtype = "rolling",
                            roll_window = 15,
                            roll_trim = 2,
                            time_col = "time")

ggplot(data = df_unavg, aes(x = vco2, y = ve)) +
    geom_point() +
    geom_point(data = df_avg, aes(x = vco2, y = ve), color = "red")

ggplot(data = df_unavg, aes(x = as.numeric(time), y = ve)) +
    geom_point() +
    geom_point(data = df_avg, aes(x = time, y = ve), color = "red")

ggplot(data = df_unavg, aes(x = as.numeric(time), y = vo2_abs)) +
    geom_point() +
    geom_point(data = df_avg, aes(x = time, y = vo2_abs), color = "red")


.data <- df_avg
.x <- "vo2_abs"
.y <- "vco2"
x_vt2 <- "vco2"
y_vt2 <- "ve"

# breakpoint(.data = df_avg,
#            x_vt1 = "vo2_abs",
#            y_vt1 = "vco2",
#            algorithm_vt1 = "jm",
#            x_vt2 = "vco2",
#            y_vt2 = "ve",
#            algorithm_vt2 = "jm",
#            vo2 = "vo2_abs",
#            bps = "both")

bp_dat <- jm(.data = df_avg, x_vt2, y_vt2, vo2 = "vo2_abs", bp = "vt2")
bp_dat$breakpoint_data %>% View

ggplot(data = bp_dat$fitted_vals, aes(x = vco2, y = ve)) +
    geom_line() +
    geom_vline(xintercept = bp_dat$breakpoint_data$vco2) +
    theme_bw()


ggplot(data = df_avg, aes(x = vco2, y = ve)) +
    geom_point(alpha = 0.5, color = "brown") +
    geom_line(data = bp_dat$fitted_vals, aes(x = vco2, y = ve)) +
    theme_bw() +
    geom_vline(xintercept = bp_dat$breakpoint_data$vco2) +
    ggtitle()





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
