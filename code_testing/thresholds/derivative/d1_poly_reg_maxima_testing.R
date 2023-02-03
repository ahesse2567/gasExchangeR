library(gasExchangeR)
library(tidyverse)
library(janitor)
library(broom)
library(Deriv)
library(Ryacas)
library(readxl)
library(devtools)

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

df_avg <- avg_exercise_test(df_unavg, type = "breath", subtype = "rolling",
                            time_col = "time", roll_window = 15)

ggplot(data = df_avg, aes(x = vo2)) +
    geom_point(aes(y = ve_vo2), color = "purple") +
    geom_point(aes(y = ve_vco2), color = "green") +
    theme_minimal()

vt2_dat <- d1_poly_reg_maxima(.data = df_avg, .x = "vo2", .y = "ve_vco2", bp = "vt2")
vt2_dat

vent_eq_plot <- df_avg %>%
    arrange(vo2, time) %>%
    ggplot(aes(x = vo2)) +
    geom_point(aes(y = ve_vo2), color = "purple") +
    geom_point(aes(y = ve_vco2), color = "green") +
    theme_minimal() +
    geom_vline(xintercept = vt2_dat$breakpoint_data$vo2) +
    geom_line(aes(y = vt2_dat$lm_poly$fitted.values))
vent_eq_plot

vt1_dat <- df_avg %>%
    # filter(vo2 < vt2_dat$breakpoint_data$vo2) %>%
    d1_poly_reg_maxima(.x = "vo2", .y = "ve_vo2", bp = "vt1", degree = 5)
vt1_dat
