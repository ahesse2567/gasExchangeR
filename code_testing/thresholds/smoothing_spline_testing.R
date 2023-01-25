library(gasExchangeR)
library(tidyverse)
library(janitor)
library(broom)
library(splines)
library(zoo)

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
    mutate(ve_vo2 = ve*1000 / vo2_abs,
           ve_vco2 = ve*1000 / vco2,
           excess_co2 = (vco2^2 / vo2_abs) - vco2) %>%
    filter(!is.na(time)) %>%
    filter(stage >=6) %>%
    relocate(time, speed, grade)
df_unavg

# plot raw VO2 vs. time data
ggplot(data = df_unavg, aes(x = time, y = vo2_rel)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_line(alpha = 0.5) +
    theme_minimal() +
    ggtitle("Raw Data")

df_avg <- avg_exercise_test(df_unavg, type = "breath", subtype = "rolling",
                            time_col = "time", roll_window = 9, roll_trim = 4)

ggplot(data = df_avg, aes(x = time, y = vo2_rel)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_line(alpha = 0.5) +
    theme_minimal() +
    ggtitle("Averaged Data")

.data <- df_avg
.x <- "time"
.y <- "ve_vco2"
time <- "time"
degree <- 3
alpha_linearity = 0.05

smoothing_spline <- smooth.spline(x = .data[[.x]], y = .data[[.y]], cv = TRUE)
predict(smoothing_spline)

zoo::rollmean(x = predict(smoothing_spline, x = .data[[.x]], deriv = 2)$y,
              k = 10, fill = NA)


ggplot(data = .data, aes(x = .data[[.x]], y = .data[[.y]])) +
    geom_line(color = "green") +
    geom_line(aes(y = smoothing_spline$y)) +
    geom_line(aes(y = predict(smoothing_spline, x = .data[[.x]], deriv = 1)$y),
              linetype = "dashed") +
    geom_line(aes(y = zoo::rollmean(x = predict(smoothing_spline, x = .data[[.x]], deriv = 2)$y,
                                    k = 10, fill = NA)),
              linetype = "dotted") +
    theme_minimal()

