library(gasExchangeR)
library(tidyverse)
library(devtools)
library(janitor)
library(broom)
library(Deriv)
library(Ryacas)

# big thanks to this stack overflow thread: https://stackoverflow.com/questions/40438195/function-for-polynomials-of-arbitrary-order-symbolic-method-preferred
# splineDesign function to get higher order spline functions?

# df_raw <- read_delim("inst/extdata/Aaron_VO2max.txt",
#                      delim = "\t",
#                      na = c("", "NA", " "),
#                      show_col_types = FALSE)

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

df_avg <- avg_exercise_test(df_unavg, type = "time", subtype = "bin", bin_w = 10)

# plot raw VO2 and VCO2 vs. time data
ggplot(data = df_unavg, aes(x = time, y = vo2_abs)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_line(alpha = 0.5, color = "red") +
    geom_point(aes(y = vco2), alpha = 0.5, color = "blue") +
    geom_line(aes(y = vco2), alpha = 0.5, color = "blue") +
    theme_minimal() +
    ggtitle("Raw Data")

# plot averaged VO2 and VCO2 vs. time data
ggplot(data = df_avg, aes(x = time, y = vo2_abs)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_line(alpha = 0.5, color = "red") +
    geom_point(aes(y = vco2), alpha = 0.5, color = "blue") +
    geom_line(aes(y = vco2), alpha = 0.5, color = "blue") +
    theme_minimal() +
    ggtitle("Averaged Data")

ggplot(data = df_avg, aes(x = vo2_abs, y = vco2)) +
    geom_point(color = "blue", alpha = 0.5) +
    theme_minimal()

ggplot(data = df_avg, aes(x = time)) +
    geom_point(aes(y = ve_vo2), color = "purple", alpha = 0.5) +
    geom_line(aes(y = ve_vo2), color = "purple", alpha = 0.5) +
    geom_point(aes(y = ve_vco2), color = "green", alpha = 0.5) +
    geom_line(aes(y = ve_vco2), color = "green", alpha = 0.5) +
    theme_minimal() +
    ggtitle("Ventilatory Equivalents")

ggplot(data = df_avg, aes(x = time, y = excess_co2)) +
    geom_point(alpha = 0.5) +
    geom_line() +
    theme_minimal() +
    ggtitle("Excess CO2")


.x <- "vco2"
.y <- "ve"
time <- "time"
.data <- df_avg
degree <- NULL
alpha_linearity = 0.05
bp = "vt2"
degree <- 1

spline_bp_dat <- spline_bp(.data = df_avg, .x = "vco2", .y = "ve",
          bp = "vt2", vo2 = "vo2_abs")
spline_bp_dat$breakpoint_data
spline_bp_dat$bp_plot
