library(gasExchangeR)
library(devtools)
library(tidyverse)
library(devtools)
library(janitor)

file_lines <- readLines("inst/extdata/Anton_vo2max.txt")
df_raw <- read.table(textConnection(file_lines[-2]), header = TRUE, sep="\t")

df_unavg <- df_raw %>%
    as_tibble() %>%
    clean_names() %>%
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

df_avg <- avg_exercise_test(df_unavg, type = "breath", subtype = "rolling",
                            time_col = "time", roll_window = 15)

ggplot(data = df_avg, aes(x = time)) +
    geom_point(aes(y = vo2), color = "red", alpha = 0.5) +
    geom_point(aes(y = vco2), color = "blue", alpha = 0.5) +
    theme_minimal()

ggplot(data = df_avg, aes(x = vo2, y = vco2)) +
    geom_point(color = "blue", alpha = 0.5) +
    ggtitle("V-slope") +
    theme_minimal()

ggplot(data = df_avg, aes(x = vco2, y = ve)) +
    geom_point(color = "orange", alpha = 0.5) +
    ggtitle("Respiratory Compensation") +
    theme_minimal()

ggplot(data = df_avg, aes(x = time)) +
    geom_point(aes(y = ve_vco2), color = "green", alpha = 0.5) +
    geom_line(aes(y = ve_vco2), color = "green", alpha = 0.5) +
    geom_point(aes(y = ve_vo2), color = "purple", alpha = 0.5) +
    geom_line(aes(y = ve_vo2), color = "purple", alpha = 0.5) +
    ggtitle("Ventilatory Equivalents") +
    theme_minimal()

vt1_only_algs <- c("v-slope", "d1_crossing", )

algorithms <- c("jm", "orr", "dmax", "spline_bp",
                "d2_inflection", "d1_poly_reg_maxima",
                "d2_reg_spline_maxima")

bp_dat <- breakpoint(.data = df_avg, method = "excess_co2",
                     algorithm_vt1 = "spline_bp",
                     algorithm_vt2 = "d2_poly_reg_maxima",
                     x_vt2 = "vo2", y_vt2 = "ve_vco2", truncate = TRUE,
                     front_trim_vt1 = 90)

bp_dat <- breakpoint(.data = df_avg, method = "excess_co2",
                     algorithm_vt1 = "spline_bp",
                     algorithm_vt2 = "d2_poly_reg_maxima",
                     x_vt2 = "vo2", y_vt2 = "ve_vco2", truncate = TRUE,
                     front_trim_vt1 = 90)

bp_dat <- breakpoint(.data = df_avg, x_vt1 = "time", y_vt1 = "vo2",
                     algorithm_vt1 = "d1_crossing",
                     algorithm_vt2 = "d2_poly_reg_maxima",
                     x_vt2 = "vo2", y_vt2 = "ve_vco2", truncate = FALSE,
                     front_trim_vt1 = 0)

bp_dat <- breakpoint(.data = df_avg, x_vt1 = "vo2", y_vt1 = "ve",
                     algorithm_vt1 = "d2_inflection",
                     algorithm_vt2 = "d2_inflection",
                     x_vt2 = "vco2", y_vt2 = "ve", truncate = TRUE,
                     front_trim_vt1 = 0)

bp_dat <- breakpoint(.data = df_avg, x_vt1 = "vo2", y_vt1 = "vco2",
                     algorithm_vt1 = "dmax",
                     algorithm_vt2 = "dmax",
                     x_vt2 = "vco2", y_vt2 = "ve", truncate = TRUE)

bp_dat <- breakpoint(.data = df_avg, x_vt1 = "vo2", y_vt1 = "vco2",
                     algorithm_vt1 = "jm",
                     algorithm_vt2 = "jm",
                     x_vt2 = "vco2", y_vt2 = "ve", truncate = TRUE)

bp_dat <- breakpoint(.data = df_avg, x_vt1 = "vo2", y_vt1 = "vco2",
                     algorithm_vt1 = "v-slope",
                     algorithm_vt2 = "d2_reg_spline_maxima",
                     x_vt2 = "vo2", y_vt2 = "ve_vco2", truncate = TRUE)

bp_dat <- breakpoint(.data = df_avg, x_vt1 = "vo2", y_vt1 = "vco2",
                     algorithm_vt1 = "orr",
                     algorithm_vt2 = "orr",
                     x_vt2 = "vco2", y_vt2 = "ve", truncate = TRUE)

print(bp_dat$bp_dat, width = Inf)

bp_dat$vt1_dat$bp_plot
bp_dat$vt2_dat$bp_plot

ggplot(data = df_avg, aes(x = vo2)) +
    geom_point(aes(y = ve_vco2), color = "green", alpha = 0.5) +
    # geom_line(aes(y = ve_vco2), color = "green", alpha = 0.5) +
    geom_point(aes(y = ve_vo2), color = "purple", alpha = 0.5) +
    # geom_line(aes(y = ve_vo2), color = "purple", alpha = 0.5) +
    ggtitle("Ventilatory Equivalents") +
    geom_vline(xintercept = bp_dat$bp_dat$vo2) +
    theme_minimal()

ggplot(data = df_avg, aes(x = vo2)) +
    geom_point(aes(y = excess_co2), alpha = 0.5) +
    ggtitle("Excess CO2") +
    geom_vline(xintercept = bp_dat$bp_dat$vo2) +
    theme_minimal()

ggplot(data = df_avg, aes(x = vo2, y = vco2)) +
    geom_point(color = "blue", alpha = 0.5) +
    geom_vline(xintercept = bp_dat$bp_dat$vo2) +
    ggtitle("V-slope") +
    theme_minimal()

