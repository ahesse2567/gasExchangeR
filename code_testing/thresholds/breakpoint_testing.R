library(gasExchangeR)
library(devtools)
library(tidyverse)
library(janitor)
library(readxl)

df_unavg <- read_xlsx("../gasExchangeR_validation/data/processed/rand_15_cpet_exercisethresholds/mar22_117_post_gxt.xlsx") %>%
    rename(vo2 = vo2_abs,
           vo2_kg = vo2)

# file_lines <- readLines("inst/extdata/Anton_vo2max.txt")
# df_raw <- read.table(textConnection(file_lines[-2]), header = TRUE, sep="\t")
#
# df_unavg <- df_raw %>%
#     as_tibble() %>%
#     clean_names() %>%
#     separate(`time`, into = c("m1", "s1"), sep = ":") %>%
#     separate(ex_time, into = c("m2", "s2"), sep = ":") %>%
#     separate(time_clock,
#              into = c("h3", "m3", "s3"),
#              sep = ":") %>%
#     mutate(across(where(is.character), as.numeric)) %>%
#     mutate(time = (m1*60 + s1), .keep = "unused") %>%
#     mutate(ex_time = (m2*60 + s2 ), .keep = "unused") %>%
#     mutate(clock_time = hms::hms(s3, m3, h3), .keep = "unused") %>%
#     relocate(contains("time")) %>%
#     filter(!is.na(ex_time)) %>%
#     filter(speed >= 4.5 & ex_time >= 750) %>%
#     select(-time) %>%
#     rename(time = ex_time,
#            vo2_kg = vo2,
#            vo2 = vo2_1,
#            ve = ve_btps) %>%
#     mutate(ve_vo2 = ve / vo2 * 1000,
#            ve_vco2 = ve/vco2*1000,
#            excess_co2 = vco2^2 / vo2 - vco2) %>%
#     ventilatory_outliers(plot_outliers = FALSE)

df_avg <- avg_exercise_test(df_unavg, type = "time", subtype = "bin",
                            time_col = "time", bin_w = 15)

ggplot(data = df_avg, aes(x = time)) +
    geom_point(aes(y = vo2, color = "vo2"), alpha = 0.5) +
    geom_line(aes(y = vo2, color = "vo2"), alpha = 0.5) +
    geom_point(aes(y = vco2, color = "vco2"), alpha = 0.5) +
    geom_line(aes(y = vco2, color = "vco2"), alpha = 0.5) +
    scale_color_manual(name = NULL,
                       values = c("vo2" = "red", "vco2" = "blue")) +
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
    geom_point(aes(y = ve_vco2, color = "ve_vco2"), alpha = 0.5) +
    geom_line(aes(y = ve_vco2, color = "ve_vco2"), alpha = 0.5) +
    geom_point(aes(y = ve_vo2, color = "ve_vo2"), alpha = 0.5) +
    geom_line(aes(y = ve_vo2, color = "ve_vo2"), alpha = 0.5) +
    ggtitle("Ventilatory Equivalents") +
    scale_color_manual(name = NULL,
                       values = c("ve_vco2" = "purple", "ve_vo2" = "green")) +
    theme_minimal()

debug(gasExchangeR::orr)
# debug(d2_reg_spline_maxima)
# Rprof()
breakpoint(df_avg,
           algorithm_vt1 = "v-slope",
           x_vt1 = "vo2",
           y_vt1 = "vco2",
           algorithm_vt2 = "d2_inflection",
           x_vt2 = "vco2",
           y_vt2 = "ve",
           bp = "both",
           truncate = TRUE,
           pos_slope_after_bp = TRUE,
)
# Rprof(NULL)
prof_summary <- summaryRprof()

prof_summary$by.self

breakpoint(df_avg,
           algorithm_vt1 = "orr",
           x_vt1 = "vo2",
           y_vt1 = "vco2",
           algorithm_vt2 = "orr",
           x_vt2 = "vco2",
           y_vt2 = "ve",
           bp = "both",
           truncate = TRUE,
           pos_slope_after_bp = TRUE,
)

undebug(breakpoint)
undebug(d2_reg_spline_maxima)


bp_dat <- breakpoint(.data = df_avg, method = "excess_co2",
                     x_vt1 = "vo2",
                     algorithm_vt1 = "d2_poly_reg_maxima",
                     algorithm_vt2 = "d2_poly_reg_maxima",
                     x_vt2 = "ve", y_vt2 = "ve_vco2", truncate = FALSE,
                     front_trim_vt1 = 90)
bp_dat$vt1_dat$bp_plot
bp_dat$vt2_dat$bp_plot

debug(add_threshold_lines)
undebug(add_threshold_lines)
debug(breakpoint)

bp_dat <- breakpoint(.data = df_avg,
                     algorithm_vt1 = "d2_reg_spline_maxima",
                     x_vt1 = "vo2", y_vt1 = "vco2",
                     algorithm_vt2 = "d2_reg_spline_maxima",
                     x_vt2 = "vo2", y_vt2 = "vco2",
                     bp = "both", truncate = FALSE)
bp_dat$bp_dat
bp_dat$bp_plots$vt1_plot
bp_dat$bp_plots$vt2_plot
bp_dat$vt1_dat$bp_plot
bp_dat$vt2_dat$bp_plot


nine_panel_plot(.data = df_avg,
                vt1_dat = bp_dat$vt1_dat,
                vt2_dat = bp_dat$vt2_dat,
                vo2 = "vo2",
                vco2 = "vco2",
                ve = "ve",
                hr = "hr",
                time = "time",
                vt = "vt_btps")



# debug(d2_reg_spline_maxima)
bp_dat <- breakpoint(.data = df_avg, x_vt1 = "vo2", y_vt1 = "vco2",
                     algorithm_vt1 = "d2_reg_spline_maxima",
                     algorithm_vt2 = "d2_reg_spline_maxima",
                     x_vt2 = "vco2", y_vt2 = "ve",
                     pos_change_vt2 = TRUE, truncate = TRUE)
bp_dat$vt1_dat$bp_plot
bp_dat$vt2_dat$bp_plot

bp_dat <- breakpoint(.data = df_avg, x_vt1 = "vo2", y_vt1 = "ve",
                     algorithm_vt1 = "d2_inflection",
                     algorithm_vt2 = "d2_inflection",
                     x_vt2 = "vco2", y_vt2 = "ve", truncate = TRUE)
bp_dat$vt1_dat$bp_plot
bp_dat$vt2_dat$bp_plot

bp_dat <- breakpoint(.data = df_avg, x_vt1 = "vo2", y_vt1 = "vco2",
                     algorithm_vt1 = "dmax",
                     algorithm_vt2 = "dmax",
                     x_vt2 = "vco2", y_vt2 = "ve", truncate = TRUE,
                     ordering = "time")
bp_dat$vt1_dat$bp_plot
bp_dat$vt2_dat$bp_plot
# debug(check_if_determinant_bp)
bp_dat <- breakpoint(.data = df_avg, x_vt1 = "vo2", y_vt1 = "vco2",
                     algorithm_vt1 = "v-slope",
                     algorithm_vt2 = "spline_bp",
                     x_vt2 = "vco2", y_vt2 = "ve", truncate = TRUE,
                     ordering = "time")
bp_dat$bp_dat
bp_dat$vt1_dat$bp_plot
bp_dat$vt2_dat$bp_plot

bp_dat <- breakpoint(.data = df_avg, x_vt1 = "vo2", y_vt1 = "vco2",
                     algorithm_vt1 = "v-slope",
                     algorithm_vt2 = "d2_reg_spline_maxima",
                     x_vt2 = "vo2", y_vt2 = "ve_vco2", truncate = TRUE,
                     ordering = "time")
bp_dat$vt1_dat$bp_plot
bp_dat$vt2_dat$bp_plot
debug(orr)
# undebug(orr)
bp_dat <- breakpoint(.data = df_avg, x_vt1 = "vo2", y_vt1 = "vco2",
                     algorithm_vt1 = "orr",
                     algorithm_vt2 = "orr",
                     x_vt2 = "vco2", y_vt2 = "ve", truncate = TRUE,
                     ordering = "time")
bp_dat$vt1_dat$bp_plot
bp_dat$vt2_dat$bp_plot

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

