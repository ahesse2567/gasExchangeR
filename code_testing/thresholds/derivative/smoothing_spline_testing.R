library(gasExchangeR)
library(tidyverse)
library(janitor)
library(broom)
library(splines)
library(zoo)

# Deprase the function and see where they set the order
# probably go to their GitHub
# show_methods()
# Maybe just go into the package folder on my computer
# .libpaths(function)
# Probably okay to add 5th order smoothing spline later

# I think I need to use a 5th deriv smoothing spline a la Sherrill et al. (1990)
# I'm starting to get skeptical that Sherrill et al. (1990) actually used true
# "smoothing splines". Sherrill et al. (1990) specifically mentioned a 5th order
# polynomial spline
# Also, realize that Sherrill's study specifically used VCO2 vs. VO2

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

ggplot(data = df_avg, aes(x = time, y = vo2_kg)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_line(alpha = 0.5) +
    theme_minimal() +
    ggtitle("Averaged Data")

.data <- df_avg
.x <- "vo2"
.y <- "vco2"
time <- "time"
degree <- 3
alpha_linearity = 0.05

smoothing_spline <- smooth.spline(x = .data[[.x]], y = .data[[.y]], cv = TRUE)
# predict(smoothing_spline)

.data %>%
    arrange(.data[[.x]], .data[[time]]) %>%
    ggplot(aes(x = .data[[.x]], y = .data[[.y]])) +
    geom_point(color = "blue", alpha = 0.5) +
    geom_line(aes(y = smoothing_spline$y)) +
    # geom_line(aes(y = predict(smoothing_spline, x = .data[[.x]], deriv = 1)$y),
    #           linetype = "dashed") +
    # geom_line(aes(y = zoo::rollmean(x = predict(smoothing_spline, x = .data[[.x]], deriv = 2)$y,
    #                                 k = 10, fill = NA)),
    #           linetype = "dotted") +
    theme_minimal()

ggplot(data = .data, aes(x = .data[[.x]], y = .data[[.y]])) +
    geom_line(aes(y = predict(smoothing_spline, x = .data[[.x]], deriv = 1)$y),
              linetype = "dotted") +
    # geom_line(aes(y = zoo::rollmean(x = predict(smoothing_spline,
    #                                             x = .data[[.x]],
    #                                             deriv = 2)$y,
    #                                 k = 10, fill = NA)), linetype = "dotted",
    #           color = "red") +
    # ggtitle("2nd derivative smoothing spline") +
    theme_minimal()

ggplot(data = .data, aes(x = .data[[.x]], y = .data[[.y]])) +
    geom_line(aes(y = zoo::rollmean(x = predict(smoothing_spline,
                                                x = .data[[.x]],
                                                deriv = 2)$y,
                                    k = 10, fill = NA)), linetype = "dotted",
              color = "red") +
    ggtitle("2nd derivative smoothing spline - 10-point rolling average") +
    theme_minimal()

# compare to The Respiratory Compensation Point is Not a Valid Surrogate for Critical Power. That paper has a MUCH smoother curve.



