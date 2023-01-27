library(gasExchangeR)
library(tidyverse)
library(janitor)
library(broom)
library(Deriv)
library(Ryacas)
library(readxl)
library(devtools)

# for derivative methods, do users need to choose how they want the function
# created, and then what derivative method to use? e.g. they say they want
# to use a polynomial regression, regression splines, or smoothing splines.
# then, they say if they want to use 2nd derivative maxima, 1st derivative crossing
# or 2nd derivative inflection

# d1 crossing is most similar to Wisén and Wohlfart (2004)
# this is a VT1 method ONLY

df_colnames <- read_xlsx("inst/extdata/Foreman_Ramp_10_18_2022.xlsx",
                    n_max = 0) %>%
    colnames()
df_raw <- read_xlsx("inst/extdata/Foreman_Ramp_10_18_2022.xlsx",
                    col_names = FALSE, skip = 3) %>%
    set_names(df_colnames)

df_unavg <- df_raw %>%
    clean_names() %>%
    rename(time = t) %>%
    filter(phase == "EXERCISE") %>%
    relocate(time, speed, grade)

df_avg <- avg_exercise_test(df_unavg, type = "breath", subtype = "rolling",
                            time_col = "time", roll_window = 9, roll_trim = 4)

# plot VO2 vs. time data
ggplot(data = df_avg, aes(x = time, y = vo2)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_line(alpha = 0.5) +
    geom_point(aes(y = vco2), alpha = 0.5, color = "blue") +
    geom_line(aes(y = vco2), alpha = 0.5) +
    theme_minimal()

# find 1st derivative of vo2 vs. time
# find 1st derivative of vco2 vs. time
# find where they cross

d1_crossing <- function(.data,
                        .x,
                        .y,
                        degree = NULL,
                        vo2 = "vo2",
                        vco2 = "vco2",
                        ve = "ve",
                        time = "time",
                        alpha_linearity = 0.05, # change to just alpha?
                        bp) {

