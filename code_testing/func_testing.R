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

ggplot(data = df_unavg, aes(x = vco2, y = ve)) +
    geom_point(color = "orange", alpha = 0.5) +
    theme_bw()

ggplot(data = df_unavg, aes(x = time, y = vo2_abs)) +
    geom_point(color = "red", alpha = 0.5) +
    theme_bw()

df_avg <- avg_exercise_test(df_unavg, type = "breath", subtype = "rolling",
                  time_col = "time", roll_window = 9, roll_trim = 4)

ggplot(data = df_avg, aes(x = time, y = ve)) +
    geom_point(color = "orange", alpha = 0.5) +
    theme_bw()

ggplot(data = df_avg, aes(x = vco2, y = ve)) +
    geom_point(color = "brown", alpha = 0.5) +
    theme_bw()

rc_stuff <- breakpoint(df_avg, algorithm_vt1 = "jm", x_vt1 = "vo2", y_vt1 = "vco2",
           algorithm_rc = "jm", x_rc = "vo2_abs", y_rc = "vco2", vo2 = "vo2_abs")

rc_stuff

ggplot(data = df_avg, aes(x = vo2_abs, y = vco2)) +
    geom_point(color = "blue", alpha = 0.5) +
    theme_bw() +
    geom_point(data = rc_stuff[[2]], aes(x = x, y = y_hat)) +
    geom_vline(xintercept = rc_stuff[[1]][["vo2_abs"]])

# ggplot(data = df_avg, aes(x = vco2, y = ve)) +
#     geom_point(color = "orange", alpha = 0.5) +
#     theme_bw() +
#     geom_point(data = rc_stuff[[2]], aes(x = x, y = y_hat))


