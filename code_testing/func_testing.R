library(gasExchangeR)
library(tidyverse)

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

ggplot(data = df_unavg, aes(x = time, y = vo2_abs)) +
    geom_point(color = "red", alpha = 0.5) +
    theme_bw()

df_avg <- avg_exercise_test(df_unavg, type = "breath", subtype = "rolling",
                  time_col = "time", roll_window = 9, trim = 4)

ggplot(data = df_avg, aes(x = time, y = vo2_abs)) +
    geom_point(color = "red", alpha = 0.5) +
    theme_bw()

breakpoint(df_avg, algorithm_aert = "jm", x_aert = "ve_vo2", y_aert = "time",
           algorithm_rc = "jm", x_rc = "ve_vco2", y_rc = "time", vo2 = "vo2_abs")
