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

df_avg <- avg_exercise_test(df_unavg, type = "breath", subtype = "rolling",
                  time_col = "time", roll_window = 9, roll_trim = 4)

ggplot(df_avg, aes(x = time, y = petco2)) +
    geom_point()


orr(.data = df_avg, .x = "vco2", .y = "ve",
    vo2 = "vo2_abs", vco2 = "vco2", ve = "ve",
    time = "time", alpha_linearity = 0.05, bp = "vt1")

breakpoint(.data = df_avg,
           x_vt1 = "vo2_abs",
           y_vt1 = "vco2",
           algorithm_vt1 = "orr",
           x_vt2 = "vco2",
           y_vt2 = "ve",
           algorithm_vt2 = "orr",
           vo2 = "vo2_abs",
           bps = "both")
