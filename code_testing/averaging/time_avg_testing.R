library(devtools)
library(gasExchangeR)
library(tidyverse)

df_raw <- read_csv("inst/extdata/mar16_101_pre.csv")
df_unavg <- df_raw %>%
    rename_with(.fn = tolower) %>%
    rename(vo2_rel = vo2,
           vo2_abs = vo2.1,
           ve = ve.btps) %>%
    trim_pre_post(intensity_col = "speed",
                  pre_ex_intensity = 5.9) %>%
    select(time, speed, grade, vo2_rel, vo2_abs, vco2, ve, peto2, petco2)

ggplot(data = df_unavg, aes(x = time, y = vo2_abs)) +
    geom_point()

df_avg <- avg_exercise_test(df_unavg,
                            type = "time",
                            subtype = "bin",
                            bin_w = 15)
df_avg
ggplot(data = df_avg, aes(x = time, y = vo2_rel)) +
    geom_point()
