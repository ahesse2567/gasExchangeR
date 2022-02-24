library(devtools)
library(gasExchangeR)
library(tidyverse)

df_raw <- read_csv("inst/extdata/mar16_101_pre.csv")
df <- df_raw %>%
    rename_with(.fn = tolower) %>%
    rename(vo2_rel = vo2,
           vo2_abs = vo2.1,
           ve = ve.btps) %>%
    trim_rest_rec(intensity_col = "speed") %>%
    select(time, speed, grade, vo2_rel, vo2_abs, vco2, ve, peto2, petco2)

df_avg <- avg_exercise_test(df,
                            type = "breath",
                            subtype = "bin_roll",
                            roll_window = 20,
                            bin_w = 4)
df_avg
ggplot(data = df_avg, aes(x = time, y = vo2_rel)) +
    geom_point()
