library(gasExchangeR)
library(tidyverse)
library(devtools)

df_raw <- read_csv("inst/extdata/37818_vo2max_unavg.csv",
                   show_col_types = FALSE)

df_unavg <- df_raw %>%
    rename_all(.funs = tolower) %>%
    rename(vo2_rel = `vo2...4`,
           vo2_abs = `vo2...5`,
           ve = `ve btps`,
           vt = `vt btps`) %>%
    select(time, speed, grade, vo2_rel, vo2_abs, vco2, ve, vt) %>%
    mutate(ve_vo2 = ve*1000 / vo2_abs,
           ve_vco2 = ve*1000 / vco2) %>%
    trim_pre_post(intensity_col = "grade",
                  pre_ex_intensity = 300.1,
                  post_ex_intensity = 300.1)

df_unavg <- df_unavg %>%
    mutate(progress = normalize(df_unavg$vo2_abs) +
               normalize(df_unavg$vco2) +
               normalize(df_unavg$ve))

ggplot(data = df_unavg, aes(x = time)) +
    geom_point(aes(y = vo2_abs), color = "red", alpha = 0.5) +
    geom_point(aes(y = vco2), color = "blue", alpha = 0.5) +
    theme_bw()

ggplot(data = df_unavg, aes(x = time)) +
    geom_point(aes(y = ve), color = "orange", alpha = 0.5) +
    theme_bw()

ggplot(data = df_unavg, aes(x = time)) +
    geom_point(aes(y = progress)) +
    theme_bw()

