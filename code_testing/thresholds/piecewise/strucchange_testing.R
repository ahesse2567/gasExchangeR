library(strucchange)
library(gasExchangeR)
library(tidyverse)
library(devtools)
library(lubridate)

df_raw <- read_csv("inst/extdata/37818_vo2max_unavg.csv",
                   show_col_types = FALSE)
df_unavg <- df_raw %>%
    janitor::clean_names() %>%
    rename(vo2_rel = vo2_4,
           vo2_abs = vo2_5,
           ve = `ve_btps`,
           vt = `vt_btps`) %>%
    select(time, speed, grade, vo2_rel, vo2_abs, vco2, ve) %>%
    trim_pre_post(intensity_col = "grade",
                  pre_ex_intensity = 300.1,
                  post_ex_intensity = 300.1) %>%
    filter(speed >= 3.6) %>%
    mutate(ve_vo2 = ve / (vo2_abs/1000),
           ve_vco2 = ve / (vco2/1000),
           time = as.numeric(ms(str_remove(as.character(time), ":00")))) %>%
    ventilatory_outliers(outlier_cols = "vo2_abs", plot_outliers = FALSE)

df_avg <- avg_exercise_test(df_unavg,
                            method = "time",
                            calc_type = "bin",
                            bin_w = 10,
                            time_col = "time")

fs <- Fstats(df_avg$vo2_abs ~ df_avg$ve_vo2, from = 0.1)
plot(fs)
bps <- breakpoints(fs)
bps$breakpoints
plot(df_avg$vo2_abs, df_avg$ve_vo2)
abline(v = df_avg$vco2[bps$breakpoints])

data("Nile")
plot(Nile)

fs.nile <- Fstats(Nile ~ 1)
plot(fs.nile)
sctest(fs.nile)
lines(breakpoints(fs.nile))
