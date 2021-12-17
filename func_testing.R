library(tidyverse)

df <- read_csv("inst/extdata/mar16_101_pre.csv")

df <- df %>%
    trim_rest_rec(intensity_col = "speed", start_intensity = 3) %>%
    # vo2max_window() %>%
    # x_breath_mean(b = 4) %>%
    select(time, speed, grade, vo2, vo2.1, vco2, ve.btps, peto2, petco2) %>%
    rename(vo2_rel = vo2,
           vo2_abs = vo2.1,
           ve = ve.btps)

ggplot(data = df, aes(x = time, y = vo2_rel)) +
    geom_point() +
    theme_bw()

avg_exercise_test(df, type = "breath", subtype = "rolling", mos = "mean")
avg_exercise_test(df, type = "breath", subtype = "rolling", mos = "mean", trim = 4)
avg_exercise_test(df, type = "breath", subtype = "rolling", mos = "median", trim = 0)

avg_exercise_test(df, type = "breath", subtype = "bin", mos = "mean")


#################################
## Digital filter option

butter_lowpass <- function(cutoff, fs, order = 5){
    nyq <- 0.5 * fs # nyquist frequency is half the sampling rate (fs) b/c you need
    # at a minimum two data points per wave in order to construct the wave
    normal_cutoff <- cutoff / nyq
    bf <- signal::butter(n = order, W = normal_cutoff, type = "low", plane = "z")
    bf
}

butter_lowpass_filter <- function(data, cutoff, fs, order = 5){
    # browser()
    bf <- butter_lowpass(cutoff, fs, order=order)
    filtered_data <- signal::filter(bf, data)
    filtered_data
}

data_num <- df %>% # coerce to numeric b/c time may not be of another class
    dplyr::mutate(dplyr::across(where(purrr::negate(is.character)),
                                as.numeric))

fs <- 1 # sampling rate in Hz
cutoff <- 0.04 # low-pass filter cutoff (sampling rate?)
order <- 3

bf <- butter_lowpass(cutoff = cutoff, fs = fs, order = order)
bf
head(data_num$vo2_rel)
head(signal::filter(bf, data_num$vo2_rel))

out <- purrr::map(.x = data_num,
           .f = function(.x, bf) signal::filter(bf, .x), bf = bf)
head(out$vo2_rel)
head(signal::filter(bf, data_num$vo2_rel))

out$time

out$vo2_rel

signal::filter(bf, data_num$vo2_rel)

df

df_digi_filt <- avg_exercise_test(df, type = "digital")

head(df_digi_filt$vo2_rel)
ggplot(data = df, aes(x = time, y = vo2_rel)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    geom_line(aes(y = df_digi_filt$vo2_rel), color = "red")

