library(tidyverse)
library(gasExchangeR)

# load data
file_lines <- readLines("inst/extdata/Anton_vo2max.txt")
df_raw <- read.table(textConnection(file_lines[-2]), header = TRUE, sep="\t")

# basic cleaning, renaming variables
df_unavg <- df_raw %>%
    as_tibble() %>%
    janitor::clean_names() %>%
    separate(`time`, into = c("m1", "s1"), sep = ":") %>%
    separate(ex_time, into = c("m2", "s2"), sep = ":") %>%
    separate(time_clock,
             into = c("h3", "m3", "s3"),
             sep = ":") %>%
    mutate(across(where(is.character), as.numeric)) %>%
    mutate(time = (m1*60 + s1), .keep = "unused") %>%
    mutate(ex_time = (m2*60 + s2 ), .keep = "unused") %>%
    mutate(clock_time = hms::hms(s3, m3, h3), .keep = "unused") %>%
    relocate(contains("time")) %>%
    filter(!is.na(ex_time)) %>%
    filter(speed >= 4.5 & ex_time >= 750) %>%
    select(-c(time, hrr)) %>%
    rename(time = ex_time,
           vo2_kg = vo2,
           vo2 = vo2_1,
           ve = ve_btps) %>%
    mutate(time = time - min(time),
           ve_vo2 = ve / vo2 * 1000,
           ve_vco2 = ve/vco2*1000,
           excess_co2 = vco2^2 / vo2 - vco2)

# remove outliers
df_unavg_no_outliers <- ventilatory_outliers(df_unavg,
                                             outlier_cols = "vo2",
                                             time = "time",
                                             sd_lim = 3,
                                             width = 5,
                                             mos = "mean",
                                             align = "center",
                                             use_global_sd = TRUE,
                                             global_sd_mos = "median",
                                             exclude_test_val = TRUE,
                                             remove_outliers = TRUE,
                                             max_passes = Inf,
                                             plot_outliers = FALSE)

df_avg <- avg_exercise_test(df_unavg_no_outliers,
                            method = "time",
                            calc_type = "bin",
                            bin_w = 10)

ggplot(data = df_unavg_no_outliers, aes(x = time, y = vo2)) +
    geom_point(color = "red", alpha = 0.5) +
    geom_point(data = df_avg, aes(x = time, y = vo2)) +
    theme_bw()

# undebug(vo2_plateau.vo2max_neighbor)

vo2_plateau(df_avg,
            method = "vo2max_neighbor")


