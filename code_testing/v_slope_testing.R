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

.data = df_avg
x_vt1 = "vo2_abs"
y_vt1 = "vco2"
algorithm_vt1 = "v-slope"
x_vt2 = "vco2"
y_vt2 = "ve"
algorithm_vt2 = "v-slope"
vo2 = "vo2_abs"
bps = "both"
.x <- "vco2"
.y <- "ve"
slope_change_lim <- 0.1

breakpoint(.data = df_avg,
           x_vt1 = "vo2_abs",
           y_vt1 = "vco2",
           algorithm_vt1 = "v-slope",
           x_vt2 = "vco2",
           y_vt2 = "ve",
           algorithm_vt2 = "v-slope",
           vo2 = "vo2_abs",
           bps = "both")

num_to_na <- function(x) {
    if(is.numeric(x)) {
        x <- NA
    }
    x
}

num_to_na(bp_dat$bp)

map_df(bp_dat, num_to_na)

na_if(is.numeric(bp_dat$ve), NA)

bp_dat[is.character(bp_dat)]

ss <- loop_v_slope(df_avg, "vo2_abs", "vco2")
plot(ss)
bp_idx <- v_slope(df_avg, "vo2_abs", "vco2")

df_left <- df_avg[1:bp_idx,]
df_right <- df_avg[(bp_idx+1):nrow(df_avg),]
lm_left <- lm(vco2 ~ 1 + vo2_abs, data = df_left)
lm_right <- lm(vco2 ~ 1 + vo2_abs, data = df_right)

plot_data_left <- tibble(vo2_abs = df_left$vo2_abs,
                         vco2 = predict(lm_left))
plot_data_right <- tibble(vo2_abs = df_right$vo2_abs,
                         vco2 = predict(lm_right))

ggplot(data = df_avg, aes(x = vo2_abs, y = vco2)) +
    geom_point(color = "purple", alpha = 0.5) +
    # geom_vline(xintercept = df_avg$vo2_abs[bp_idx]) +
    geom_line(data = plot_data_left, aes(x = vo2_abs, y = vco2)) +
    geom_line(data = plot_data_right, aes(x = vo2_abs, y = vco2)) +
    theme_bw()
    # geom_smooth(method = "lm", se = FALSE)

