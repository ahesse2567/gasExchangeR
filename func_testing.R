library(gasExchangeR)
library(tidyverse)

df <- read_csv("inst/extdata/mar16_101_pre.csv")

df <- df %>%
    trim_rest_rec(intensity_col = "speed") %>%
    # vo2max_window() %>%
    # x_breath_mean(b = 4) %>%
    select(time, speed, grade, vo2, vo2.1, vco2, ve.btps, peto2, petco2) %>%
    rename(vo2_rel = vo2,
           vo2_abs = vo2.1,
           ve = ve.btps) %>%
    mutate(ve_vo2 = ve*1000 / vo2_abs,
           ve_vco2 = ve*1000 / vco2)

ggplot() +
    geom_line(data = df[(nrow(df)-60):nrow(df),],
              aes(x = as.numeric(time), y = vo2_abs),
              color = "red",
              alpha = 0.5) +
    theme_bw() +
    ggtitle("Last 60 seconds, Unaveraged, Raw") +
    xlim(c(790, 925)) +
    ylim(c(2100, 4500))

df_knn <- exercise_outliers(df, method = "knn", keep_outliers = FALSE)

ggplot() +
    geom_line(data = df_knn[(nrow(df_knn)-60):nrow(df_knn),],
               aes(x = as.numeric(time), y = vo2_abs),
               color = "red",
               alpha = 0.5) +
    theme_bw() +
    ggtitle("Last 60 seconds, Unaveraged, Outliers Removed") +
    xlim(c(790, 925)) +
    ylim(c(2100, 4500))

df_avg <- avg_exercise_test(df_knn, type = "breath", subtype = "rolling")

ggplot() +
    geom_line(data = df_avg[(nrow(df_avg)-60):nrow(df_avg),],
               aes(x = as.numeric(time), y = vo2_abs),
               color = "red",
               alpha = 0.5) +
    theme_bw() +
    ggtitle("Last 60 seconds, Averaged") +
    xlim(c(790, 925)) +
    ylim(c(2100, 4500))

#### Plateau
df_avg
last_x_s <- 30
time_col <- "time"
vo2_col <- "vo2_abs"
delta_vo2 = 50
alpha <- 0.05

eot_idx <- which((df_avg[[time_col]] - max(df_avg[[time_col]]) + last_x_s) >= 0)
eot <- df_avg[eot_idx,]

ggplot(data = eot, aes(x = time, y = vo2_abs)) +
    geom_point(color = "red") +
    geom_line() +
    geom_smooth(method = "lm", se = FALSE) +
    theme_bw()

form <- as.formula(paste(vo2_col, "~ 1 +", time_col))
lm <- lm(form, data = eot)
broom::tidy(lm)

if(broom::tidy(lm)[["p.value"]][2] > alpha) {
    plateau <- TRUE
} else if(broom::tidy(lm)[["p.value"]][2] < alpha & coef(lm)[2]*60 < 0){
    plateau <- TRUE
} else {
    plateau <- FALSE
}
plateau

out <- tibble(plateau = plateau,
              vo2_time_slope_min = coef(lm)[2]*60,
              p.value = broom::tidy(lm)[["p.value"]][2])
out
vo2_plateau(df_avg, method = "zero_slope", vo2_col = "vo2_abs",
            time_col = "time", last_x_s = 35)


vo2max_idx <- which.max(df_avg[[vo2_col]])
vo2max_neighbor_diff <- df_avg[[vo2_col]][vo2max_idx] -
    df_avg[[vo2_col]][vo2max_idx - 1]

if(vo2max_neighbor_diff < delta_vo2) {
    plateau <- TRUE
} else {
    plateau <- FALSE
}

out <- tibble(plateau = plateau,
              vo2max_neighbor_diff = vo2max_neighbor_diff)
out


max(df_avg$vo2_abs) / max(df_avg$grade)

ggplot(data = df_avg, aes(x = time)) +
    geom_point(aes(y = vo2_abs, color = plateau), alpha = 0.5) +
    geom_line(aes(y = speed*(max(vo2_abs) / max(speed))), color = "green") +
    geom_line(aes(y = grade*(max(vo2_abs) / max(grade))), color = "blue") +
    scale_y_continuous(
        name = "Absolute VO2",
        sec.axis = sec_axis(trans = ~./(max(df_avg$vo2_abs) / max(df_avg$grade)),
                            name = "Speed/Grade")
        )+
    theme_bw() +
    ggtitle("Only using the difference from data point to data point says\nthere's a plateau for almost the entire test when the diff is >= 50")
# I think the issue is that this is a rolling average, whereas the traditional
# vo2 plateau criteria uses time-bin averagin


###### Last 60 s
last_min_idx <- which((df_avg$time - max(df_avg$time) + 60) >= 0)
last_min <- df_avg[last_min_idx,]

ggplot(data = last_min, aes(x = time)) +
    geom_point(aes(y = vo2_abs, color = plateau), alpha = 0.5) +
    # geom_line(aes(y = speed*(max(vo2_abs) / max(speed))), color = "green") +
    # geom_line(aes(y = grade*(max(vo2_abs) / max(grade))), color = "blue") +
    scale_y_continuous(
        name = "Absolute VO2",
        sec.axis = sec_axis(trans = ~./(max(df_avg$vo2_abs) / max(df_avg$grade)),
                            name = "Speed/Grade")
    )+
    theme_bw() +
    geom_smooth(method = "lm", aes(y = vo2_abs)) +
    ggtitle("Last 60 s with diff() showing when ΔVO2 <50 mL/O2/min\nfrompoint to point")

lm_last_min <- lm(vo2_abs ~ 1 + time, data = last_min)
summary(lm_last_min)
coef(lm_last_min)[2] * 60

# do you think they took the average


###### Last 30 s
last_30s_idx <- which((df_avg$time - max(df_avg$time) + 30) >= 0)
last_30s <- df_avg[last_30s_idx,]

ggplot(data = last_30s, aes(x = time)) +
    geom_point(aes(y = vo2_abs, color = plateau), alpha = 0.5) +
    # geom_line(aes(y = speed*(max(vo2_abs) / max(speed))), color = "green") +
    # geom_line(aes(y = grade*(max(vo2_abs) / max(grade))), color = "blue") +
    scale_y_continuous(
        name = "Absolute VO2",
        sec.axis = sec_axis(trans = ~./(max(df_avg$vo2_abs) / max(df_avg$grade)),
                            name = "Speed/Grade")
    )+
    theme_bw() +
    geom_smooth(method = "lm", aes(y = vo2_abs), se = FALSE) +
    ggtitle("Last 60 s with diff() showing when ΔVO2 <50 mL/O2/min\nfrompoint to point")

lm_last_30s <- lm(vo2_abs ~ 1 + time, data = last_30s)
summary(lm_last_30s)
coef(lm_last_30s)[2] * 60



.x <- df_avg$vo2_abs
plateau_func <- function(.x, cutoff = 150) {
    # browser()
    diff_x <- diff(.x)
    max_x <- max(.x)
    idx_max_x <- which(.x == max(.x))
    plateau <- logical(length = length(diff_x))
    for(i in 1:length(plateau)) {
        temp <- max_x - .x[i]
        if(temp < cutoff) {
            plateau[i] <- TRUE
        }
    }
    plateau <- c(NA, plateau)
    plateau
}

plateau <- plateau_func(df_avg$vo2_abs, cutoff = 150)
df_avg$plateau <- plateau

ggplot(data = df_avg, aes(x = time)) +
    geom_point(aes(y = vo2_abs, color = plateau), alpha = 0.5) +
    geom_line(aes(y = speed*(max(vo2_abs) / max(speed))), color = "green") +
    geom_line(aes(y = grade*(max(vo2_abs) / max(grade))), color = "blue") +
    scale_y_continuous(
        name = "Absolute VO2",
        sec.axis = sec_axis(trans = ~./(max(df_avg$vo2_abs) / max(df_avg$speed)),
                            name = "Speed/Grade")
    ) +
    theme_bw()

round(df_avg$grade / 0.5) * 0.5

grades <- c(0, 1, 2.5, 4, 5.5)

x <- df_avg$grade - grades

grade_factors <- function(.x, grades = c(0, 1, 2.5, 4.0, 5.5, 7.0)) {
    # browser()
    new_grades <- numeric(length = length(.x))
    for(i in 1:length(new_grades)) {
        temp <- abs(.x[i] - grades)
        new_grades[i] <- grades[which.min(temp)]
    }
    new_grades
}

actual_grades <- grade_factors(df_avg$grade)
df_avg <- df_avg %>% mutate(actual_grades = actual_grades)

lm1 <- lm(vo2_abs ~ 1 + as.factor(actual_grades), data = df_avg)
summary(lm1)

vo2_by_grade <- df_avg %>%
    group_by(actual_grades) %>%
    summarize(avg_vo2_abs = mean(vo2_abs))
diff(vo2_by_grade$avg_vo2_abs)

######
# The paper V̇O2max, protocol duration, and the V̇O2 plateau
# uses the last 30 s of a test and finds the VO2 vs. time slope
# to find the presence of plateau

last_min_idx <- which((df_avg$time - max(df_avg$time) + 60) >= 0)
last_min <- df_avg[last_min_idx,]
last_min$vo2_abs[nrow(last_min)] - last_min$vo2_abs[1] # 85, suggesting
# that there is a plateau if you use the original Taylor <= 150 mL, but not if you
# use the Robers <= 50 mL

last_min$vo2_abs[nrow(last_min)] - last_min$vo2_abs[30]
# if you use just the last 30 seconds, then we meet the <= 50 mL criteria

ggplot(data = df_avg, aes(x = time)) +
    geom_point(aes(y = peto2), color = "red")

ggplot(data = df_avg, aes(x = time)) +
    geom_point(aes(y = petco2), color = "blue")

(max(df_avg$vco2) / max(df_avg$ve))

ggplot(data = df_avg[df_avg[["speed"]] > 5.9,], aes(x = time)) +
    geom_point(aes(y = vo2_abs), color = "red") +
    geom_point(aes(y = vco2), color = "blue") +
    geom_point(aes(y = ve*(max(df_avg$vco2) / max(ve))), color = "green") +
    scale_y_continuous(
        name = "ve",
        sec.axis = sec_axis(trans = ~./(max(df_avg$vco2) / max(df_avg$ve)),
                            name = "VE")
    ) +
    theme_bw()




