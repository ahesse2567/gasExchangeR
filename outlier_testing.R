library(gasExchangeR)
library(tidyverse)

df <- read_csv("inst/extdata/mar16_101_pre.csv")

df <- df %>%
    trim_rest_rec(intensity_col = "speed", start_intensity = 3) %>%
    # vo2max_window() %>%
    # x_breath_mean(b = 4) %>%
    select(time, speed, grade, vo2, vo2.1, vco2, ve.btps, peto2, petco2) %>%
    rename(vo2_rel = vo2,
           vo2_abs = vo2.1,
           ve = ve.btps) %>%
    mutate(ve_vo2 = ve*1000 / vo2_abs,
           ve_vco2 = ve*1000 / vco2) %>%
    filter(speed > 5.6)

ggplot(data = df, aes(x = time)) +
    geom_point(aes(y = vo2_rel)) +
    theme_bw()

# library(ccTools)
# findOutlier(df, Q = 0.05, outliers = T)

df_knn <- exercise_outliers(
    df = df,
    method = "knn",
    vars = c("time", "speed", "grade", "vo2_abs", "vco2", "ve"),
    cutoff = 95,
    # passes = 1,
    k = 5,
    dist_method = "euclidean",
    p = 1.5,
    keep_outliers = TRUE)

ggplot(data = df_knn, aes(x = time, y = vo2_rel)) +
    geom_point(aes(color = outlier)) +
    theme_bw()


library(broom)
######################### Using prediction intervals
sd_to_pct <- function(sd) {
    pct <- abs(1 - pnorm(sd) * 2)
    pct
}

lm1 <- lm(vo2_abs ~ 1 + time, data = df)
au <- augment(lm1)

pred <- predict(lm1, interval = "prediction", level = sd_to_pct(2))
au <- bind_cols(au, .lwr = pred[,2], .upr = pred[,3])

au <- au %>%
    mutate(outlier = vo2_abs < .lwr | vo2_abs > .upr)

ggplot(data = au, aes(x = time)) +
    geom_point(aes(y = vo2_abs, color = outlier)) +
    geom_smooth(method = "lm", aes(y = vo2_abs)) +
    geom_line(aes(y = .lwr), linetype = "dashed") +
    geom_line(aes(y = .upr), linetype = "dashed") +
    theme_bw()

df <- bind_cols(df, au["outlier"])

ggplot(data = df, aes(x = time)) +
    geom_point(aes(y = vco2, color = outlier)) +
    geom_smooth(method = "lm", aes(y = vco2)) +
    theme_bw()

ggplot(data = df, aes(x = time)) +
    geom_point(aes(y = ve, color = outlier)) +
    geom_smooth(method = "lm", aes(y = ve)) +
    theme_bw()
