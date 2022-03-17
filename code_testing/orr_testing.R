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

ggplot(data = df_unavg, aes(x = vco2, y = ve)) +
    geom_point(color = "orange", alpha = 0.5) +
    theme_bw()

df_avg <- avg_exercise_test(df_unavg, type = "breath", subtype = "rolling",
                  time_col = "time", roll_window = 9, roll_trim = 4)

ggplot(data = df_avg, aes(x = vo2_abs, y = vco2)) +
    geom_point(color = "purple", alpha = 0.5) +
    theme_bw()

.data <- df_avg
.x <- "vo2_abs"
.y <- "vco2"

lm_simple <- lm(vco2 ~ 1 + vo2_abs, data = df_avg)
summary(lm_simple)
anova(lm_simple)
# deviance(lm_simple)
RSS_simple <- sum(resid(lm_simple)^2)

# logLik(lm_simple)

ss_both <- loop_orr(df_avg, .x = .x, .y = .y)
plot(ss_both)
ss_min_idx <- which.min(ss_both)

df_left <- df_avg[1:ss_min_idx,]
df_right <- df_avg[ss_min_idx:nrow(df_avg),]

lm_left <- lm(vco2 ~ 1 + vo2_abs, data = df_left)
lm_right <- lm(vco2 ~ 1 + vo2_abs, data = df_right)

RSS_two <- sum(resid(lm_left)^2) + sum(resid(lm_right)^2)
MSE_two <- RSS_two / (nrow(df_avg) - 4) # -4 b/c estimating 4 parameters
F_stat <- (RSS_simple - RSS_two) / (2 * MSE_two)

pf(F_stat, df1 = 2, df2 = nrow(df_avg) - 4, lower.tail = FALSE)

pred_data <- tibble(vo2_abs = c(df_left$vo2_abs, df_right$vo2_abs),
                    vco2 = c(lm_left$fitted.values, lm_right$fitted.values))

plot_data_left <- tibble(vo2_abs = seq(min(df_avg$vo2_abs), max(df_avg$vo2_abs))) %>%
    mutate(vco2_left = predict(lm_left, newdata = .))

plot_data_right <- tibble(vo2_abs = seq(min(df_avg$vo2_abs), max(df_avg$vo2_abs))) %>%
    mutate(vco2_right = predict(lm_right, newdata = .))


ggplot(data = df_avg, aes(x = vo2_abs, y = vco2)) +
    geom_point(color = "purple", alpha = 0.5) +
    theme_bw() +
    geom_smooth(method = "lm", se = FALSE) +
    geom_line(data = pred_data[1:ss_min_idx,],
              aes(x = vo2_abs, y = vco2)) +
    geom_line(data = pred_data[ss_min_idx:nrow(pred_data),],
              aes(x = vo2_abs, y = vco2))

