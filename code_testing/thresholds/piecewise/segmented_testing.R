library(gasExchangeR)
library(tidyverse)
library(devtools)
library(readxl)
library(segmented)

# df_unavg <- read_csv(
#     file.path("inst/extdata/anton_vo2max_clean.csv"),
#     show_col_types = FALSE)

# df_unavg <- read_xlsx(
#     file.path("../gasExchangeR_validation/data/interim/rand_15_cpet_exercisethresholds/mar22_131_post_gxt.xlsx")) %>%
#     janitor::clean_names() %>%
#     rename(vo2 = vo2_abs,
#            vo2_kg = vo2) %>%
#     filter(time < max(time) - 10)

df_colnames <- read_xlsx("inst/extdata/Foreman_Ramp_10_18_2022.xlsx",
                         n_max = 0) %>%
    colnames()
df_raw <- read_xlsx("inst/extdata/Foreman_Ramp_10_18_2022.xlsx",
                    col_names = FALSE, skip = 3) %>%
    set_names(df_colnames)

df_unavg <- df_raw %>%
    janitor::clean_names() %>%
    rename(time = t) %>%
    mutate(time = lubridate::minute(time) * 60 +
               lubridate::second(time),
           excess_co2 = (vco2^2 / vo2) - vo2) %>%
    filter(phase == "EXERCISE") %>%
    relocate(time, speed, grade) %>%
    ventilatory_outliers()

# remove outliers
df_unavg_no_outliers <- ventilatory_outliers(df_unavg,
                                             outlier_cols = "vo2",
                                             time = "time",
                                             sd_lim = 4,
                                             width = 5,
                                             mos = "mean",
                                             align = "center",
                                             use_global_sd = TRUE,
                                             global_sd_mos = "median",
                                             exclude_test_val = TRUE,
                                             remove_outliers = TRUE,
                                             max_passes = Inf,
                                             plot_outliers = FALSE)

df_avg <- avg_exercise_test(df_unavg, method = "time", calc_type = "bin",
                            time_col = "time", bin_w = 10)

ggplot(data = df_avg, aes(x = time, y = vo2)) +
    geom_point(alpha = 0.5, color = "red") +
    geom_point(aes(y = vco2), alpha = 0.5, color = "blue") +
    theme_bw()

ggplot(data = df_avg, aes(x = vo2, y = excess_co2)) +
    geom_point(alpha = 0.5, color = "red") +
    theme_bw()

ggplot(data = df_avg, aes(x = vco2, y = ve)) +
    geom_point(alpha = 0.5, color = "orange") +
    theme_bw()

ggplot(data = df_avg, aes(x = vo2, y = ve)) +
    geom_point(alpha = 0.5, color = "orange") +
    theme_bw()

ggplot(data = df_avg, aes(x = vo2)) +
    geom_point(aes(y = ve_vo2), alpha = 0.5, color = "green") +
    # geom_point(aes(y = ve_vco2), alpha = 0.5, color = "purple") +
    geom_smooth(data = df_avg[df_avg$time >= min(df_avg$time) + 60,],
                aes(x = vo2, y = ve_vo2),
                method = "lm",
                formula = y ~ poly(x, 5)) +
    theme_bw()

ggplot(data = df_avg[df_avg$time >= min(df_avg$time) + 60,],
       aes(x = vo2)) +
    geom_point(data = df_avg, aes(x = vo2, y = pet_co2),
               alpha = 0.5, color = "purple") +
    theme_bw() +
    geom_smooth(aes(x = vo2, y = pet_co2),
                method = "lm",
                formula = y ~ poly(x, 5))


reg_data <- df_avg[df_avg$time >= min(df_avg$time) + 60, ]

# start with simple linear regression
lm_simple_rcp <- lm(ve ~ 1 + vco2, data = reg_data)
plot(df_avg$vco2, df_avg$ve)
abline(lm_simple_rcp)
summary(lm_simple)
RSS_lm_simple_rcp <- sum(lm_simple_rcp$residuals^2)

# find respiratory compensation point w/ VE vs. VCO2
lm_segmented_rcp <- segmented(lm_simple_rcp,
                              seg.Z = ~ vco2,
                              npsi = 1)
lm_segmented_rcp
plot(lm_segmented_rcp)
points(df_avg$vco2, df_avg$ve)
abline(lm_simple_rcp, lwd = 2)
lines(lm_segmented_rcp)
points(lm_segmented_rcp)
RSS_lm_segmented_rcp <- sum(lm_segmented_rcp$residuals^2)

MSE_segmented_RCP <- RSS_lm_segmented_rcp / (nrow(reg_data) - 4)

# determine if RCP is determinant
crit_F <- qf(0.95, 1, nrow(reg_data) - 4, lower.tail = TRUE)
crit_F

f_stat <- (RSS_lm_simple_rcp - RSS_lm_segmented_rcp) /
    (2 * MSE_segmented_RCP)
f_stat
p_val_F <- stats::pf(f_stat, df1 = 2, df2 = nrow(reg_data) - 4,
                     lower.tail = FALSE)
p_val_F

# anova(lm_simple_rcp, lm_segmented_rcp)

trunc_idx <- which.min(abs(reg_data$vco2 - lm_segmented_rcp$psi[1,2]))
# vt1_reg_data <- reg_data[1:trunc_idx,]
vt1_reg_data <- reg_data
vo2_rcp <- reg_data[[trunc_idx, "vo2"]]
pct_vo2_max_rcp <- vo2_rcp / max(df_avg$vo2) * 100
pct_vo2_max_rcp

lm_simple_vt1 <- lm(vco2 ~ 1 + vo2, data = vt1_reg_data)
lm_segmented_vt1 <- segmented(lm_simple_vt1,
                              seg.Z = ~ vo2,
                              npsi = 1)
anova(lm_simple_vt1, lm_segmented_vt1)

vo2_vt1 <- lm_segmented_vt1$psi[1,2]
pct_vo2_max_vt1 <- vo2_vt1 / max(df_avg$vo2) * 100
pct_vo2_max_vt1

plot(lm_segmented_vt1)
lines(lm_segmented_vt1)
points(lm_segmented_vt1)
points(vt1_reg_data$vo2, vt1_reg_data$vco2)

plot(df_avg$vo2, df_avg$vco2)
# abline(v = vo2_rcp)
abline(v = vo2_vt1)

plot(df_avg$vo2, df_avg$ve_vo2)
abline(v = vo2_rcp)
abline(v = vo2_vt1)

plot(df_avg$vco2, df_avg$ve)
abline(v = vo2_rcp)
abline(v = vo2_vt1)

plot(df_avg$vo2, df_avg$ve_vco2)
abline(v = vo2_rcp)
abline(v = vo2_vt1)



# comparing 1 vs. two breakpoints w/ VE vs. VO2

lm_simple <- lm(ve ~ 1 + vo2, data = reg_data)
lm_segmented_1 <- segmented(lm_simple,
                            seg.Z = ~ vo2,
                            npsi = 1)
plot(lm_segmented_1)
points(df_avg$vo2, df_avg$ve)
lines(lm_segmented_1)

RSS_lm_segmented_1 <- sum(lm_segmented_1$residuals^2)
crit_F <- qf(0.95, 1, lm_segmented_1$df.residual, lower.tail = TRUE)

lm_segmented_2 <- segmented(lm_simple,
                            seg.Z = ~ vo2,
                            npsi = 2)
plot(lm_segmented_2)
RSS_lm_segmented_2 <- sum(lm_segmented_2$residuals^2)

MSE_lm_segmented_2 <- RSS_lm_segmented_2 / lm_segmented_2$df.residual

f_stat <- (RSS_lm_segmented_1 - RSS_lm_segmented_2) /
    (2 * MSE_lm_segmented_2)

p_val_F <- stats::pf(f_stat, df1 = 2, df2 = lm_segmented_1$df.residual,
                     lower.tail = FALSE)
p_val_F


# points(reg_data$vo2, lm_segmented_2$fitted.values, type = "p", col = "red")
plot(lm_segmented_2)
points(df_avg$vo2, df_avg$ve)
lines(lm_segmented_2)
confint(lm_segmented_2)

points(reg_data$vco2, lm_segmented_rcp$fitted.values, type = "l")
points(lm_segmented_rcp)
lines(lm_segmented_rcp)


