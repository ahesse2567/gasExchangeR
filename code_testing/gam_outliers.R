library(gasExchangeR)
library(tidyverse)
library(gam)
library(mgcv)

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

normalize <- function(.x) {
    out <- (.x - min(.x)) / (max(.x) - min(.x))
    out
}

df_unavg <- df_unavg %>%
    mutate(progress = normalize(df_unavg$vo2_abs) +
               normalize(df_unavg$vco2) +
               normalize(df_unavg$ve))

ggplot(data = df_unavg, aes(x = time)) +
    geom_point(aes(y = normalize(progress)), alpha = 0.5) +
    geom_smooth(aes(y = normalize(progress)), se = F) +
    theme_bw()


ggplot(data = df_unavg, aes(x = time)) +
    geom_point(aes(y = normalize(vo2_abs)), color = "red", alpha = 0.5) +
    geom_smooth(aes(y = normalize(vo2_abs)), color = "red", alpha = 0.5, se = F) +
    geom_point(aes(y = normalize(vco2)), color = "blue", alpha = 0.5) +
    geom_smooth(aes(y = normalize(vco2)), color = "blue", alpha = 0.5, se = F) +
    geom_point(aes(y = normalize(ve)), color = "orange", alpha = 0.5) +
    geom_smooth(aes(y = normalize(ve)), color = "orange", alpha = 0.5, se = F) +
    geom_point(aes(y = normalize(progress)), alpha = 0.5) +
    geom_smooth(aes(y = normalize(progress)), color = "black", se = F) +
    theme_bw()

?mgcv::gam




gam_prog <- mgcv::gam(vo2_abs ~  1 + s(ve) +s(vco2),
                     data = df_unavg, method = "REML")

summary(gam_prog)
gam_prog$fitted.values

df_unavg <- df_unavg %>%
    mutate(prog_fit = gam_prog$fitted.values)

ggplot(data = df_unavg, aes(x = time)) +
    geom_point(aes(y = vo2_abs), alpha = 0.5) +
    geom_line(aes(y = prog_fit)) +
    # geom_smooth(aes(y = progress), se= F) +
    theme_bw()


library(broom)
sd_to_pct <- function(sd) {
    pct <- abs(1 - pnorm(sd) * 2)
    pct
}

df_unavg$prog_fit
predict(gam_prog)

predict(gam_prog, interval = "prediction", level = sd_to_pct(2))

