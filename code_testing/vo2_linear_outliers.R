library(gasExchangeR)
library(tidyverse)
library(broom)
library(gridExtra)

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
           ve_vco2 = ve*1000 / vco2,
           time = as.numeric(time)) %>%
    trim_pre_post(intensity_col = "grade",
                  pre_ex_intensity = 300.1,
                  post_ex_intensity = 300.1)

ggplot(data = df_unavg, aes(x = time)) +
    geom_point(aes(y = vo2_rel)) +
    theme_bw()

sd_to_pct <- function(sd) {
    pct <- abs(1 - pnorm(sd) * 2)
    pct
}

lm1 <- lm(vo2_abs ~ 1 + time, data = df_unavg)
au <- augment(lm1)

sd_lim <- 2

pred <- predict(lm1, interval = "prediction", level = sd_to_pct(sd_lim))
au <- bind_cols(au, .lwr = pred[,2], .upr = pred[,3])

au <- au %>%
    mutate(outlier = vo2_abs < .lwr | vo2_abs > .upr)

df_unavg <- bind_cols(df_unavg, au["outlier"])

vo2_plot <- ggplot(data = au, aes(x = time)) +
    geom_point(aes(y = vo2_abs, color = outlier)) +
    geom_smooth(method = "lm", aes(y = vo2_abs)) +
    geom_line(aes(y = .lwr), linetype = "dashed") +
    geom_line(aes(y = .upr), linetype = "dashed") +
    theme_bw() +
    ggtitle(paste("Local outlier removal of points ±",
                  sd_lim,
                  "sd from the local mean of VO2 vs. time"))

vco2_plot <- ggplot(data = df_unavg, aes(x = time)) +
    geom_point(aes(y = vco2, color = outlier)) +
    # geom_smooth(method = "lm", aes(y = vco2)) +
    theme_bw()

ve_plot <- ggplot(data = df_unavg, aes(x = time)) +
    geom_point(aes(y = ve, color = outlier)) +
    # geom_smooth(method = "lm", aes(y = ve)) +
    theme_bw()

grid.arrange(vo2_plot, vco2_plot, ve_plot, nrow = 3)
# my initial assessment is that vo2 and vco2 actually line up pretty well
# but there are points with ve that I would have flagged as an outlier. I also
# think that sometimes the assumption of linearity will not hold for vo2 vs. time
# for tests where someone has a long vo2 plateau.


ggplot(data = df_unavg, aes(x = time)) +
    geom_point(aes(y = vt, color = outlier)) +
    theme_bw()
