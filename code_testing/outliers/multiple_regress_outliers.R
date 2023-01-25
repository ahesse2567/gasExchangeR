library(gasExchangeR)
library(tidyverse)
library(car)
# library(broom)
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

time_mod <- lm(time ~ 1 + vo2_abs + vco2 + vt, data = df_unavg)
# autoplot(time_mod)
# influenceIndexPlot(time_mod)
plot(time_mod)

sd_to_pct <- function(sd) {abs(1 - pnorm(sd) * 2)}
sd_lim <- 2

pred <- predict(time_mod, interval = "prediction", level = sd_to_pct(sd_lim))
df_unavg <- df_unavg %>%
    mutate(.lwr = pred[,2],
           .upr = pred[,3],
           outlier = time < .lwr | vo2_abs > .upr)

# normalize <- function(.x) {(.x - min(.x)) / (max(.x) - min(.x))}

vo2_plot <- ggplot(data = df_unavg, aes(x = time)) +
    geom_point(aes(y = vo2_abs, color = outlier), alpha = 0.5) +
    theme_bw()

vco2_plot <- ggplot(data = df_unavg, aes(x = time)) +
    geom_point(aes(y = vco2, color = outlier), alpha = 0.5) +
    theme_bw()

vt_plot <- ggplot(data = df_unavg, aes(x = time)) +
    geom_point(aes(y = vt, color = outlier), alpha = 0.5) +
    theme_bw()

grid.arrange(vo2_plot, vco2_plot, vt_plot, nrow = 3)
