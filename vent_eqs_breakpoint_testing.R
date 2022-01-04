library(gasExchangeR)
library(tidyverse)
library(splines)
library(AICcmodavg)
library(ggfortify)
library(car)
library(broom)

df_raw <- read_csv("inst/extdata/37818_vo2max_unavg.csv")
df_unavg <- df_raw %>%
    rename_with(.fn = tolower) %>%
    rename(vo2_rel = vo2...4,
           vo2_abs = vo2...5,
           ve = `ve btps`,
           vt = `vt btps`) %>%
    select(time, speed, grade, vo2_rel, vo2_abs, vco2, ve) %>%
    trim_rest_rec(intensity_col = "grade",
                  pre_ex_intensity = 300.1,
                  end_ex_intensity = 300.1) %>%
    mutate(ve_vo2 = ve / (vo2_abs/1000),
           ve_vco2 = ve / (vco2/1000))

ggplot(data = df_unavg, aes(x = time)) +
    geom_point(aes(y = ve_vo2), color = "orange") +
    geom_point(aes(y = ve_vco2), color = "green") +
    theme_bw()

knot_max <- 100
bs_mod_list <- vector(mode = "list", length = length(seq_along(3:knot_max)))
for (i in 3:knot_max) {
    spline_mod_name <- paste0("bs.", i)
    bs_mod_list[[i-2]] <- lm(ve_vo2 ~ 1 + bs(time,
                                             df = i + 3,
                                             degree = 3),
                             data = df_unavg)
    names(bs_mod_list)[i-2] <- spline_mod_name
}

bs_mod_comp <- aictab(cand.set = bs_mod_list, modnames = names(bs_mod_list)) %>%
    as_tibble() %>%
    dplyr::select(1:4)

best_bs <- bs_mod_comp %>%
    dplyr::filter(AICc == min(AICc))

best_bs_idx <- which(names(bs_mod_list) == best_bs$Modnames)

autoplot(bs_mod_list[best_bs_idx], which = 1:6, ncol = 3, label.size = 3)
influenceIndexPlot(bs_mod_list[[best_bs_idx]])
outlierTest(bs_mod_list[[best_bs_idx]])

plot_data = tibble(
    time = seq(from = min(df_avg$time),
                 to = max(df_avg$time),
                 by = 1))

pred <- predict(bs_mod_list[[best_bs_idx]],
                interval = "prediction") %>%
    as_tibble()
bs_mod_list[[best_bs_idx]]$fitted.values

pred <- pred %>%
    mutate(outlier = if_else(fit < lwr | fit > upr, TRUE, FALSE))

pred$fit > pred$upr

plot_data <- df_unavg %>%
    select(time, ve_vo2) %>%
    bind_cols(pred)

ggplot(data = plot_data, aes(x = time, y = ve_vo2)) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = fit)) +
    geom_line(aes(y = lwr), linetype = "dashed") +
    geom_line(aes(y = upr), linetype = "dashed") +
    theme_bw()

ggplot(data = plot_data, aes(x = time, y = yhat)) +
    geom_line(color = "black") +
    geom_line(aes(y = .pred_lwr), linetype = "dashed", color = "purple") +
    geom_line(aes(y = .pred_upr), linetype = "dashed", color = "purple") +
    geom_point(data = df_unavg, aes(x = time, y = ve_vo2),
               color = "green", alpha = 0.5) +
    theme_bw()

ggplot(data = plot_data, aes(x = time, y = yhat)) +
    geom_line() +
    geom_point(data = df_unavg, aes(x = time, y = ve_vo2), color = "orange") +
    theme_bw()

tidy(lm_ve_vo2)
lm_ve_vo2$coefficients
