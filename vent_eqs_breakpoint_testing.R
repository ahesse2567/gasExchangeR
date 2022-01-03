library(gasExchangeR)
library(tidyverse)
library(splines)
library(AICcmodavg)
library(ggfortify)
library(car)

df <- read_csv("inst/extdata/37818_vo2max_unavg.csv")
df <- df %>%
    rename_with(.fn = tolower) %>%
    rename(vo2_rel = vo2...4,
           vo2_abs = vo2...5,
           ve = `ve btps`,
           vt = `vt btps`) %>%
    select(time, speed, grade, vo2_rel, vo2_abs, vco2, ve) %>%
    filter(speed > 3.5) %>% # after beginning steady state portion
    slice(1:which.max(ve))

df <- df %>%
    mutate(ve_vo2 = ve / (vo2_abs/1000),
           ve_vco2 = ve / (vco2/1000))

df_avg <- avg_exercise_test(df,
                            type = "breath",
                            subtype = "rolling",
                            roll_window = 9,
                            trim = 2)

ggplot(data = df_avg, aes(x = time)) +
    geom_point(aes(y = ve_vo2), color = "orange") +
    geom_point(aes(y = ve_vco2), color = "green") +
    theme_bw()

knot_max <- 100
bs_mod_list <- vector(mode = "list", length = length(seq_along(3:knot_max)))
for (i in 3:knot_max) {
    spline_mod_name <- paste0("bs.", i)
    spline_mod <- lm(ve_vo2 ~ 1 + bs(time, df = i + 3), data = df_avg)
    assign(spline_mod_name, spline_mod)
    bs_mod_list[[i - 2]] <- spline_mod
}

bs_mod_names <- character(length = length(3:knot_max))
for (i in 3:knot_max) {
    bs_mod_names[i - 2] <- paste("BS, ", i, " knots (df = ", i + 3, ")",
                                 sep = "")
}

bs_mod_comp <- aictab(cand.set = bs_mod_list, modnames = bs_mod_names) %>%
    as_tibble() %>%
    dplyr::select(1:4)

best_bs <- bs_mod_comp <- aictab(cand.set = bs_mod_list,
                                 modnames = bs_mod_names) %>%
    as_tibble() %>%
    dplyr::filter(AICc == min(AICc)) %>%
    dplyr::select(1:3)
best_bs
autoplot(bs.42, which = 1:6, ncol = 3, label.size = 3)
influenceIndexPlot(bs.42)
outlierTest(bs.42)

ggplot(data = df_avg, aes(x = time, y = ve_vo2)) +
    geom_point(color = "orange") +
    geom_smooth(method = "lm", formula = y ~ bs(x, 45)) +
    theme_bw()

plot_data = tibble(
    time = seq(from = min(df_avg$time),
                 to = max(df_avg$time),
                 by = 1)) %>%
    mutate(
        yhat = predict(lm_ve_vo2.6, newdata = .),
    )

ggplot(data = plot_data, aes(x = time, y = yhat)) +
    geom_line() +
    geom_point(data = df_avg, aes(x = time, y = ve_vo2), color = "orange") +
    theme_bw()

tidy(lm_ve_vo2)
lm_ve_vo2$coefficients
