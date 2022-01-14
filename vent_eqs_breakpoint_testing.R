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
    trim_pre_post(intensity_col = "grade",
                  pre_ex_intensity = 300.1,
                  post_ex_intensity = 300.1) %>%
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
bs_mod_comp

best_bs <- bs_mod_comp %>%
    dplyr::filter(AICc == min(AICc))

best_bs_idx <- which(names(bs_mod_list) == best_bs$Modnames)

autoplot(bs_mod_list[best_bs_idx], which = 1:6, ncol = 3, label.size = 3)
influenceIndexPlot(bs_mod_list[[best_bs_idx]])
outlierTest(bs_mod_list[[best_bs_idx]])

pred <- predict(bs_mod_list[[best_bs_idx]],
                interval = "prediction") %>%
    as_tibble()

plot_data <- df_unavg %>%
    select(time, ve_vo2) %>%
    bind_cols(pred) %>%
    mutate(outlier = ve_vo2 < lwr | ve_vo2 > upr)

ggplot(data = plot_data, aes(x = time, y = ve_vo2)) +
    geom_point(aes(color = outlier), alpha = 0.5) +
    geom_line(aes(y = fit)) +
    geom_line(aes(y = lwr), linetype = "dashed") +
    geom_line(aes(y = upr), linetype = "dashed") +
    theme_bw()



## with k-fold CV
library(boot)
library(lmvar)

# using glm()
set.seed(1234987)
k <- 10
knot_max <- 100
cv_error <- numeric(length = length(0:knot_max))
cv_mod_list <- vector(mode = "list", length = length(0:knot_max))
degree = 3
for(i in 0:knot_max) {
    spline_mod_name <- paste0("bs.", i)
    cv_mod_list[[i+1]] <- glm(ve_vo2 ~ 1 + bs(as.numeric(time),
                                             df = i + degree,
                                             degree = degree),
                             data = df_unavg,
                             family = )
    names(cv_mod_list)[i+1] <- spline_mod_name
    cv_error[i+1] <- cv.glm(df_unavg, glmfit = cv_mod_list[[i+1]], K = k)$delta[1]
    # is there a cv.lm somewhere?
}

# using lm()
set.seed(1234987)
k <- 10
knot_max <- 100
cv_error <- numeric(length = length(0:knot_max))
cv_mod_list <- vector(mode = "list", length = length(0:knot_max))
degree = 3
# this is REALLY slow compared to using glm and cv.glm()
for(i in 0:knot_max) {
    spline_mod_name <- paste0("bs.", i)
    cv_mod_list[[i + 1]] <- lm(
        ve_vo2 ~ 1 + bs(
            as.numeric(time),
            df = i + degree,
            degree = degree
        ),
        data = df_unavg,
        x = TRUE,
        y = TRUE
    )
    names(cv_mod_list)[i + 1] <- spline_mod_name
    cv_error[i + 1] <- lmvar::cv.lm(cv_mod_list[[i + 1]],
                                    k = k)$MSE$mean
}

cv_error
knot_tab <- tibble(knots = 0:knot_max,
                   df = 0:knot_max + degree,
                   cv_error = cv_error)
knot_tab
plot(knot_tab$cv_error)

best_cv <- which.min(cv_error)
best_cv # actually
knot_tab[best_cv,]

autoplot(cv_mod_list[best_cv], which = 1:6, ncol = 3, label.size = 3)
influenceIndexPlot(cv_mod_list[[best_cv]])
outlierTest(cv_mod_list[[best_cv]])
predict(cv_mod_list[[best_cv]], se.fit = TRUE)

## prediction errors are not built into a glm model
pred_cv <- predict(cv_mod_list[[best_cv]],
                interval = "prediction") %>%
    as_tibble()

plot_data <- df_unavg %>%
    select(time, ve_vo2) %>%
    bind_cols(pred_cv) %>%
    mutate(outlier = ve_vo2 < lwr | ve_vo2 > upr)

ggplot(data = plot_data, aes(x = time, y = ve_vo2)) +
    geom_point(aes(color = outlier), alpha = 0.5) +
    geom_line(aes(y = fit)) +
    geom_line(aes(y = lwr), linetype = "dashed") +
    geom_line(aes(y = upr), linetype = "dashed") +
    theme_bw()

