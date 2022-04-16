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

## with k-fold CV
library(boot)
library(lmvar)

# using glm()
set.seed(1234987)
k <- 5
knot_max <- nrow(df_unavg) / 2
cv_error <- numeric(length = length(0:knot_max))
cv_mod_list <- vector(mode = "list", length = length(0:knot_max))
degree = 3
for(i in 0:knot_max) {
    spline_mod_name <- paste0("bs.", i)
    cv_mod_list[[i+1]] <- glm(ve_vo2 ~ 1 + bs(as.numeric(time),
                                             df = i + degree,
                                             degree = degree),
                             data = df_unavg)
    names(cv_mod_list)[i+1] <- spline_mod_name
    cv_error[i+1] <- cv.glm(df_unavg, glmfit = cv_mod_list[[i+1]], K = k)$delta[1]
    # is there a cv.lm somewhere?
}

knot_tab <- tibble(knots = 0:knot_max,
                   df = 0:knot_max + degree,
                   cv_error = cv_error)
knot_tab
plot(knot_tab$cv_error)
plot(knot_tab$cv_error, ylim = c(1, 2.5))


best_cv <- which.min(cv_error)
best_cv # actually
knot_tab[best_cv,]

# I got around the lack of a good cv.lm() function by just redoing
# the best glm as an lm because lm is just the gaussian case of glm
# this allowed me to make a prediction interval
best_cv_mod <- lm(ve_vo2 ~ 1 + bs(as.numeric(time),
                                  df = knot_tab[[best_cv,"df"]],
                                  degree = degree),
                  data = df_unavg)

# autoplot(best_cv_mod, which = 1:6, ncol = 3, label.size = 3)
# influenceIndexPlot(cv_mod_list[[best_cv]])
# outlierTest(cv_mod_list[[best_cv]])
# predict(cv_mod_list[[best_cv]], se.fit = TRUE)

## prediction errors are not built into a glm model
sd_to_pct <- function(sd) {
    pct <- abs(1 - pnorm(sd) * 2)
    pct
}

sd_lim <- 3

pred_cv <- predict(best_cv_mod,
                interval = "prediction",
                level = sd_to_pct(sd_lim)) %>%
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
    theme_bw() +
    ggtitle(paste0("Local outlier removal of points ± ",
                   sd_lim,
                   " sd in non-linear\nrelationship using a b-spline model ",
                   k,
                   "-fold cv, and ",
                   best_cv - 1,
                   " knots"))

# now, refit the spline model but with the outliers removed


#### now using the one-standard-error rule

se <- function(.x) {sqrt(var(.x) / length(.x))}
se(cv_error)

plot(knot_tab$cv_error, ylim = c(1, 2.5))

cv_error_df <- tibble(cv_error = cv_error)

which((cv_error - min(cv_error) - se(cv_error)) < 0)[1]
one_se_cv_err <- which((cv_error - min(cv_error) - se(cv_error)) < 0)[1]

ggplot(data = cv_error_df, aes(x = as.integer(rownames(cv_error_df)),
                               y = cv_error)) +
    geom_point() +
    ylim(c(0.25, 2.5)) + # major outliers so I adjusted the axes
    theme_bw() +
    geom_hline(yintercept = min(cv_error), color = "red") +
    geom_hline(yintercept = min(cv_error) - se(cv_error),
               color = "red", linetype = "dashed") +
    geom_hline(yintercept = min(cv_error) + se(cv_error),
               color = "red", linetype = "dashed") +
    geom_vline(xintercept = one_se_cv_err, color = "blue")

knot_tab[one_se_cv_err,]

# I got around the lack of a good cv.lm() function by just redoing
# the best glm as an lm because lm is just the gaussian case of glm
# this allowed me to make a prediction interval
one_se_cv_mod <- lm(ve_vo2 ~ 1 + bs(as.numeric(time),
                                  df = knot_tab[[one_se_cv_err,"df"]],
                                  degree = degree),
                  data = df_unavg)

# autoplot(one_se_cv_mod, which = 1:6, ncol = 3, label.size = 3)
# influenceIndexPlot(cv_mod_list[[one_se_cv_err]])
# outlierTest(cv_mod_list[[one_se_cv_err]])
# predict(cv_mod_list[[one_se_cv_err]], se.fit = TRUE)

## prediction errors are not built into a glm model
sd_to_pct <- function(sd) {
    pct <- abs(1 - pnorm(sd) * 2)
    pct
}

lims <- 1:5
map(lims, sd_to_pct)

sd_lim <- 2

pred_one_se_cv <- predict(one_se_cv_mod,
                            interval = "prediction",
                            level = sd_to_pct(sd_lim)) %>%
    as_tibble()

plot_data_one_se <- df_unavg %>%
    select(time, ve_vo2) %>%
    bind_cols(pred_one_se_cv) %>%
    mutate(outlier = ve_vo2 < lwr | ve_vo2 > upr)

ggplot(data = plot_data_one_se, aes(x = time, y = ve_vo2)) +
    geom_point(aes(color = outlier), alpha = 0.5) +
    geom_line(aes(y = fit)) +
    geom_line(aes(y = lwr), linetype = "dashed") +
    geom_line(aes(y = upr), linetype = "dashed") +
    theme_bw() +
    ggtitle(paste0("Local outlier removal of points ± ",
                   sd_lim,
                   " sd in non-linear\nrelationship using a b-spline model ",
                   k,
                   "-fold cv, and ",
                   one_se_cv_err - 1,
                   " knots"))

### show multiple prediction intervals at once


pred_one_se_cv_3 <- predict(one_se_cv_mod,
                   interval = "prediction",
                   level = sd_to_pct(3)) %>%
    as_tibble()

plot_data_one_se_3 <- df_unavg %>%
    select(time, ve_vo2) %>%
    bind_cols(pred_one_se_cv_3) %>%
    mutate(outlier_3 = ve_vo2 < lwr | ve_vo2 > upr) %>%
    rename(fit3 = fit,
           lwr3 = lwr,
           upr3 = upr)

pred_one_se_cv_4 <- predict(one_se_cv_mod,
                            interval = "prediction",
                            level = sd_to_pct(4)) %>%
    as_tibble()

plot_data_one_se_4 <- df_unavg %>%
    select(time, ve_vo2) %>%
    bind_cols(pred_one_se_cv) %>%
    mutate(outlier_4 = ve_vo2 < lwr | ve_vo2 > upr) %>%
    rename(fit4 = fit,
           lwr4 = lwr,
           upr4 = upr)

plot_data_one_se_comb <- full_join(plot_data_one_se_3, plot_data_one_se_4) %>%
    mutate(outlier = case_when(outlier_3 == TRUE & outlier_4 == TRUE ~ "4sd",
                               outlier_3 == FALSE & outlier_4 == FALSE ~ "Not",
                               outlier_3 == FALSE & outlier_4 == TRUE ~ "3sd"))

ggplot(data = plot_data_one_se_comb, aes(x = time, y = ve_vo2)) +
    geom_point(aes(color = outlier), alpha = 0.5) +
    geom_line(aes(y = fit3)) +
    geom_line(aes(y = lwr3), linetype = "dashed") +
    geom_line(aes(y = upr3), linetype = "dashed") +
    geom_line(aes(y = fit4)) +
    geom_line(aes(y = lwr4), linetype = "dashed") +
    geom_line(aes(y = upr4), linetype = "dashed") +
    theme_bw()


# +
#     ggtitle(paste0("Local outlier removal of points ± ",
#                    sd_lim,
#                    " sd in non-linear\nrelationship using a b-spline model ",
#                    k,
#                    "-fold cv, and ",
#                    one_se_cv_err - 1,
#                    " knots"))
