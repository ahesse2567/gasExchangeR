library(tidyverse)
library(gasExchangeR)

x <- 1:100
s1 <- if_else(x <= round(length(x)/3), 0, x - round(length(x)/3))
s2 <- if_else(x <= round(length(x)*2/3), 0, x - round(length(x)*2/3))
y <- x + 4*s1 + 6*s2 + rnorm(x, sd = 10)

test_dat <- tibble(x = x, s1 = s1, s2 = s2, y = y)
# test_dat %>% View
ggplot(data = test_dat, aes(x = x, y = y)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    geom_smooth(method = "lm")

lm_simple <- lm(y ~ 1 + x, data = test_dat)
plot(lm_simple)

sign_runs <- sign(residuals(lm_simple)) %>% rle()
sign_runs


plot(residuals(lm_simple))
abline(h = 0)

ss <- loop_jm(.data = test_dat, .x = "x", .y = "y")
min_ss_idx <- which.min(ss)
plot(ss_models)

df_left <- test_dat[1:min_ss_idx,] # x < x0
df_right <- test_dat[(min_ss_idx):nrow(test_dat),] # x >= x0

x_knot <- test_dat[["x"]][min_ss_idx]

df_right <- df_right %>%
    dplyr::mutate(s1 = df_right[["x"]] - x_knot)

lm_left <- stats::lm(df_left[["y"]] ~ 1 + df_left[["x"]], data = df_left)
# according to the JM algorithm, the right regression line will have a constant equal to b0 + b1*x0
b0_plus_b1x0 <- lm_left$coefficients[1] + lm_left$coefficients[2] * test_dat[["x"]][min_ss_idx]

# for some reason, using I() function to try to force the regression to use b0_plus_b1x0 as an intercept doesn't work, while the offset argument does. The slope coefficient (b3) is the same, but the predicted values are not.
lm_right <- stats::lm(df_right[["y"]] ~ 0 + s1, data = df_right,
                      offset = rep(b0_plus_b1x0, nrow(df_right)))

resids <- c(residuals(lm_left), residuals(lm_right))

ggplot(data = test_dat, aes(x = x, y = y)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    geom_smooth(method = "lm") +
    geom_line(data = tibble(x = x, y_hat = c(lm_left$fitted.values[],
                                       lm_right$fitted.values[-1])),
               aes(y = y_hat))

plot(resids)
abline(h = 0)

# there's gotta be some sort of rle for how long the sign of the residuals
# is consistent

