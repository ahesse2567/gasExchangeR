library(gasExchangeR)
library(tidyverse)

df_raw <- read_csv("inst/extdata/mar16_101_pre.csv")

df_unavg <- df_raw %>%
    trim_rest_rec(intensity_col = "speed") %>%
    rename(vo2_rel = vo2,
           vo2_abs = vo2.1,
           ve = ve.btps) %>%
    # vo2max_window() %>%
    # x_breath_mean(b = 4) %>%
    select(time, speed, grade, vo2_rel, vo2_abs, vco2, ve, peto2, petco2) %>%
    mutate(ve_vo2 = ve*1000 / vo2_abs,
           ve_vco2 = ve*1000 / vco2)

df_avg <- df_unavg %>%
    avg_exercise_test(type = "breath",
                      subtype = "rolling",
                      roll_window = 15)

ggplot(data = df_avg, aes(x = time, y = ve)) +
    geom_point(color = "orange", alpha = 0.5) +
    theme_bw()

df_inc <- df_avg %>%
    filter(speed > 5.6)

ggplot(data = df_inc, aes(x = time)) +
    geom_point(aes(y = vo2_abs/1000), color = "red", alpha = 0.5) +
    geom_point(aes(y = vco2/1000), color = "blue", alpha = 0.5) +
    # geom_point(aes(y = ve), color = "orange", alpha = 0.5) +
    theme_bw()

ggplot(data = df_inc, aes(x = time)) +
    geom_point(aes(y = petco2), color = "blue", alpha = 0.5) +
    geom_point(aes(y = ve_vco2), color = "green", alpha = 0.5) +
    ylim(c(20, 45)) +
    theme_bw()

## Nile data with one breakpoint: the annual flows drop in 1898
## because the first Ashwan dam was built
data("Nile")
plot(Nile)

## F statistics indicate one breakpoint
fs.nile <- Fstats(Nile ~ 1)
plot(fs.nile)
breakpoints(fs.nile)
lines(breakpoints(fs.nile))

## or
bp.nile <- breakpoints(Nile ~ 1)
summary(bp.nile)

## the BIC also chooses one breakpoint
plot(bp.nile)
breakpoints(bp.nile)

## fit null hypothesis model and model with 1 breakpoint
fm0 <- lm(Nile ~ 1)
fm1 <- lm(Nile ~ breakfactor(bp.nile, breaks = 1))
plot(Nile)
lines(ts(fitted(fm0), start = 1871), col = 3)
lines(ts(fitted(fm1), start = 1871), col = 4)
lines(bp.nile)

## confidence interval
ci.nile <- confint(bp.nile)
ci.nile
lines(ci.nile)

### strucchange with exercise data
plot(df_inc$ve_vco2)
fs.ve_vco2 <- Fstats(df_inc$ve_vco2 ~ 1)
plot(fs.ve_vco2)
breakpoints(fs.ve_vco2)
lines(breakpoints(fs.ve_vco2))


bp.ve_vco2 <- breakpoints(df_inc$ve_vco2 ~ 1)
summary(bp.ve_vco2)
plot(bp.ve_vco2)
breakpoints(bp.ve_vco2)

fm0 <- lm(df_inc$ve_vco2 ~ 1)
bp.ve_vco2 <- breakpoints(df_inc$ve_vco2 ~ 1)
fm1 <- lm(df_inc$ve_vco2 ~ breakfactor(bp.ve_vco2, breaks = 4))
plot(df_inc$ve_vco2)
lines(ts(fitted(fm0)), col = 3)
lines(ts(fitted(fm1)), col = 4)
lines(bp.ve_vco2)


### with petco2
plot(df_inc$petco2)
fs.petco2 <- Fstats(df_inc$petco2 ~ 1)
plot(fs.petco2)
breakpoints(fs.petco2)
lines(breakpoints(fs.petco2))


bp.petco2 <- breakpoints(df_inc$petco2 ~ 1)
summary(bp.petco2)
plot(bp.petco2)
breakpoints(bp.petco2)

fm0 <- lm(df_inc$petco2 ~ 1)
bp.petco2 <- breakpoints(df_inc$petco2 ~ 1)
fm1 <- lm(df_inc$petco2 ~ breakfactor(bp.petco2, breaks = 2))
plot(df_inc$petco2)
lines(ts(fitted(fm0)), col = 3)
lines(ts(fitted(fm1)), col = 4)
lines(bp.ve_vco2)
