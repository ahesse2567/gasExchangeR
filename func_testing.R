library(gasExchangeR)
library(tidyverse)

df_raw <- read_csv("inst/extdata/mar16_101_pre.csv")

df_unavg <- df_raw %>%
    rename_all(.funs = tolower) %>%
    rename(vo2_rel = vo2,
           vo2_abs = vo2.1,
           ve = ve.btps) %>%
    trim_pre_post(intensity_col = "speed") %>%
    trim_pre_post(intensity_col = "speed", pre_ex_intensity = 3) %>%
    trim_pre_post(intensity_col = "speed", pre_ex_intensity = 5.6) %>%
    select(time, speed, grade, vo2_rel, vo2_abs, vco2, ve, peto2, petco2) %>%
    mutate(ve_vo2 = ve*1000 / vo2_abs,
           ve_vco2 = ve*1000 / vco2)

ggplot(data = df_unavg, aes(x = time, y = vo2_abs)) +
    geom_point(color = "red", alpha = 0.5) +
    theme_bw()

df_avg <- avg_exercise_test(df_unavg, type = "breath", subtype = "bin",
                  time_col = "time", bin_w = 10, trim = 2)

ggplot(data = df_avg, aes(x = time, y = vo2_abs)) +
    geom_point(color = "red", alpha = 0.5) +
    theme_bw()

ggplot(data = df_avg, aes(x = vco2, y = ve)) +
    geom_point(color = "brown", alpha = 0.5) +
    theme_bw() +
    geom_smooth(method = "lm")

lm1 <- lm(time ~ petco2 + 1, data = df_avg)
lm1
out <- segmented::segmented(lm1)
out
plot(out)
points(x = df_avg$vco2, y = df_avg$ve_vco2)


ffset.seed(12)
xx<-1:100
zz<-runif(100)
yy<-2+1.5*pmax(xx-35,0)-1.5*pmax(xx-70,0)+15*pmax(zz-.5,0)+rnorm(100,0,2)
dati<-data.frame(x=xx,y=yy,z=zz)
out.lm<-lm(y~x,data=dati)

ggplot(data = dati, aes(x = x, y = y)) +
    geom_point()

#the simplest example: the starting model includes just 1 covariate
#.. and 1 breakpoint has to be estimated for that
o<-segmented(out.lm) #1 breakpoint for x
class(o)
plot(o)
