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

bcp_ve <- bcp(y = cbind(df_inc$ve_vco2, df_inc$petco2))
plot(bcp_ve)

## there are wayyyy too many changepoints using bcp

#### bcp package for beysian change point
library(bcp)

##### univariate sequential data #####
# an easy problem with 2 true change points
set.seed(5)
x <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
bcp.1a <- bcp(x)
plot(bcp.1a, main="Univariate Change Point Example")
legacyplot(bcp.1a)

# a hard problem with 1 true change point
set.seed(5)
x <- rep(c(0,1), each=50)
y <- x + rnorm(50, sd=1)
bcp.1b <- bcp(y)
plot(bcp.1b, main="Univariate Change Point Example")

##### multivariate sequential data #####
# an easy problem in k=3 dimensions
set.seed(5)
x <- rnorm(6, sd=3)
y <- rbind(cbind(rnorm(50, x[1]), rnorm(50, x[2]), rnorm(50, x[3])),
           cbind(rnorm(50, x[4]), rnorm(50, x[5]), rnorm(50, x[6])))
bcp.2a <- bcp(cbind(y))
plot(bcp.2a, main="Multivariate (k=3) Change Point Example")
plot(bcp.2a, separated=TRUE, main="Multivariate (k=3) Change Point Example")

# a harder problem in k=5 dimensions
set.seed(5)
means1 <- rep(0, 5)
means2 <- rep(1, 5)
x <- rbind(matrix(rep(means1, each=50), nrow=50),
           matrix(rep(means2, each=50), nrow=50))
y <- x + rnorm(length(x), sd=1)
bcp.2b <- bcp(cbind(y))
plot(bcp.2b, main="Multivariate (k=5) Change Point Example")

##### linear models with sequential data #####
# 1 true change point at location 50; the predicting variable x is not related to location
x <- rnorm(100)
b <- rep(c(3,-3), each=50)
y <- b*x + rnorm(100)
bcp.3a <- bcp(y, x)
# in the two plots that follow, the location IDs are used as the plot characters
par(mfrow=c(1,2))
plot(y ~ x, type="n", main="Linear Regression: Raw Data")
text(x, y, as.character(1:100), col=(b/3)+2)
plot(y ~ x, type="n", main="Linear Regression: Posterior Means")
text(x, bcp.3a$posterior.mean[,1], as.character(1:100), col=(b/3)+2)
plot(bcp.3a, main="Linear Regression Change Point Example")

# 1 true change point at location 50; the predicting variable x is equal to location
x <- 1:100
b <- rep(c(3,-3), each=50)
y <- b*x + rnorm(100, sd=50)
bcp.3b <- bcp(y, x)
plot(bcp.3b, main="Linear Regression Change Point Example")

