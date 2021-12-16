library(tidyverse)

df <- read_csv("inst/extdata/mar16_101_pre.csv")

df <- read_csv(test_to_read) %>%
    start_end_test() %>%
    vo2max_window() %>%
    # x_breath_mean(b = 4) %>%
    select(time, speed, grade, vo2, vo2.1, vco2, ve.btps, peto2, petco2) %>%
    rename(vo2_rel = vo2,
           vo2_abs = vo2.1,
           ve = ve.btps)

ggplot(data = df, aes(x = time, y = vo2_rel)) +
    geom_point() +
    theme_bw()

tr_mean <- function(.v, x, y) {
    # idx_col/decision_col = "ve_vo2_prod" should be implemented in the future
    # because Breeze uses the ve * vco2 to determine which rows to remove
    stopifnot(is.numeric(.v) | is.logical(.v),
              x >= 1 & x %% 1 ==0 & x %% 2 == 1,
              y >= 1 & y %% 1 ==0 & y %% 2 == 1 & y - x >= 2)
    # browser()
    num2remove <- y - x
    for (i in 1:(num2remove / 2)) {
        min_val <- which(.v == min(.v, na.rm = TRUE))
        if (length(min_val) > 1) { # if two values are tied, remove one at random
            min_val <- sample(min_val, 1)
        }
        .v <- .v[-min_val]
        max_val <- which(.v == max(.v, na.rm = TRUE))
        if (length(max_val) > 1) { # if two values are tied, remove one at random
            max_val <- sample(max_val, 1)
        }
        .v <- .v[-max_val]
        out <- mean(.v)
    }
    out
}

avg_exercise_test(df, type = "breath", subtype = "rolling")

my_func <- function() {
    browser()
    avg_exercise_test(df, type = "breath", subtype = "rolling")
}

