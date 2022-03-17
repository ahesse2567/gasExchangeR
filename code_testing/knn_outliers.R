library(gasExchangeR)
library(tidyverse)

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


# library(ccTools)
# findOutlier(df, Q = 0.05, outliers = T)

df_knn <- exercise_outliers(
    df = df_unavg,
    method = "knn",
    vars = c("time", "speed", "grade", "vo2_abs", "vco2", "vt"),
    cutoff = 95,
    # passes = 1,
    k = 5,
    dist_method = "euclidean",
    p = 1.5,
    keep_outliers = TRUE)

vo2_plot <- ggplot(data = df_knn, aes(x = time)) +
    geom_point(aes(y = vo2_abs, color = outlier), alpha = 0.5) +
    theme_bw()

vco2_plot <- ggplot(data = df_knn, aes(x = time)) +
    geom_point(aes(y = vco2, color = outlier), alpha = 0.5) +
    theme_bw()

ve_plot <- ggplot(data = df_knn, aes(x = time)) +
    geom_point(aes(y = ve, color = outlier), alpha = 0.5) +
    theme_bw()

grid.arrange(vo2_plot, vco2_plot, ve_plot, nrow = 3)


