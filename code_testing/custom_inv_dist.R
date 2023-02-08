library(tidyverse)

normalize01 <- function(x, na.rm = TRUE) {
    return((x- min(x)) /(max(x)-min(x)))
}

dist_idx_conv <- function (i, j, dist_obj) {
    # convert a 2D to a 1D index
    # i = row (remember this starts at 2), j = column
    if (!inherits(dist_obj, "dist")) stop("please provide a 'dist' object")
    n <- attr(dist_obj, "Size")
    valid <- (i >= 1) & (j >= 1) & (i > j) & (i <= n) & (j <= n)
    k <- (2 * n - j) * (j - 1) / 2 + (i - j)
    k[!valid] <- NA_real_
    k
}

set.seed(34234)
tib <- tibble(x = seq(1, 10, by = 0.1)) %>%
    mutate(y = 0.5*x^2 - x + rnorm(n = n(), sd = 3),
           z = -1 * y + x^1.5)
tib

new_x <- 8
new_y <- 40

ggplot(data = tib, aes(x = x, y = y)) +
    geom_point() +
    geom_point(data = tibble(x = new_x, y = new_y),
               aes(x = x, y = y), color = "red")

ggplot(data = tib, aes(x = x, y = z)) +
    geom_point() +
    geom_vline(xintercept = new_x, color = "red")


# calculate distance from new point (1.5, 3)
d <- bind_rows(tibble(x = new_x, y = new_y),
               tib %>%
                   select(x, y)) %>%
    mutate(across(everything(), normalize)) %>%
    dist()

# add distances to tibble
tib <- tib %>%
    mutate(normalized_dist = d[dist_idx_conv(2:attr(d, "Size"), 1, d)],
           inv_weight = (normalized_dist / max(normalized_dist))^-1,
           pct_tot_weight = inv_weight / sum(inv_weight),
           weighted_contrib_z = pct_tot_weight * z)

new_z <- tib %>%
    summarize(new_val = sum(weighted_contrib_z)) %>%
    pull()
new_z

ggplot(data = tib, aes(x = x, y = z)) +
    geom_point() +
    geom_point(data = tibble(x = new_x, y = new_z),
               aes(x = new_x, y = new_z), color = "red") +
    geom_vline(xintercept = new_x, color = "red")

# interesting. When you make somewhat of an intentional outlier, this code
# gives you a reasonable value in the z dimension. I bet this would improve
# if there were more dimensions of data to consider.

