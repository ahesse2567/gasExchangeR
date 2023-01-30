library(tidyverse)
# moving average filter
# in filter(), the filter parameter is the vector of filter coefficients
# each coefficient is the 'weight' of how much each value contributes to the
# filtered value

# see http://www.cookbook-r.com/Manipulating_data/Calculating_a_moving_average/

set.seed(993)
x <- 1:300
y <- sin(x/20) + rnorm(300,sd=.1)
y[251:255] <- NA

# Plot the unsmoothed data (gray)
plot(x, y, type="l", col=grey(.5))
# Draw gridlines
grid()

f20 <- rep(1/20, 20)
f20

y_lag <- stats::filter(y, f20, sides=1)
lines(x, y_lag, col="red")

# Smoothed symmetrically:
# average of current sample, 10 future samples, and 10 past samples (blue) (therefore rep amount is 21)
f21 <- rep(1/21,21)
f21
#>  [1] 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905
#>  [8] 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905
#> [15] 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905
y_sym <- stats::filter(y, f21, sides=2)
lines(x, y_sym, col="blue")


library(zoo)

#### Rolling but exclude middle value
# https://stackoverflow.com/questions/36603765/calculating-moving-average-excluding-middle-value-in-r

roll_func_exclude_middle <- function(x, f = "mean") {
    # browser()
    # if window is even, return average b/c there is no middle value
    if(length(x) %% 2 == 0L) {return(do.call(what = f, args = list(x = x)))}
    # -ceiling(0.5*length(x)) identifies and excludes the middle value
    if(length(x) %% 2 == 1L) {
        return(do.call(what = f,
                       args = list(x = x[-ceiling(0.5*length(x))])))
               }
}

# Fake data
set.seed(5)
values = cumsum(sample(1:5, size = 10, replace = TRUE))
values[4] <- 40
values
align <- "center"
width <-  5

rolling_mean <- rollapply(values, width = width, FUN = roll_func_exclude_middle,
                          f = mean, fill = NA, align = align)
rolling_sd <- rollapply(values, width = width, FUN = roll_func_exclude_middle,
                          f = sd, fill = NA, align = align)
sd_lim <- 3
rolling_upper_lim <- rolling_mean + rolling_sd * sd_lim
# rolling_upper_lim
rolling_lower_lim <- rolling_mean - rolling_sd * sd_lim
# rolling_lower_lim

tibble(x = values,
       rolling_mean = rolling_mean,
       rolling_sd = rolling_sd,
       lower_lim = rolling_lower_lim,
       upper_lim = rolling_upper_lim) %>%
    mutate(outlier = x < rolling_lower_lim | x > rolling_upper_lim)

# find values that are between the upper and lower limits
# also keep NA values at the beginning and the end
outlier_idx <- values < rolling_lower_lim |
    values > rolling_upper_lim # what if there's internal NA values?
outlier_idx

n_reassign <- if_else(align == "center", floor(width / 2), width - 1)

if(align == "center") {
    idx <- c(1:n_reassign,
             (length(outlier_idx) - n_reassign + 1):length(outlier_idx))
    outlier_idx[idx] <- FALSE
} else if (align == "left") {
    idx <- (length(outlier_idx) - n_reassign + 1):length(outlier_idx)
    outlier_idx[idx] <- FALSE
} else if (align == "right") {
    idx <- 1:n_reassign
    outlier_idx[idx] <- FALSE
}

values[outlier_idx] <- NA
values



# calculate rolling middle-excluded mean
# calculate rolling middle-excluded sd
# specify SD limit
# remove values outside mean ± SD limit
# needs a data frame or tibble
# specify which column you're using to remove outliers
# method argument, default == "rolling_breath" or something

library(readxl)
df_colnames <- read_xlsx("inst/extdata/Foreman_Ramp_10_18_2022.xlsx",
                         n_max = 0) %>%
    colnames()
df_raw <- read_xlsx("inst/extdata/Foreman_Ramp_10_18_2022.xlsx",
                    col_names = FALSE, skip = 3) %>%
    set_names(df_colnames)

df_unavg <- df_raw %>%
    janitor::clean_names() %>%
    rename(time = t) %>%
    mutate(time = lubridate::minute(time) * 60 + lubridate::second(time)) %>%
    filter(phase == "EXERCISE") %>%
    relocate(time, speed, grade)

.df <- df_unavg
outlier_col = "vo2"
width = 5
sd_lim = 3
mos = "mean"
align = "center"
exclude_test_val = TRUE
remove_outliers = TRUE
max_passes = Inf

ventilatory_outliers <- function(.df,
                                 outlier_col = "vo2",
                                 sd_lim = 3, # add option for conf level later
                                 width = 5,
                                 mos = "mean",
                                 align = "center",
                                 exclude_test_val = TRUE,
                                 remove_outliers = TRUE,
                                 max_passes = Inf) {
    # browser()
    # create key for flagging all outliers later
    .df <- .df %>%
        mutate(key = paste0(.df[[outlier_col]],
                            "-",
                            1:length(.df[[outlier_col]])))

    # set proper func if excluding the value to consider
    func <- if_else(exclude_test_val, "roll_func_exclude_middle", mos)

    # create data frame copy
    copy_df <- .df
    any_outliers <- TRUE
    n_passes <- 0
    while(any_outliers & n_passes < max_passes) {
        rolling_mos <- rollapply(data = .df[[outlier_col]],
                                 width = width,
                                 FUN = func,
                                 f = mos,
                                 fill = NA,
                                 align = align)
        rolling_sd <- rollapply(data = .df[[outlier_col]],
                                width = width,
                                FUN = sd,
                                fill = NA,
                                align = align)
        rolling_lower_lim <- rolling_mos - rolling_sd * sd_lim
        rolling_upper_lim <- rolling_mos + rolling_sd * sd_lim

        outlier_idx <- .df[[outlier_col]] < rolling_lower_lim |
            .df[[outlier_col]] > rolling_upper_lim

        n_reassign <- if_else(align == "center", floor(width / 2), width - 1)

        # TODO Should we try to use a method to determine if the head
        # and tail of the values are outliers? E.g., use the nearby values
        # to see if they are outliers? However, we need a way to keep the test
        # values from being used to calculate the mean and sd

        if(align == "center") {
            idx <- c(1:n_reassign,
                     (length(outlier_idx) - n_reassign + 1):length(outlier_idx))
            outlier_idx[idx] <- FALSE
        } else if (align == "left") {
            idx <- (length(outlier_idx) - n_reassign + 1):length(outlier_idx)
            outlier_idx[idx] <- FALSE
        } else if (align == "right") {
            idx <- 1:n_reassign
            outlier_idx[idx] <- FALSE
        }

        copy_df <- copy_df %>%
            mutate(outlier = outlier_idx)



    }

    # check which keys are still in the original .df

    if(remove_outliers) {
        .df <- .df %>%
            filter(!outlier)
    }

    .df

}

df_outliers <- df_unavg %>%
    # select(-phase) %>%
    ventilatory_outliers(outlier_col = "vo2",
                         sd_lim = 2, exclude_center_value = TRUE,
                         remove_outliers = FALSE)
