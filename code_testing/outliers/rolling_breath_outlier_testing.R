library(tidyverse)
library(readxl)
library(devtools)
library(gasExchangeR)

#### Rolling but exclude middle value
# https://stackoverflow.com/questions/36603765/calculating-moving-average-excluding-middle-value-in-r

# much inspiration for this modeled after VO2Fitting
# https://www.mdpi.com/2075-4663/7/2/31
# Zacca, R., Azevedo, R., Figueiredo, P., Vilas-Boas, J. P., Castro, F. A. D. S., Pyne, D. B., & Fernandes, R. J. (2019). VO2FITTING: a free and open-source software for modelling oxygen uptake kinetics in swimming and other exercise modalities. Sports, 7(2), 31.


# method argument, default == "rolling_breath" or something??

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

df_unavg_no_outliers <- df_unavg %>%
    ventilatory_outliers(exclude_test_val = TRUE, width = 5,
                         align = "center", sd_lim = 3)

# df_unavg_no_outliers <- df_unavg %>%
#     ventilatory_outliers(outlier_cols = "vo2", time = "time",
#                          sd_lim = 3, width = 5, mos = "mean", align = "center",
#                          use_global_sd = TRUE, global_sd_mos = "median",
#                          exclude_test_val = TRUE, remove_outliers = FALSE,
#                          max_passes = Inf, plot_outliers = TRUE)

ggplot(data = df_unavg, aes(x = time, y = vo2)) +
    geom_point(color = "red", alpha = 0.5) +
    geom_line(color = "red", alpha = 0.5) +
    geom_point(data = df_unavg_no_outliers, color = "black", alpha = 0.5) +
    geom_line(data = df_unavg_no_outliers, color = "black", alpha = 0.5) +
    theme_minimal()





.data <- df_unavg
outlier_cols = "vo2"
time = "time"
width = 5
sd_lim = 3
mos = "mean"
align = "center"
exclude_test_val = TRUE
remove_outliers = FALSE
max_passes = Inf
plot_outliers = TRUE
use_global_sd = TRUE
global_sd_mos = "median"

# change to more general form that takes "align" variable?
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

ventilatory_outliers <- function(.data,
                                 outlier_cols = "vo2",
                                 time = "time",
                                 sd_lim = 3, # add option for conf level later
                                 width = 5,
                                 mos = "mean",
                                 align = "center",
                                 use_global_sd = TRUE,
                                 global_sd_mos = "median",
                                 exclude_test_val = TRUE,
                                 remove_outliers = TRUE,
                                 max_passes = Inf,
                                 plot_outliers = FALSE) {
    # confirm good user input
    mos <- match.arg(mos, choices = c("mean", "median"))
    align <- match.arg(align, choices = c("right", "left", "center"))
    global_sd_mos <- match.arg(global_sd_mos, choices = c("mean", "median"))

    # browser()

    # create key column for flagging all outliers later
    .data <- .data %>%
        mutate(key = paste0(.data[[outlier_col]],
                            "-",
                            1:length(.data[[outlier_col]])))

    # set proper func if excluding the value to consider
    func <- if_else(exclude_test_val, "roll_func_exclude_middle", mos)

    # create data frame copy
    copy_df <- .data
    any_outliers <- TRUE
    # if(is.null(max_passes))
    n_passes <- 0

    # TODO do I use lead/lag to allow exclude_test_val to work with
    # right and left-aligned rolling averages?

    while(any_outliers & n_passes < max_passes) {
        rolling_mos <- rollapply(data = copy_df[[outlier_col]],
                                 width = width,
                                 FUN = func,
                                 f = mos,
                                 fill = NA,
                                 align = align)
        rolling_sd <- rollapply(data = copy_df[[outlier_col]],
                                width = width,
                                FUN = sd,
                                fill = NA,
                                align = align)
        if(use_global_sd) {
            rolling_sd <- do.call(global_sd_mos,
                                 args = list(x = rolling_sd, na.rm = TRUE))
        }

        lower_lim <- rolling_mos - rolling_sd * sd_lim
        upper_lim <- rolling_mos + rolling_sd * sd_lim

        outlier_idx <- copy_df[[outlier_col]] < lower_lim |
            copy_df[[outlier_col]] > upper_lim

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

        if(plot_outliers) {
            outlier_plot <- ggplot(data = copy_df, aes(x = .data[[.time]],
                                       y = .data[[outlier_col]])) +
                geom_point(aes(color = outlier)) +
                geom_line(aes(y = lower_lim), linetype = "dashed") +
                geom_line(aes(y = upper_lim), linetype = "dashed") +
                theme_minimal()
            plot(outlier_plot)
        }

        copy_df <- copy_df %>%
            filter(!outlier)

        n_passes <- n_passes + 1
        if(all(!outlier_idx)) any_outliers <- FALSE

    }

    # check which keys are still in the original .data
    outliers <- if_else(.data$key %in% copy_df$key, FALSE, TRUE)
    .data <- .data %>%
        mutate(outlier = outliers)

    if(remove_outliers) {
        .data <- .data %>%
            filter(!outlier) %>%
            select(-outlier)
    }

    if(any(outliers)) {
        event <- if_else(remove_outliers, "removed", "detected")
        outs <- paste(which(outliers), collapse = ", ")
        print(glue::glue("Outliers {event} at indicies {outs}"))
    }

    .data

}

df_no_outliers <- df_unavg %>%
    ventilatory_outliers(outlier_cols = "vo2", time = "time",
                         sd_lim = 3, width = 5, mos = "mean", align = "center",
                         use_global_sd = TRUE, global_sd_mos = "median",
                         exclude_test_val = TRUE, remove_outliers = FALSE,
                         max_passes = Inf, plot_outliers = FALSE)



############ Testing just the rollapply components


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
upper_lim <- rolling_mean + rolling_sd * sd_lim
# upper_lim
lower_lim <- rolling_mean - rolling_sd * sd_lim
# lower_lim

tibble(x = values,
       rolling_mean = rolling_mean,
       rolling_sd = rolling_sd,
       lower_lim = lower_lim,
       upper_lim = upper_lim) %>%
    mutate(outlier = x < lower_lim | x > upper_lim)

# find values that are between the upper and lower limits
# also keep NA values at the beginning and the end
outlier_idx <- values < lower_lim |
    values > upper_lim # what if there's internal NA values?
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



