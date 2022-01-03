library(gasExchangeR)
library(tidyverse)

df <- read_csv("inst/extdata/37818_vo2max_unavg.csv")

trim_rest_rec <- function(.data,
                          intensity_col,
                          pre_ex_intensity,
                          end_intensity = 0) {
    browser()
    start <- which(dplyr::lag(diff(.data[[intensity_col]])) == pre_ex_intensity)
    # I think there's an issue here if the speed/grade doesn't repeat
    end <- which(diff(.data[[intensity_col]]) < end_intensity)
    if(length(end) == 0) { # test was terminated before clicking recovery on the computer
        end <- nrow(.data) # instead, set end to last data point in file
    }
    .data <- .data[start:end,]
    .data
}

df <- df %>%
    rename_with(.fn = tolower) %>%
    rename(vo2_rel = vo2...4,
           vo2_abs = vo2...5,
           ve = `ve btps`,
           vt = `vt btps`) %>%
    select(time, speed, grade, vo2_rel, vo2_abs, vco2, ve)

intensity_col <- "speed"
pre_ex_intensity <- 0
diff(df[[intensity_col]])
lag(df[[intensity_col]])


which(df[[intensity_col]] == pre_ex_intensity)
which(diff(df[[intensity_col]]) == pre_ex_intensity)

which(dplyr::lag(diff(df[[intensity_col]])) == pre_ex_intensity)

start <- which(dplyr::lag(diff(df[[intensity_col]])) == pre_ex_intensity)
# I think there's an issue here if the speed/grade doesn't repeat
end <- which(diff(.data[[intensity_col]]) < end_intensity)
if(length(end) == 0) { # test was terminated before clicking recovery on the computer
    end <- nrow(.data) # instead, set end to last data point in file
}
.data <- .data[start:end,]
.data


mara_df <- read_csv("inst/extdata/mar16_101_pre.csv")

mara_df <- mara_df %>%
    rename(vo2_rel = vo2,
           vo2_abs = vo2.1,
           ve = ve.btps) %>%
    select(time, speed, grade, vo2_rel, vo2_abs, vco2, ve, peto2, petco2) %>%
    mutate(ve_vo2 = ve*1000 / vo2_abs,
           ve_vco2 = ve*1000 / vco2)

mara_df2 <- mara_df[-nrow(mara_df),]
mara_df3 <- rbind(mara_df, mara_df[nrow(mara_df),])

pre_ex_intensity <- 0
x <- mara_df2[[intensity_col]] == pre_ex_intensity
rle_intensity <- rle(x)
start_idx <- rle_intensity$lengths[1] + 1 # index 1 gets first change,
start_idx

end_ex_intensity <- 0
y <- mara_df2[[intensity_col]] == end_ex_intensity
y
# diff(y)
# rev(y)
rle_y <- rle(rev(y))
rle_y
from_end_idx <- rle_y$lengths[1]
end_idx <- nrow(mara_df2) - from_end_idx
end_idx
if(end_idx < start_idx) {
    end_idx <- nrow(mara_df2)
}
# rev(diff(y))
# rle_y <- rle(rev(diff(y)))
# rle_y
cumsum(rle_y$lengths)


if(rle_y$lengths[1] == 1) {
    end_idx <- length(y) - 1
} else {
    end_idx <- length(y) - (rle_y$lengths[1] + 1)
}



rle_intensity <- df[[intensity_col]] %>%
    diff() %>%
    rle()
start_dix <- rle_intensity$lengths[1] + 2 # index 1 gets first change,
# but needs to add 2 to account for diff and for how rle gives you the index
# of each run. That is, if the first run is 5, those 5 points contributed
# to the pre-exercise data.

end_intensity <-  1.7
end_idx <- min(which(df[[intensity_col]] == end_intensity)) - 1

trim_rest_rec <- function(.data,
                          intensity_col,
                          pre_ex_intensity = 0,
                          end_ex_intensity = 0) {
    rle_start_intensity <- rle(.data[[intensity_col]] == pre_ex_intensity)
    start_idx <- rle_start_intensity$lengths[1] + 1 # index 1 gets first change,
    # but needs to add 1 to account for how rle gives you the index
    # of each run.

    rle_end_intensity <- rle(rev(.data[[intensity_col]] == end_ex_intensity))
    from_end_idx <- rle_end_intensity$lengths[1]
    end_idx <- nrow(.data) - from_end_idx
    if(end_idx < start_idx) {
        end_idx <- nrow(.data)
    }
    .data <- .data[start_idx:end_idx,]
    .data
}

tr_df <- trim_rest_rec(df,
                       intensity_col = "speed",
                       pre_ex_intensity = 0,
                       end_ex_intensity = 1.7)
