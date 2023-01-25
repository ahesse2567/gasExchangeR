#' Delete any 'resting/pre-exercise' and 'recovery/post-exercise' portions from a graded exercise test
#'
#' @param .data Unaveraged, breath-by-breath exercise test.
#' @param intensity_col Speed, wattage, or grade values that help determine start and end of a test.
#' @param pre_ex_intensity Intensity of 'pre-exercise' data.
#' @param post_ex_intensity Intensity that denotes the test has ended or entered the final recovery phase.
#'
#' @return
#' @export
#'
#' @examples
#' df <- data.frame(time = c(14, 18, 61, 78, 88, 100, 120, 150, 220, 231),
#' speed = c(0, 0, 0, 0, 3, 3, 3, 7.4, 7.4, 0))
#' trim_pre_post(.data = df, intensity_col = "speed")

# honestly I don't love how this function works rn
# my intuition is that I should enter the speeds below or above
# that should be excluded, not the speed that should be included...

trim_pre_post <- function(.data,
                          intensity_col,
                          pre_ex_intensity = 0,
                          post_ex_intensity = 0) {
    # browser()
    rle_start_intensity <- rle(.data[[intensity_col]] == pre_ex_intensity)
    start_idx <- rle_start_intensity$lengths[1] + 1 # index 1 gets first change,
    # but needs to add 1 to account for how rle gives you the index
    # of each run.

    rle_end_intensity <- rle(rev(.data[[intensity_col]] == post_ex_intensity))
    from_end_idx <- rle_end_intensity$lengths[1]
    end_idx <- nrow(.data) - from_end_idx
    if(end_idx < start_idx) {
        end_idx <- nrow(.data)
    }
    .data <- .data[start_idx:end_idx,]
    .data
}

