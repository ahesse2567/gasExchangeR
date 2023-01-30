#' Identify and remove outliers from exercise testing data
#'
#' @param df Gas exchange data
#' @param method Outlier identification method. Choose between KNN or SD prediction interval.
#' @param vars Variables names to include in the distance matrix for KNN
#' @param k How many neighbors for KNN
#' @param cutoff The quantile of data points to keep with KNN or the SD prediction interval level.
#' @param passes How many times to loop through the KNN outlier method. If ]
#' @param dist_method The distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given. See stats::dist.
#' @param p The power of the Minkowski distance. See stats::dist.
#' @param keep_outliers Should the function only identify outliers (\code{TRUE}) or remove them (\code{FALSE})?
#'
#' @return
#' @export
#'
#' @details The knn method normalizes the values that are passed to the distance matrix before calculating the distance matrix. They are normalized via \code{(.x - min(.x)) / (max(.x) - min(.x))}.
#'
#' @examples
<<<<<<< Updated upstream
#' # Add later
exercise_outliers <- function(
    df,
    method = c("knn"),
    vars = c("time", "speed", "grade", "vo2_abs", "vco2", "ve"),
    k = 5,
    cutoff = 95,
    passes,
    dist_method = "euclidean",
    p = 1.5,
    keep_outliers = TRUE) {
    stopifnot(!missing(df),
              k >= 1 & k %% 1 == 0,
              cutoff > 0 & cutoff <= 100,
              !missing(method))
    method = match.arg(method)
    class(df) <- append(class(df), method)
    UseMethod("exercise_outliers", df)
}
=======
#'
#' # TODO add an example
#'
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
        mutate(key = paste0(.data[[outlier_cols]],
                            "-",
                            1:length(.data[[outlier_cols]])))

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
        rolling_mos <- zoo::rollapply(data = copy_df[[outlier_cols]],
                                 width = width,
                                 FUN = get(func),
                                 f = mos,
                                 fill = NA,
                                 align = align)
        rolling_sd <- zoo::rollapply(data = copy_df[[outlier_cols]],
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

        outlier_idx <- copy_df[[outlier_cols]] < lower_lim |
            copy_df[[outlier_cols]] > upper_lim

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
            outlier_plot <- ggplot(data = copy_df, aes(x = .data[[time]],
                                                       y = .data[[outlier_cols]])) +
                geom_point(aes(color = outlier)) +
                geom_line(aes(y = lower_lim), linetype = "dashed") +
                geom_line(aes(y = upper_lim), linetype = "dashed") +
                ggtitle(glue::glue("Pass {n_passes + 1}, {length(which(outlier_idx))} outliers removed.")) +
                theme_minimal()
            plot(outlier_plot)
        }

        copy_df <- copy_df %>%
            filter(!outlier)

        n_passes <- n_passes + 1
        if(all(!outlier_idx)) any_outliers <- FALSE
>>>>>>> Stashed changes

#' @export
exercise_outliers.knn <- function(
    df,
    method = c("knn"),
    vars = c("time", "speed", "grade", "vo2_abs", "vco2", "ve"),
    k = 5,
    cutoff = 95,
    passes,
    dist_method = "euclidean",
    p = 1.5,
    keep_outliers = TRUE) {

    dist_method <- match.arg(dist_method,
                             choices = c("euclidean", "maximum", "manhattan",
                                         "canberra", "binary", "minkowski"))
    copy <- df
    cutoff <- round(cutoff,0)
    if(missing(passes)) {
        # passes <- round((100-cutoff)/100 * nrow(df)) # rm one per pass
        passes <- round(100 - cutoff) # rm more than one per pass
    } else if ((100 - cutoff) %% passes != 0) {
        # change the passes and cutoff if they aren't divisible
        pass_mult <- round((100 - cutoff) / passes, 0)
        change_to <- round((100 - cutoff) / pass_mult)
        diff <- change_to - passes
        passes <- passes + diff
        cutoff <- passes*pass_mult
    }
    quant <- (100-((100-cutoff)/passes))/100

    for(i in 1:passes) {
        d <- copy[vars] %>%
            dplyr::summarize_all(as.numeric) %>%
            dplyr::summarize_all(normalize) %>%
            stats::dist(method = dist_method) %>%
            as.matrix() %>%
            tibble::as_tibble() %>%
            purrr::map(sort) %>%
            dplyr::bind_rows()
        closest_k <- as.numeric(d[k+1,])
        copy <- tibble::tibble(copy, closest_k = closest_k)
        copy <- copy %>%
            dplyr::mutate(knn_outlier = closest_k > stats::quantile(closest_k,
                                                                    quant))
        copy <- copy %>%
            dplyr::filter(knn_outlier == FALSE) %>%
            dplyr::select(-closest_k)
    }
    cols_to_use <- intersect(colnames(df), colnames(copy))
    common <- match(do.call("paste", df[, cols_to_use]),
                    do.call("paste", copy[, cols_to_use]))
    df <- df %>%
        dplyr::mutate(outlier = dplyr::if_else(is.na(common) %in% as.integer(rownames(df)),
                                 TRUE,
                                 FALSE))
    if(keep_outliers == FALSE) {
        df <- df %>%
            dplyr::filter(outlier == FALSE)
    }
    df
}

#' @keywords internal
normalize <- function(.x) {
    out <- (.x - min(.x)) / (max(.x) - min(.x))
    out
}
