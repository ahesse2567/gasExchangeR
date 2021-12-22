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
    # browser()

    dist_method <- match.arg(dist_method,
                             choices = c("euclidean", "maximum", "manhattan",
                                         "canberra", "binary", "minkowski"))
    copy <- df
    cutoff <- round(cutoff,0)
    if(missing(passes)) {
        passes <- round(100 - cutoff)
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
        copy <- tibble::tibble(copy, closest_k = closest_k,)
        copy <- copy %>%
            dplyr::mutate(knn_outlier = closest_k > stats::quantile(closest_k, quant))
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
