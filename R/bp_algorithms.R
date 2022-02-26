#' Breakpoint Algorithm
#'
#' @param .data Gas exchange data.
#' @param algorithm The specific breakpoint algorithm you want to use to find a given threshold. Choices include \code{dmax}, \code{dmax_mod}, \code{jm}, \code{orr}, \code{v-slope}, \code{v_slope_simple}, and \code{splines}.
#' @param .x The x-axis variable
#' @param .y the y-axis variable
#' @param vo2 The name of the \code{vo2} variable
#' @param vco2 The name of the \code{vco2} variable
#' @param ve The name of the \code{ve} variable
#' @param time The name of the \code{time} variable
#'
#' @return
#' @export
#'
#' @references
#' Jones, R. H., & Molitoris, B. A. (1984). A statistical method for determining the breakpoint of two lines. Analytical Biochemistry, 141(1), 287–290. https://doi.org/10.1016/0003-2697(84)90458-5
#'
#' @examples
#' # Come up with an example later
#'
bp_algorithm <- function(.data,
                         algorithm,
                         .x,
                         .y,
                         vo2 = "vo2",
                         vco2 = "vco2",
                         ve = "ve",
                         time = "time") {
    stopifnot(!missing(.data), !missing(.x), !missing(.y), !missing(algorithm))
    # browser()

    algorithm <- match.arg(algorithm,
                           choices = c("dmax", "dmax_mod", "jm",
                                       "orr", "v-slope", "v_slope_simple",
                                       "splines")) # several ok = FALSE?

    class(.data) <- append(class(.data), algorithm) # could be useful to
    UseMethod("bp_algorithm", .data)
}

bp_algorithm.jm <- function(.data,
                            algorithm,
                            .x,
                            .y,
                            vo2 = "vo2",
                            vco2 = "vco2",
                            ve = "ve",
                            time = "time") {
    # browser()
    min_ss_idx <- which.min(loop(.data,
                          algorithm,
                          .x,
                          .y,
                          vo2 = vo2,
                          vco2 = vco,
                          ve = ve,
                          time = time))

    df_left <- .data[1:min_ss_idx,] # x < x0
    df_right <- .data[(min_ss_idx+1):nrow(.data),] # x >= x0

    x_knot <- .data[[.x]][min_ss_idx+1]

    df_right <- df_right %>%
        mutate(s1 = df_right[[.x]] - x_knot)

    lm_left <- lm(df_left[[.y]] ~ 1 + df_left[[.x]], data = df_left)
    # according to the JM method, the right regression line will have a constant equal to b0 + b1*x0
    b0_plus_b1x0 <- lm_left$coefficients[1] + lm_left$coefficients[2] * .data[[.x]][min_ss_idx+1]

    # for some reason, using I() function to try to force the regression to use b0_plus_b1x0 as an intercept doesn't work, while the offset argument does. The slope coefficient (b3) is the same, but the predicted values are not.
    lm_right <- lm(df_right[[.y]] ~ 0 + s1, data = df_right,
                   offset = rep(b0_plus_b1x0, nrow(df_right)))

    y_hat_left <- tibble(x = df_left[[.x]],
                         y_hat = lm_left$fitted.values,
                         method = "jm")
    y_hat_right <- tibble(x = df_right[[.x]],
                          y_hat = lm_right$fitted.values,
                          method = "jm")
    pred <- bind_rows(y_hat_left, y_hat_right)
    pct_slope_change <- 100*(lm_right$coefficients[1] - lm_left$coefficients[2]) /
        lm_left$coefficients[2]

    jm_row <- .data[min_ss_idx+1,] %>%
        select(time, vo2, .x, .y) %>%
        mutate(algorithm = algorithm)

    return(list(jm_row, pred, pct_slope_change))

}

#' @keywords internal
loop <- function(.data,
                 algorithm,
                 .x,
                 .y,
                 vo2 = "vo2",
                 vco2 = "vco2",
                 ve = "ve",
                 time = "time") {
    loop <- match.arg(algorithm, choices = c("jm", "orr", "v-slope",
                                         "v_slope_simple", "splines"))
    class(.data) <- append(class(.data), loop)
    UseMethod("loop", .data)
}

#' @keywords internal
loop.jm <- function(.data,
                    algorithm,
                    .x,
                    .y,
                    vo2 = "vo2",
                    vco2 = "vco2",
                    ve = "ve",
                    time = "time") {
    # browser()
    # rearranging data by the x variables was maybe a bad idea?
    # .data <- .data %>%
    #     arrange(.data[[.x]], .data[[time]])

    ss_both <- numeric(length = nrow(.data)-2)
    for(i in 2:(nrow(.data)-1)) {
        # algorithm data into left and right halves. x0 = i+1
        .data_left <- .data[1:i,] # x < x0
        .data_right <- .data[(i+1):nrow(.data),] # x >= x0

        # the JM formula is y = b0 + b1*x0 + b3(x-x0). Therefore, we will create a new column in .data_right that is equal to x - x0.
        .data_right <- .data_right %>%
            mutate(s1 = .data_right[[.x]] - .data[[.x]][i+1])

        .data_right$s1 = .data_right[[.x]] - .data[[.x]][i+1]

        lm_left <- lm(.data_left[[.y]] ~ 1 + .data_left[[.x]], data = .data_left)

        # according to the JM method, the right regression line will have a constant equal to b0 + b1*x0
        b0_plus_b1x0 <- lm_left$coefficients[1] + lm_left$coefficients[2] * .data[[.x]][i+1]

        # for some reason, using I() function to try to force the regression to use b0_plus_b1x0 as an intercept doesn't work, while the offset argument does. The slope coefficient (b3) is the same, but the predicted values are not.
        lm_right <- lm(.data_right[[.y]] ~ 0 + s1, data = .data_right,
                       offset = rep(b0_plus_b1x0, nrow(.data_right)))

        # This should be the same as using the offset argument, but it's not
        # lm_right <- lm(I(.data[[.y]] - b0_plus_b1x0) ~ 0 + .data[[.x]], data = .data_right)

        ss_left <- sum((lm_left$residuals)^2)
        ss_right <- sum((lm_right$residuals)^2)
        ss_both[i-1] <- (ss_left + ss_right)

        #ss_both <- c(ss_both, ss_left + ss_right)
    }

    ss_both
}
