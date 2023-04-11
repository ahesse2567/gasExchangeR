#' Construct Wasserman's Nine Panel Plot
#'
#' This is a famous representation of a cardiopulmonary exercise test often used in clinical practice.
#'
#' @param .data Gas exchange data frame or tibble
#' @param vt1_dat A data frame or tibble of values at VT1
#' @param vt2_dat A data frame or tibble of values at VT2/RC
#' @param vo2 Name of the `VO2` variable
#' @param vco2 Name of the `VCO2` variable
#' @param ve Name of the `pulmonary ventilation (VE)` variable
#' @param time Name of the `time` variable
#' @param hr Name of the `heart rate` variable
#' @param speed Name of the `speed` variable
#' @param grade Name of the `grade` variable
#' @param watts Name of the `watts` variable
#' @param vt Name of the `tidal volume (Vt)` variable
#' @param rer Name of the `respiratory exchange ratio (RER)` variable
#' @param peto2 Name of the `end tidal oxygen (PetO2)` variable
#' @param petco2 Name of the `end tidal carbon dioxide (PetCO2)` variable
#'
#' @returns a ggplot object
#' @export
#'
#' @examples
#' # TODO write an example
nine_panel_plot <- function(.data,
                            vt1_dat,
                            vt2_dat,
                            vo2 = "vo2",
                            vco2 = "vco2",
                            ve = "ve",
                            time = "time",
                            hr = "hr",
                            speed = "speed",
                            grade = "grade",
                            watts = "watts",
                            vt = "vt",
                            rer = "rer",
                            peto2 = "peto2",
                            petco2 = "petco2") {

    # plots <- vector(mode = "list", length = 9)
    # add main_x var that could be time or VO2?

    ve_vs_time_plot <- make_ve_vs_time_plot(
        .data = .data,
        vt1_dat = vt1_dat,
        vt2_dat = vt2_dat,
        ve = ve,
        time = time)

    hr_o2pulse_vs_time_plot <- make_hr_o2pulse_vs_time_plot(
        .data = .data,
        vt1_dat = vt1_dat,
        vt2_dat = vt2_dat,
        vo2 = vo2,
        hr = hr,
        time = time)

    vo2_vco2_intensity_plot <-  make_vo2_vco2_intensity_plot(.data,
                                                             vt1_dat,
                                                             vt2_dat,
                                                             vo2 = "vo2",
                                                             vco2 = "vco2",
                                                             time = "time",
                                                             speed = "speed",
                                                             grade = "grade",
                                                             watts = "watts")

    ve_vs_vco2_plot <- make_ve_vs_vco2_plot(.data = .data,
                                            vt1_dat = vt1_dat,
                                            vt2_dat = vt2_dat,
                                            vco2 = vco2,
                                            ve = ve)

    v_slope_hr_plot <- make_v_slope_hr_plot(.data = .data,
                                            vt1_dat = vt1_dat,
                                            vt2_dat = vt2_dat,
                                            vo2 = vo2,
                                            vco2 = vco2,
                                            hr = hr)

    vent_eqs_plot <- make_vent_eqs_vs_time_plot(.data = .data,
                                                vt1_dat,
                                                vt2_dat,
                                                vo2 = vo2,
                                                vco2 = vco2,
                                                ve = ve,
                                                time = time)

    vt_vs_ve_plot <- make_vt_vs_ve_plot(.data = .data,
                                        vt1_dat,
                                        vt2_dat,
                                        ve = ve,
                                        vt = vt)

    rer_vs_time_plot <- make_rer_vs_time_plot(.data = .data,
                                              vt1_dat,
                                              vt2_dat,
                                              time = time,
                                              rer = rer)

    end_tidals_vs_time_plot <- make_end_tidals_vs_time_plot(.data = .data,
                                                            vt1_dat,
                                                            vt2_dat,
                                                            time = time,
                                                            peto2 = peto2,
                                                            petco2 = petco2)

    plots <- gridExtra::grid.arrange(ve_vs_time_plot,
                                     hr_o2pulse_vs_time_plot,
                                     vo2_vco2_intensity_plot,
                                     ve_vs_vco2_plot,
                                     v_slope_hr_plot,
                                     vent_eqs_plot,
                                     vt_vs_ve_plot,
                                     rer_vs_time_plot,
                                     end_tidals_vs_time_plot,
                                     nrow = 3, ncol = 3)
    plots
}

####### Helper functions
#' Add threshold lines to a graph
#'
#' @keywords internal
#' @noRd
add_threshold_lines <- function(plt, x_var, vt1_dat = NULL, vt2_dat = NULL) {
    # extract information from input plot

    # leaving vt1_dat and vt2_dat equal to NULL lets you create the nine
    # panel plot without first needing to find any breakpoints
    # browser()
    plt_data <- ggplot2::ggplot_build(plt)

    # calculate horizontal adjustment factor based on range of x data
    h_just_factor <- plt_data$layout$panel_params[[1]]$x.range %>%
        diff() %>%
        {. * 0.015} # use 1% of x range

    # calculate height of threshold labels
    max_y <- purrr::map_dbl(plt_data[["data"]], ~ max(.[["y"]], na.rm = TRUE)) %>%
        max()
    # max_y <- numeric(length = length(plt_data[["data"]]))
    # for(i in seq_along(plt_data[["data"]])) {
    #     max_y[i] <- plt_data[["data"]][[i]]$y %>% max()
    # }
    # max_y <- max(max_y)

    # conditionally add VT1 line
    if(!is.null(vt1_dat) &&
       length(!is.na(vt1_dat[["breakpoint_data"]][[x_var]]) > 0) &&
       any(!is.na(vt1_dat[["breakpoint_data"]][[x_var]], na.rm = TRUE)) &&
       any(vt1_dat[["breakpoint_data"]][["determinant_bp"]], na.rm = TRUE)) {

        thresh_line_vt1 <- vt1_dat$breakpoint_data %>%
            dplyr::filter(est_ci == "estimate") %>%
            {.[[x_var]]}

        plt <- plt +
            ggplot2::geom_vline(xintercept = thresh_line_vt1) +
            ggplot2::geom_text(ggplot2::aes(
                x = thresh_line_vt1 - h_just_factor,
                y = max_y,
                label = "VT1"), angle = 90)
    }

    # conditionally add VT2 line
    if(!is.null(vt2_dat) &&
       length(!is.na(vt2_dat[["breakpoint_data"]][[x_var]]) > 0) &&
       any(!is.na(vt2_dat[["breakpoint_data"]][[x_var]], na.rm = TRUE)) &&
       any(vt2_dat[["breakpoint_data"]][["determinant_bp"]], na.rm = TRUE)) {

        thresh_line_vt2 <- vt2_dat$breakpoint_data %>%
            dplyr::filter(est_ci == "estimate") %>%
            {.[[x_var]]}

        plt <- plt +
            ggplot2::geom_vline(xintercept = thresh_line_vt2) +
            ggplot2::geom_text(ggplot2::aes(
                x = thresh_line_vt2  - h_just_factor,
                y = max_y,
                label = "VT2"), angle = 90)
    }

    plt
}

###### Plotting functions
# later add a method to recognize a breakpoint object?
#'
#'
#' @keywords internal
#' @noRd
make_ve_vs_time_plot <- function(.data,
                                 vt1_dat,
                                 vt2_dat,
                                 ve = "ve",
                                 time = "time") {

    out_plot <- ggplot2::ggplot(data = .data,
                                       ggplot2::aes(
                                           x = hms::hms(.data[[time]]),
                                           y = .data[[ve]])) +
        ggplot2::geom_point(color = "orange") +
        ggplot2::geom_line(color = "orange") +
        ggplot2::xlab("Time") +
        ggplot2::ylab("VE (L/min)") +
        ggplot2::theme_bw()

    out_plot <- add_threshold_lines(out_plot, vt1_dat, vt2_dat, time)
    out_plot
}

#' @keywords internal
#' @noRd
make_hr_o2pulse_vs_time_plot <- function(.data,
                                         vt1_dat,
                                         vt2_dat,
                                         vo2 = "vo2",
                                         hr = "hr",
                                         time = "time") {
    # check if heart rate column exists
    if (!hr %in% colnames(.data)) {
        message("Heart rate ", hr, " not found in input data frame.")
        return(ggplot2::ggplot(data = .data,ggplot2::aes(x = .data[[time]])))
    }

    scale_factor <-
        max(.data[[hr]]) / (max(.data[[vo2]]) / max(.data[[hr]]))

    out_plot <- ggplot2::ggplot(
        .data,
        ggplot2::aes(x = hms::hms(.data[[time]]))) +
        ggplot2::geom_point(ggplot2::aes(y = .data[[hr]],
                                         shape = "HR", color = "HR")) +
        ggplot2::geom_line(ggplot2::aes(y = .data[[hr]], color = "HR")) +
        ggplot2::geom_point(
            ggplot2::aes(y = (.data[[vo2]] / .data[[hr]]) * scale_factor,
                         shape = "O2Pulse", color = "O2Pulse")) +
        ggplot2::geom_line(
            ggplot2::aes(y = (.data[[vo2]] / .data[[hr]]) * scale_factor,
                         color = "O2Pulse")) +
        ggplot2::scale_y_continuous(
            name = "Heart Rate (BPM)",
            sec.axis = ggplot2::sec_axis(trans = ~ . / scale_factor,
                                         name = "Oxygen pulse (mL/min/beat)")
        ) +
        ggplot2::xlab("Time") +
        ggplot2::scale_color_manual(name = NULL, values = c("darkred", "black"),
                                    labels = c("Heart Rate", "Oxygen pulse"),
                                    guide = ggplot2::guide_legend(override.aes = list(
                                        shape = c(15, 16)))) +
        ggplot2::scale_shape_manual(name = NULL, values = c(15, 16),
                                    labels = c("Heart Rate", "Oxygen pulse"),
                                    guide = ggplot2::guide_legend(override.aes = list(
                                        color = c("darkred", "black")))) +
        ggplot2::theme_bw()

    out_plot <- add_threshold_lines(out_plot, vt1_dat, vt2_dat, time)
    out_plot
}

#' @keywords internal
#' @noRd
make_vo2_vco2_intensity_plot <- function(.data,
                                         vt1_dat,
                                         vt2_dat,
                                         vo2 = "vo2",
                                         vco2 = "vco2",
                                         time = "time",
                                         speed = "speed",
                                         grade = "grade",
                                         watts = "watts") {
    # make vo2 and vco2 vs. time plot
    # add speed, grade, and watts if they exist?

    out_plot <- ggplot2::ggplot(data = .data, ggplot2::aes(x = .data[[time]])) +
        ggplot2::geom_point(ggplot2::aes(y = .data[[vo2]], color = "vo2")) +
        ggplot2::geom_point(ggplot2::aes(y = .data[[vco2]], color = "vco2")) +
        ggplot2::geom_line(ggplot2::aes(y = .data[[vo2]], color = "vo2")) +
        ggplot2::geom_line(ggplot2::aes(y = .data[[vco2]], color = "vco2")) +
        ggplot2::scale_color_manual(name = NULL,
                                    values = c("vo2" = "red", "vco2" = "blue")) +
        ggplot2::xlab("Time") +
        ggplot2::ylab("VO2 & VCO2") +
        ggplot2::theme_bw()

    out_plot <- add_threshold_lines(out_plot, vt1_dat, vt2_dat, time)

    # it's going to be very difficult to plot both speed and grade on this
    # plot with their own y-axis scaling
    out_plot
}

#' @keywords internal
#' @noRd
make_ve_vs_vco2_plot <- function(.data,
                                 vt1_dat,
                                 vt2_dat,
                                 vco2 = "vco2",
                                 ve = "ve") {

    out_plot <- ggplot2::ggplot(data = .data,
                                ggplot2::aes(
                                    x = .data[[vco2]],
                                    y = .data[[ve]])) +
        geom_point(color = "brown") +
        ggplot2::xlab("VCO2") +
        ggplot2::ylab("VE (L/min)") +
        ggplot2::theme_bw()

    out_plot <- add_threshold_lines(out_plot, vt1_dat, vt2_dat, vco2)
    out_plot
}

#' @keywords internal
#' @noRd
make_v_slope_hr_plot <- function(.data,
                                 vt1_dat,
                                 vt2_dat,
                                 vo2 = "vo2",
                                 vco2 = "vco2",
                                 hr = "hr") {

    out_plot <- ggplot2::ggplot(data = .data, ggplot2::aes(x = .data[[vo2]])) +
        ggplot2::geom_point(ggplot2::aes(y = .data[[vco2]]), color = "blue") +
        ggplot2::xlab("VO2") +
        ggplot2::ylab("VCO2") +
        ggplot2::theme_bw()

    # conditional second axis??
    # check if heart rate column exists
    if (hr %in% colnames(.data)) {

        scale_factor <-
            (max(.data[[vco2]]) / max(.data[[hr]]))

        out_plot <- out_plot +
            ggplot2::geom_point(ggplot2::aes(y = .data[[vco2]], color = "vco2")) +
            ggplot2::geom_point(ggplot2::aes(y = .data[[hr]] * scale_factor,
                                             color = "hr")) +
            ggplot2::scale_y_continuous(
                name = "VCO2",
                sec.axis = ggplot2::sec_axis(trans = ~ . / scale_factor,
                                             name = "HR (BPM)")) +
            ggplot2::scale_color_manual(name = NULL,
                                        values = c("vco2" = "blue",
                                                   "hr" = "black"))
        # add other symbols
    }

    out_plot <- add_threshold_lines(out_plot, vt1_dat, vt2_dat, vo2)
    out_plot
}

#' @keywords internal
#' @noRd
make_vent_eqs_vs_time_plot <- function(.data,
                                       vt1_dat,
                                       vt2_dat,
                                       vo2 = "vo2",
                                       vco2 = "vco2",
                                       ve = "ve",
                                       time = "time") {
    # calculate ventilatory equivalents in case they aren't already in the df
    .data <- .data %>%
        dplyr::mutate(ve_vo2 = !!rlang::sym(ve) / !!rlang::sym(vo2) * 1000,
               ve_vco2 = !!rlang::sym(ve) / !!rlang::sym(vco2) * 1000)


    out_plot <- ggplot2::ggplot(data = .data,
                                ggplot2::aes(x = hms::hms(.data[[time]]))) +
        ggplot2::geom_point(ggplot2::aes(y = ve_vo2, color = "ve_vo2")) +
        ggplot2::geom_line(ggplot2::aes(y = ve_vo2, color = "ve_vo2")) +
        ggplot2::geom_point(ggplot2::aes(y = ve_vco2, color = "ve_vco2")) +
        ggplot2::geom_line(ggplot2::aes(y = ve_vco2, color = "ve_vco2")) +
        ggplot2::scale_color_manual(name = NULL,
                                    values = c("ve_vo2" = "purple",
                                               "ve_vco2" = "lightgreen")) +
        ggplot2::ylab("Ventilatory Equivalents") +
        ggplot2::xlab(stringr::str_to_sentence(time)) +
        ggplot2::theme_bw()

    out_plot <- add_threshold_lines(out_plot, vt1_dat, vt2_dat, time)
    out_plot
}

#' @keywords internal
#' @noRd
make_vt_vs_ve_plot <- function(.data,
                               vt1_dat,
                               vt2_dat,
                               ve = "ve",
                               vt = "vt") {

    if (!hr %in% colnames(.data)) {
        message("Heart rate ", hr, " not found in input data frame.")
        return(ggplot2::ggplot(data = .data,ggplot2::aes(x = .data[[ve]])))
    }

    out_plot <- ggplot2::ggplot(data = .data,
                                ggplot2::aes(x = .data[[ve]], y = .data[[vt]])) +
        ggplot2::geom_point(color = "darkslategray3") +
        ggplot2::xlab(paste(stringr::str_to_upper(ve), "(L/min)")) +
        ggplot2::ylab(stringr::str_to_title(vt)) +
        ggplot2::theme_bw()

    out_plot <- add_threshold_lines(out_plot, vt1_dat, vt2_dat, ve)
    out_plot
}

#' @keywords internal
#' @noRd
make_rer_vs_time_plot <- function(.data,
                                  vt1_dat,
                                  vt2_dat,
                                  time = "time",
                                  rer = "rer") {

    out_plot <- ggplot2::ggplot(data = .data,
                                ggplot2::aes(x = hms::hms(.data[[time]]),
                                             y = .data[[rer]])) +
        ggplot2::geom_point(color = "lightsalmon3") +
        ggplot2::geom_line(color = "lightsalmon3") +
        ggplot2::xlab(stringr::str_to_title(time)) +
        ggplot2::ylab(stringr::str_to_upper(rer)) +
        ggplot2::theme_bw() +
        ggplot2::scale_y_continuous(limits = c(0.65, max(.data[[rer]])))

    out_plot <- add_threshold_lines(out_plot, vt1_dat, vt2_dat, time)
    out_plot
}

#' @keywords internal
#' @noRd
make_end_tidals_vs_time_plot <- function(.data,
                                         vt1_dat,
                                         vt2_dat,
                                         time = "time",
                                         peto2 = "peto2",
                                         petco2 = "petco2") {

    scale_factor <- max(.data[[peto2]] / .data[[petco2]])

    out_plot <- ggplot2::ggplot(data = .data,
                                ggplot2::aes(x = hms::hms(.data[[time]]))) +
        ggplot2::geom_point(ggplot2::aes(y = .data[[peto2]], color = "PetO2")) +
        ggplot2::geom_point(ggplot2::aes(y = .data[[petco2]] * scale_factor,
                                         color = "PetCO2")) +
        ggplot2::geom_line(ggplot2::aes(y = .data[[peto2]], color = "PetO2")) +
        ggplot2::geom_line(ggplot2::aes(y = .data[[petco2]] * scale_factor,
                                         color = "PetCO2")) +
        ggplot2::scale_y_continuous(name = "PetO2",
                                    sec.axis = ggplot2::sec_axis(
                                        trans = ~ . / scale_factor,
                                    name = "PetCO2")) +
        ggplot2::scale_color_manual(name = NULL,
                           values = c("PetO2" = "mediumvioletred",
                                      "PetCO2" = "royalblue1")) +
        ggplot2::xlab(stringr::str_to_title(time)) +
        ggplot2::theme_bw()

    out_plot <- add_threshold_lines(out_plot, vt1_dat, vt2_dat, time)
    out_plot
}
