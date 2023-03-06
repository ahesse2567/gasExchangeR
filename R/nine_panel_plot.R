nine_panel_plot <- function(.data,
                            vt1_dat,
                            vt2_dat,
                            vo2 = "vo2",
                            vco2 = "vco2",
                            ve = "ve",
                            time = "time",
                            hr = "hr") {

    # plots <- vector(mode = "list", length = 9)

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

    ve_vs_vco2_plot <- make_ve_vs_vco2_plot(.data = .data,
                                            vt1_dat = vt1_dat,
                                            vt2_dat = vt2_dat,
                                            vco2 = vco2,
                                            ve = ve)

    plots <- gridExtra::grid.arrange(ve_vs_time_plot,
                                     hr_o2pulse_vs_time_plot,
                                     ve_vs_vco2_plot,
                                     nrow = 3, ncol = 3)
    plots
}

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
        ggplot2::geom_vline(xintercept = vt1_dat[["breakpoint_data"]][[time]]) +
        ggplot2::geom_vline(xintercept = vt2_dat[["breakpoint_data"]][[time]]) +
        ggplot2::theme_bw()
    out_plot
}

make_hr_o2pulse_vs_time_plot <- function(.data,
                                         vt1_dat,
                                         vt2_dat,
                                         vo2 = "vo2",
                                         hr = "hr",
                                         time = "time") {
    # check if heart rate column exists
    if (!hr %in% colnames(.data)) {
        message("Heart rate ", hr, " not found in input data frame.")
        return(ggplot(df, aes(x = time)))
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
        ggplot2::geom_vline(xintercept = vt1_dat[["breakpoint_data"]][[time]]) +
        ggplot2::geom_vline(xintercept = vt2_dat[["breakpoint_data"]][[time]]) +
        ggplot2::scale_color_manual(name = NULL, values = c("darkred", "black"),
                                    labels = c("Heart Rate", "Oxygen pulse"),
                                    guide = guide_legend(override.aes = list(
                                        shape = c(15, 16)))) +
        ggplot2::scale_shape_manual(name = NULL, values = c(15, 16),
                                    labels = c("Heart Rate", "Oxygen pulse"),
                                    guide = guide_legend(override.aes = list(
                                        color = c("darkred", "black")))) +
        ggplot2::theme_bw()

    out_plot
}

make_vo2_vco2_intensity_plot <- function(.data,
                                         vt1_dat,
                                         vt2_dat,
                                         vo2 = "vo2",
                                         vco2 = "vco2",
                                         speed = "speed",
                                         grade = "grade",
                                         watts = "watts",
                                         time = "time") {
    # make vo2 and vco2 vs. time plot
    # add speed, grade, and watts if they exist?

    out_plot <- ggplot2::ggplot(data = .data, ggplot2::aes(x = .data[[time]])) +
        ggplot2::geom_point(ggplot2::aes(y = .data[[vo2]], color = "vo2")) +
        ggplot2::geom_point(ggplot2::aes(y = .data[[vco2]], color = "vco2")) +
        ggplot2::geom_line(ggplot2::aes(y = .data[[vo2]], color = "vo2")) +
        ggplot2::geom_line(ggplot2::aes(y = .data[[vco2]], color = "vco2")) +
        ggplot2::scale_color_manual(name = NULL,
                                    values = c("vo2" = "red", "vco2" = "blue")) +
        xlab("Time") +
        ylab("VO2 & VCO2") +
        ggplot2::geom_vline(xintercept = vt1_dat[["breakpoint_data"]][[time]]) +
        ggplot2::geom_vline(xintercept = vt2_dat[["breakpoint_data"]][[time]]) +
        theme_bw()

    # it's going to be very difficult to plot both speed and grade on this
    # plot with their own y-axis scaling
    out_plot
}

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
        ggplot2::geom_vline(xintercept = vt1_dat[["breakpoint_data"]][[vco2]]) +
        ggplot2::geom_vline(xintercept = vt2_dat[["breakpoint_data"]][[vco2]]) +
        ggplot2::theme_bw()
    out_plot
}
