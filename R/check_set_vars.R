#################################################################
# This file contains internal functions that help set variables
# for the breakpoint() function
#################################################################

# Links I used to build these functions
# https://stackoverflow.com/questions/50688113/r-helper-function-that-checks-arguments-within-an-environment
# https://stackoverflow.com/questions/18541254/semi-automating-argument-validation-for-r-functions
# https://stackoverflow.com/questions/11885207/get-all-parameters-as-list
# https://stackoverflow.com/questions/15891459/get-a-list-of-all-function-parameters-from-inside-the-function

#' Check if user entered valid input, and assign function variables accordingly
#'
#' @keywords internal
#' @noRd
resolve_inputs <- function(inputs = c(as.list(environment(), all = TRUE))) {
    # check that the user input enough information to proceed

    # add inputs to local environment so I don't need to type inputs$xxx
    list2env(inputs, envir = environment())
    # Users need to input a minimum combination of methods and algorithms to continue
    stopifnot(!missing(.data),
              !all(is.null(method),
                   is.null(algorithm_vt1),
                   is.null(x_vt1),
                   is.null(y_vt1)),
              !all(is.null(method),
                   is.null(algorithm_vt2),
                   is.null(x_vt2),
                   is.null(y_vt2))
    )

    if(!is.null(method)) {
        inputs <- set_vars_by_method(inputs = inputs)
    }
    # set front_trim values to default if unspecified and using piecewise methods
    piecewise_methods <- c("jm", "orr", "v-slope", "spline_bp") # TODO add "v-slope_simple"

    if(is.null(inputs[["front_trim_vt1"]]) &
       inputs[["algorithm_vt1"]] %in% piecewise_methods) {
        inputs[["front_trim_vt1"]] <- 60
    }
    if(is.null(inputs[["front_trim_vt2"]]) &
       inputs[["algorithm_vt2"]] %in% piecewise_methods) {
        inputs[["front_trim_vt2"]] <- 60
    }
    inputs
}


#' Set x-y variables and algorithms according to user-specified methods
#'
#' @keywords internal
#' @noRd
set_vars_by_method <- function(inputs = c(as.list(environment(), all = TRUE))) {
    updated_inputs <- switch(inputs[["method"]],
                             "v-slope" = do.call(set_vars_v_slope,
                                                 args = list(inputs)),
                             "excess_co2" = do.call(set_vars_excess_co2,
                                                    args = list(inputs)),
                             "vent_eqs" = do.call(set_vars_vent_eqs,
                                                  args = list(inputs)),
                             stop("Invalid `method` value"))
    # future methods: "orr", "v-slope_simple"

    # orr: VE vs. VO2. This kinda finds RC first if the three-regression model is
    # best, but then only returns vt1. Given that, does this change the truncation
    # decision at RC?
    # V-slope

    updated_inputs
}


#' Set x-y variables and algorithms according to the v-slope method (Beaver et al., 1986)
#'
#' @keywords internal
#' @noRd
set_vars_v_slope <- function(inputs = c(as.list(environment(), all = TRUE))) {
    list2env(inputs, envir = environment())
    if(is.null(algorithm_vt1)) algorithm_vt1 <- "v-slope"
    if(is.null(x_vt1)) x_vt1 <- vo2
    if(is.null(y_vt1)) y_vt1 <- vco2
    if(is.null(algorithm_vt2)) algorithm_vt2 <- "jm" # TODO change this to the original V-slope 15% change method
    if(is.null(x_vt2)) x_vt2 <- vco2
    if(is.null(y_vt2)) y_vt2 <- ve
    if(is.null(front_trim_vt1)) front_trim_vt1 <- 60
    if(is.null(front_trim_vt2)) front_trim_vt2 <- 60

    # check if values match the original V-slope method
    if(algorithm_vt1 != "v-slope") warning("method = 'v-slope' but algorithm_vt1 is NOT set to v-slope")
    if(x_vt1 != vo2) warning("method = 'v-slope' but x variable of V-slope graph (VCO2 vs. VO2) is NOT set to vo2 variable")
    if(y_vt1 != vco2) warning("method = 'v-slope' but y variable of V-slope graph (VCO2 vs. VO2) is NOT set to vco2 variable")
    # if(algorithm_vt2 != "v-slope") warning("method = 'v-slope' but algorithm_vt1 is NOT set to v-slope") # TODO eventually change to 15% change method
    if(x_vt2 != vco2) warning("method = 'v-slope' but x variable of respiratory compensation point graph (VE vs. VCO2) is NOT set to vco2 variable")
    if(y_vt2 != ve) warning("method = 'v-slope' but y variable of respiratory compensation point graph (VE vs. VCO2) is NOT set to ve variable")
    if(front_trim_vt1 < 60) warning("method = 'v-slope' but front_trim_vt1 < 60 seconds from the beginning of the test before applying algorithms")
    if(front_trim_vt2 < 60) warning("method = 'v-slope' but front_trim_vt2 < 60 seconds from the beginning of the test before applying algorithms")

    rm("inputs") # remove old inputs
    return(c(as.list(environment(), all = TRUE))) # return new inputs
}


#' Set x-y variables according to the Excess CO2 method (Gaskill et al. 2001)
#'
#' @keywords internal
#' @noRd
set_vars_excess_co2 <- function(inputs = c(as.list(environment(),
                                                   all = TRUE))) {
    list2env(inputs, envir = environment())
    if(is.null(algorithm_vt1)) stop("The excess CO2 method must include a user-specified `algorithm_vt1`")
    if(any(is.null(algorithm_vt2),
           is.null(x_vt2),
           is.null(y_vt2)) & bp != "vt1") {
        stop("The Excess CO2 method only applies to VT1. Therefore, users must supply `algorithm_vt2`, `x_vt1`, and `y_vt2.`")
    }
    # assign null variables
    if(is.null(x_vt1)) x_vt1 <- time
    if(is.null(y_vt1)) {
        .data <- .data %>%
            dplyr::mutate(excess_co2 = .data[[vco2]]^2 /
                              .data[[vo2]] - .data[[vco2]])
        y_vt1 <- "excess_co2"
    }
    if(is.null(front_trim_vt1)) front_trim_vt1 <- 60

    # check if vco2 and vo2 are both in the same absolute units
    if(mean(.data[[vco2]] / .data[[vo2]]) > 2) {
        warning("VO2 and VCO2 columns are unlikely to both be in the same absolute units\nCheck if VO2 column indicated in arguments is absolute and that its units match VCO2")
    }

    # check for appropriate x_vt1 variables
    if(x_vt1 != time & x_vt1 != vo2) {
        warning("`x_vt1` is not set to an appropriate value for the excess co2 method. Appropriate x_vt1 values are time and VO2.")
    }
    if(front_trim_vt1 < 60) warning("method = 'excess_co2' but front_trim_vt1 < 60 seconds from the beginning of the test before applying algorithms")

    rm("inputs") # remove old inputs
    return(c(as.list(environment(), all = TRUE))) # return new inputs
}

#' Set x-y variables according to the ventilatory equivalents method
#'
#' @keywords internal
#' @noRd
set_vars_vent_eqs <- function(inputs = c(as.list(environment(),
                                                 all = TRUE))) {
    list2env(inputs, envir = environment())
    # check for algorithms
    if(is.null(algorithm_vt1)) stop("The ventilatory equivalents (`vent_eqs`) method must include a user-specified `algorithm_vt1`")
    if(is.null(algorithm_vt2)) stop("The ventilatory equivalents (`vent_eqs`) method must include a user-specified `algorithm_vt2`")
    # fill null values
    if(is.null(x_vt1)) x_vt1 <- time
    if(is.null(y_vt1)) {
        .data <- .data %>%
            dplyr::mutate(ve_vo2 = .data[[ve]] / .data[[vo2]] * 1000)
        y_vt1 <- "ve_vo2"
    }
    if(is.null(x_vt2)) x_vt2 <- time
    if(is.null(y_vt2)) {
        .data <- .data %>%
            dplyr::mutate(ve_vo2 = .data[[ve]] / .data[[vco2]] * 1000)
        y_vt2 <- "ve_vco2"
    }

    if(x_vt1 != time & x_vt1 != vo2) warning("`x_vt1` is not set to an appropriate value for the ventilatory equivalents (`vent_eqs`). Appropriate values are time and VO2.")
    if(mean(.data[[y_vt1]]) > 250 | mean(.data[[y_vt2]]) > 250) {
        warning("Abnormally high values in y variables: is minute ventilation in L/min and VO2 and VCO2 both in mL/min?")
    }

    rm("inputs") # remove old inputs
    return(c(as.list(environment(), all = TRUE))) # return new inputs
}
