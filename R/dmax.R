#' Finding a breakpoint using Cheng's Dmax algorithm
#'
#' As decribed by Cheng (1992) with modifications for brevity, a third order curvilinear regressions is first calculated. Pilot work demonstrated that a third order regression produced higher correlations coefficients than first and second order regressions, with the use of greater than a third order regression not significantly increasing the correlation coefficient. Next, a straight line formed by the two end points in each curve was drawn. Then a formula for computing the distance from point to line was used to calculate the distance from each point on the curve to the straight line. The point yielding the maximal distance (Dmax) derived from the computation was taken as the threshold.
#'
#' @param .data Gas exchange data.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param bp Is this algorithm being used to find vt1 or vt2?
#' @param ... Dot dot dot mostly allows this function to work properly if breakpoint() passes arguments that is not strictly needed by this function.
#' @param vo2 The name of the vo2 column in \code{.data}
#' @param vco2 The name of the vco2 column in \code{.data}
#' @param ve The name of the ve column in \code{.data}
#' @param time The name of the time column in \code{.data}
#' @param alpha_linearity Significance value to determine if a piecewise model explains significantly reduces the residual sums of squares more than a simpler model.
#' @param ordering Prior to fitting any functions, should the data be reordered by the x-axis variable or by time? Default is to use the current x-axis variable and use the time variable to break any ties.
#' @param pos_change Do you expect the change in slope to be positive (default) or negative? If a two-line regression explains significantly reduces the sum square error but the change in slope does not match the expected underlying physiology, the breakpoint will be classified as indeterminate.
#' @param pos_slope_after_bp Should the slope after the breakpoint be positive? Default is `TRUE`. This catches cases when the percent change in slope is positive, but the second slope is still negative. Change to `FALSE` when PetCO2 is the y-axis variable.
#'
#' @return A list including slice of the original data frame at the threshold index with new columns `algorithm`, `determinant_bp`, `pct_slope_change`, `f_stat`, and `p_val_f.` The list also includes the fitted values, the left and right sides of the piecewise regression, and a simple linear regression.
#' @importFrom rlang :=
#' @export
#'
#' @references
#' Cheng, B., Kuipers, H., Snyder, A. C., Keizer, H. A., Jeukendrup, A., & Hesselink, M. (1992). A new approach for the determination of ventilatory and lactate thresholds. International Journal of Sports Medicine, 13(7), 518–522. https://doi.org/10.1055/s-2007-1021309
#'
#' @examples
#' # TODO write examples
dmax <- function(.data,
                 .x,
                 .y,
                 bp,
                 ...,
                 vo2 = "vo2",
                 vco2 = "vco2",
                 ve = "ve",
                 time = "time",
                 ordering = c("by_x", "time"),
                 alpha_linearity = 0.05,
                 pos_change = TRUE,
                 pos_slope_after_bp = TRUE
                 ){
    # the original paper has something about 50 mL increments in O2
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))

    ordering <- match.arg(ordering, several.ok = FALSE)
    .data <- order_cpet_df(.data, .x = .x , time = time,
                           ordering = ordering)

    # Get limits of x-axis for plots
    xmin = min(.data[[.x]], na.rm = T)
    xmax = max(.data[[.x]], na.rm = T)

    # Determine 3rd order polynomial
    g.model = stats::lm(.data[[.y]] ~ 1 + poly(.data[[.x]], 3, raw = TRUE), data = .data)

    # Get min and max y hats. This will be used to create a straight line that passes through these points
    y_hat_min <- min(g.model$fitted.values)
    y_hat_max <- max(g.model$fitted.values)

    x1 = xmin
    y1 = y_hat_min
    x2 = xmax
    y2 = y_hat_max

    g = g.model$coefficients  #Get polynomial coefficients
    names(g) = 0:(length(g) - 1) #Add vco2 name to coefficients

    ##### Determine D.max: Point in interval (x1,x2) that maximizes distance from line f to polynomial g
    f     = line.eq.2points(x1, y1, x2, y2) #Line equation
    D.max = ComputeDmax(f, g, x1, x2)
    g.D.max = poly.evaluate(g, D.max) #Find y-coord on polynomial g, corresponding to D.max

    # find closest observed data point to Dmax point for plotting left and right regressions?
    dmax_idx <- .data %>%
        dplyr::mutate(dist_x_sq = (.data[[.x]] - D.max)^2,
               dist_y_sq = (.data[[.y]] - g.D.max)^2,
               sum_sq = dist_x_sq + dist_y_sq) %>%
        dplyr::select(sum_sq) %>%
        dplyr::pull() %>%
        which.min()

    # this needs to be fixed. This does NOT constrain the left line to
    # go from the first data point on the curve to the dmax point. It
    # also does not currently constrain the right line to go from the
    # dmax point to the last point on the curve
    df_left <- .data[1:dmax_idx,]
    # make linear models of the two halves
    lm_left <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::as.formula() %>%
        stats::lm(data = df_left)

    df_right <- .data[dmax_idx:nrow(.data),]
    lm_right <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::as.formula() %>%
        stats::lm(data = df_right)

    pct_slope_change <- 100*(lm_right$coefficients[2] - lm_left$coefficients[2]) /
        lm_left$coefficients[2]

    lm_simple <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::as.formula() %>%
        stats::lm(data = .data)

    pw_stats <- piecewise_stats(lm_left, lm_right, lm_simple)
    list2env(pw_stats, envir = environment())

    determinant_bp <- check_if_determinant_bp(p = pf_two,
                                              pct_slope_change = pct_slope_change,
                                              pos_change = pos_change,
                                              pos_slope_after_bp =
                                                  pos_slope_after_bp,
                                              slope_after_bp = coef(lm_right)[2],
                                              alpha = alpha_linearity)

    y_hat_left <- tibble::tibble("{.x}" := df_left[[.x]],
                         "{.y}" := lm_left$fitted.values,
                         algorithm = "dmax")
    y_hat_right <- tibble::tibble("{.x}" := df_right[[.x]],
                          "{.y}" := lm_right$fitted.values,
                          algorithm = "dmax")
    pred <- dplyr::bind_rows(y_hat_left, y_hat_right)

    # find closest actual data point to dmax point and return data
    bp_dat <- find_threshold_vals(.data = .data, thr_x = D.max, thr_y = g.D.max,
                                  .x = .x, .y = .y)
    bp_dat <- bp_dat %>%
        dplyr::mutate(bp = bp,
               algorithm = "dmax",
               x_var = .x,
               y_var = .y,
               determinant_bp = determinant_bp,
               pct_slope_change = pct_slope_change,
               f_stat = f_stat,
               p_val_f = pf_two) %>%
        dplyr::relocate(bp, algorithm, x_var, y_var, determinant_bp,
                 pct_slope_change, f_stat, p_val_f)

    bp_plot <- make_piecewise_bp_plot(.data, .x, .y, lm_left, lm_right, bp_dat)

    bp_plot <- bp_plot +
        ggplot2::geom_line(ggplot2::aes(y = g.model$fitted.values), linetype = "dashed")

    return(list(breakpoint_data = bp_dat,
                fitted_vals = pred,
                poly_model = g.model,
                dmax_point = c("x" = D.max, "y" = g.D.max),
                lm_left = lm_left,
                lm_right = lm_right,
                lm_simple = lm_simple,
                bp_plot = bp_plot))
}


# Read below for an implementation of the Dmax function by Yuri Mejia Miranda.
########################################################################################
#Goal: 1) Functions to handle polynomials such as sum, product, composition, derivative, evaluation, etc.
#      2) Dmax algorithm contained in ComputeDmax function
#      These functions are used in the calculation of Dmax and Dmax modified
#Author: Yuri Mejia Miranda, Mathematician CHDR
#Date:   20160607
########################################################################################

#distance.2points: Given two points (x1, y1) and (x2, y2), it returns their Eucllidean distance (squared)
# INPUT:
# x1,y1,x2,y2 = (numeric) x,y-coordinates of points 1 and 2, respectively
# OUTPUT:
# squared distance between (x1, y1) and (x2, y2)
#' @keywords internal
distance.2points = function(x1 ,y1, x2, y2){
    return((x1 - x2)^2 + (y1 - y2)^2)
}

#poly.sum: Gets coefficients of two polynomials and returns coefficients of their sum.
#INPUT:
# p1, p2 = (numeric vectors) contain coefficients of polynomial in increasing order
#           Example: p(x) = a0 + a1*x + a2*x^2 + ... + an*x^n is entered as c(a0, a1, a2,..., an)
#OUTPUT:
# Returns (numeric vector) with coefficients of their sum p1(x)+p2(x) in increasing order
#' @keywords internal
poly.sum = function(p1, p2){
    degree = max(length(p1), length(p2))
    p1.new = c(p1, rep(0, degree - length(p1)))
    p2.new = c(p2, rep(0, degree - length(p2)))
    sum = p1.new + p2.new
    names(sum) = 0:(length(sum)-1)
    return(sum)
}

#poly.product: Gets coefficients of two polynomials and returns coefficients of their product.
#INPUT:
# p1, p2 = (numeric vectors) contain coefficients of polynomial in increasing order
#                Example: p(x) = a0 + a1*x + a2*x^2 + ... + an*x^n is entered as c(a0, a1, a2,..., an)
#OUTPUT:
# Returns (numeric vector) with coefficients of their product p1(x)*p2(x) in increasing order
#' @keywords internal
poly.product = function(p1, p2){
    degree2 = length(p2) - 1

    product = c(p1 * p2[1], rep(0, degree2))
    if (degree2 >= 1){
        for (i in 1:degree2){
            product1 = c(rep(0, i), p1 * p2[i+1], rep(0, degree2 - i))
            product  = poly.sum(product, product1)
        }
    }
    names(product) = 0:(length(product)-1)
    return(product)
}

#poly.derivative: Gets coefficients of a polynomial and returns coefficients of its derivative.
#INPUT:
# p = (numeric vector) contains coefficients of polynomial in increasing order:
#     Example: p(x) = a0 + a1*x + a2*x^2 + ... + an*x^n is entered as c(a0, a1, a2,..., an)
#OUTPUT:
# Returns (numeric vector) with coefficients of derivative p'(x) in increasing order
#' @keywords internal
poly.derivative = function(p){
    degree = length(p) - 1
    if (degree > 0){
        derivative = p * 0:degree
        derivative = derivative[2:(degree + 1)]
    }else{
        derivative = 0
    }
    names(derivative) = 0:(length(derivative)-1)
    return(derivative)
}

#poly.power: Gets coefficients of a polynomial and integer power,
#returns coefficients of polynomial raised to power.
#INPUT:
# p = (numeric vector) contains coefficients of polynomial in increasing order:
#     Example: p(x) = a0 + a1*x + a2*x^2 + ... + an*x^n is entered as c(a0, a1, a2,..., an)
# n = (non-negative integer) power to which polynomial is raised
#OUTPUT:
# Returns (numeric vector) with coefficients of polynomial raised to power p(x)^n in increasing order
#' @keywords internal
poly.power = function(p, n){
    if(length(which(p != 0)) > 0){ #if polynomial is NOT zero polynomial
        if (n == 0){ #if power is 0
            power = 1
        }else{ #if power >= 1
            power = p
            if(n > 1){
                for(i in 2:n){
                    power = poly.product(power, p)
                }
            }
        }
    }else if(length(which(p != 0)) == 0){ #else, polynomial is zero polynomial
        power = 0
    }
    names(power) = 0:(length(power)-1)
    return(power)
}

#poly.evaluate: Gets coefficients of a polynomial p(x) and a numerical value x0, and
#               returns the  numeric value of the polynomial evaluated at this point p(x0) .
#INPUT:
# p  = (numeric vector) contains coefficients of polynomial in increasing order:
#      Example: p(x) = a0 + a1*x + a2*x^2 + ... + an*x^n is entered as c(a0, a1, a2,..., an)
# x0 = (numeric) value at which polynomial is evaluated
#OUTPUT:
# Returns (numeric) value p(x0)
#' @keywords internal
poly.evaluate = function(p, x0){
    degree = length(p) - 1
    x0.powers = x0^(0:degree)
    return( sum(p * x0.powers) )
}

#poly.composition: Gets coefficients of two polynomials and returns coefficients of their composition.
#INPUT:
# p1, p2 = (numeric vectors) contain coefficients of polynomial in increasing order
#           Example: p(x) = a0 + a1*x + a2*x^2 + ... + an*x^n is entered as c(a0, a1, a2,..., an)
#OUTPUT:
# Returns (numeric vector) with coefficients of their composition p1(p2(x)) in increasing order
#' @keywords internal
poly.composition = function(p1, p2){
    degree1 = length(p1) - 1
    compose = 0

    for(i in 0:degree1){
        compose = poly.sum(compose, p1[i+1] * poly.power(p2, i))
    }
    names(compose) = 0:(length(compose)-1)
    return(compose)
}

#line.eq.2points: Gets x,y-coordinates of two points and returns coefficients of line passing through them
#INPUT:
# x1,y1,x2,y2 = (numeric) x,y-coordinates of points 1 and 2, respectively
#OUTPUT:
# Returns (numeric vector) with coefficients of line in increasing order. Line a0 + a1*x results in c(a0, a1)
#' @keywords internal
line.eq.2points = function(x1, y1, x2, y2){
    slope = (y2-y1)/(x2-x1)
    y.intercept = y1 - slope * x1
    line = c(y.intercept, slope)
    names(line) = 0:(length(line)-1)
    return(line)
}

#my.F: Gets coefficients of line and coefficients of 3rd order polynomial and returns function F(x)
#INPUT:
# line   = (numeric vector) contains coefficients of line in increasing order
#          Example: p(x) = a0 + a1*x is entered as c(a0, a1)
# g.poly = (numeric vector) contain coefficients of 3rd order polynomial in increasing order
#          Example: p(x) = a0 + a1*x + a2*x^2 + a3*x^3 is entered as c(a0, a1, a2, a3)
#OUTPUT:
# Returns (numeric vector) with coefficients of F in increasing order.
#' @keywords internal
my.F = function(line, g.poly){
    B = line[1] #y-intercept
    m = line[2] #slope

    poly1 = 1 / (m^2 + 1) * c(-B * m, 1)
    poly2 = (m / (m^2 + 1)) * g.poly
    FF = poly.sum(poly1, poly2)

    names(FF) = 0:(length(FF)-1)
    return(FF)
}

#my.distance.derivative: Gets the coefficients of line "f", polynomial "g" and function "F" and
#returns the coefficients of derivative of Euclidean distance^2 (squared).
#INPUT:
# f  = (numeric vector) contains coefficients of line in increasing order
#      Example: Line p(x) = a0 + a1*x is entered as c(a0, a1)
# g  = (numeric vector) contains coefficients of erd order polynomial in increasing order
# FF = (numeric vector) contains coefficients of F function in increasing order
#OUTPUT:
# Returns (numeric vector) with coefficients of derivative of Distance^2 in increasing order.
#' @keywords internal
my.distance.derivative = function(f, g, FF){
    FF.deriv = poly.derivative(FF) #Derivative of F function
    f.deriv = poly.derivative(f) #Derivative of line (f)
    g.deriv = poly.derivative(g) #Derivative of 3rd order polynomial (g)
    #Derivative of Euclidean distance squared (using chain rule)
    D.deriv = 2 * poly.product( poly.sum(FF, c(0, -1)),
                                poly.sum(FF.deriv, c(-1)) ) +
        2 * poly.product( poly.sum(poly.composition(f,FF), -g),
                          poly.sum(f.deriv * FF.deriv, -g.deriv ))
    return(D.deriv)
}

#ComputeRoot: Gets coefficients of polymonial p(x) and extreme points x1 and x2, and
#returns the root of polynomial p(x) in the interval (x1,x2)
#INPUT:
# p(x)   = (numeric vector) contains coefficients of polynomial in increasing order
#           Example: p(x) = a0 + a1*x + a2*x^2 + ... + an*x^n is entered as c(a0, a1, a2,..., an)
# x1,x2  = (numeric) extreme points of interval (x1, x2) where root lies
#OUTPUT:
# Returns (numeric) root of polynomial p(x) within interval (x1,x2)
#' @keywords internal
ComputeRoot = function(p, x1, x2){
    roots = polyroot(p)
    roots.real = round(Re(roots[round(Im(roots), 5) == 0]), 5) #Get real roots
    root.index = intersect(which(roots.real > x1),
                           which(roots.real < x2)) #Get roots within interval
    root = roots.real[root.index]
    return(root)
}

#ComputeDmax: Given the line f, the polynomial g and the extreme points x1 and x2, returns Dmax.
# Algorithm is described in section 3 of documentation paper (steps 5-8).
#INPUT:
# f  = (numeric vector) contains coefficients of line in increasing order
#      Example: f(x) = a0 + a1*x is entered as c(a0, a1)
# g  = (numeric vector) contains coefficients of 3rd order polynomial in increasing order
#      Example: g(x) = a0 + a1*x + a2*x^2 + a3*x^3 is entered as c(a0, a1, a2, a3)
# x1 = (numeric) left extreme of interval where Dmax lies
# x2 = (numeric) right extreme of interval where Dmax lies
#OUTPUT:
# Returns (numeric) Dmax, which is the x-coord of the point on the polynomial curve where the maximum
# perpendicular distance between the line and the polynomial is reached in the interval (x1, x2).
#' @keywords internal
ComputeDmax = function(f, g, x1, x2){
    ##### Compute F function
    FF = my.F(line=f, g.poly=g)

    ##### Compute derivative of Euclidean distance^2
    D.deriv = my.distance.derivative(f, g, FF)

    ##### Determine D.max: Find root of derivative of distance in interval (x1,x2) that maximizes distance
    D.roots = ComputeRoot(D.deriv, x1, x2) #Get roots in interval (x1,x2)
    #Evaluate distance at roots
    D.roots.eval = sapply(D.roots,
                          function(i){
                              distance.2points(x1 = poly.evaluate(FF, i),
                                               y1 = poly.evaluate(f, poly.evaluate(FF, i)),
                                               x2 = i,
                                               y2 = poly.evaluate(g, i))
                          })
    #Get root that gives maximum distance
    D.roots.max = max(D.roots.eval)
    D.max = D.roots[which(D.roots.eval == D.roots.max)]
    return(D.max)
}
