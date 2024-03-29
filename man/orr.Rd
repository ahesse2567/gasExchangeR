% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orr.R
\name{orr}
\alias{orr}
\title{Finding a breakpoint using Orr's 'bruteforce' algorithm.}
\usage{
orr(
  .data,
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
  front_trim_vt1 = 60,
  front_trim_vt2 = 60,
  pos_change = TRUE,
  pos_slope_after_bp = TRUE,
  ci = FALSE,
  conf_level = 0.95,
  plots = TRUE
)
}
\arguments{
\item{.data}{Gas exchange data.}

\item{.x}{The x-axis variable.}

\item{.y}{the y-axis variable.}

\item{bp}{Is this algorithm being used to find vt1 or vt2?}

\item{...}{Dot dot dot mostly allows this function to work properly if breakpoint() passes arguments that is not strictly needed by this function.}

\item{vo2}{The name of the vo2 column in \code{.data}}

\item{vco2}{The name of the vco2 column in \code{.data}}

\item{ve}{The name of the ve column in \code{.data}}

\item{time}{The name of the time column in \code{.data}}

\item{ordering}{Prior to fitting any functions, should the data be reordered by the x-axis variable or by time? Default is to use the current x-axis variable and use the time variable to break any ties.}

\item{alpha_linearity}{Significance value to determine if a piecewise model explains significantly reduces the residual sums of squares more than a simpler model.}

\item{front_trim_vt1}{How much data (in seconds) to remove from the beginning of the test prior to fitting any regressions. The original V-slope paper suggests 1 minute.}

\item{front_trim_vt2}{How much data (in seconds) to remove from the beginning of the test prior to fitting any regressions. The original V-slope paper suggests 1 minute.}

\item{pos_change}{Do you expect the change in slope to be positive (default) or negative? If a two-line regression explains significantly reduces the sum square error but the change in slope does not match the expected underlying physiology, the breakpoint will be classified as indeterminate.}

\item{pos_slope_after_bp}{Should the slope after the breakpoint be positive? Default is `TRUE`. This catches cases when the percent change in slope is positive, but the second slope is still negative. Change to `FALSE` when PetCO2 is the y-axis variable.}

\item{ci}{Should the output include confidence interval data? Default is `FALSE`.}

\item{conf_level}{Confidence level to use if calculating confidence intervals.}

\item{plots}{Should this function generate plots? Set to `FALSE` to save time.}
}
\value{
A list including slice of the original data frame at the threshold index with new columns `algorithm`, `determinant_bp`, `pct_slope_change`, `f_stat`, and `p_val_f.` The list also includes the fitted values, the left and right sides of the piecewise regression, and a simple linear regression.
}
\description{
Finding a breakpoint using Orr's 'bruteforce' algorithm.
}
\details{
According to the Orr et al. (1982, p. 1350) paper, the "anaerobic threshold is then reported as the first intersection point [of the lines] of the more appropriate model." In an email to Dr. Hughson, a weakness of this method is that the intersection point hardley ever lands on an existing data point. In some cases, the intersection of the best-fit model lies beyond the range of the x-axis. Instead of returning that value, we report the solution with the lowest pooled sums of squares that satisfies the following criteria: 1) the change in slope from the left to the right regression matches the anticipated change, usually positive; 2) the slope of the right regression matches the anticipated slope, usually positive; 3) the p-value of the F-test based on the extra sums of squares principle is statistically significant; and 4) the x-coordinate of the intersection point is within the range of the observed x-values. This differs slightly from the original method, but provides potentially plausible answers when the alternative is a negative VO2 or a VO2 well beyond VO2max.

The text of the original paper states that "Regression lines are calculated for all possible divisions of the data into two \emph{contiguous} groups." We interpret contiguous to mean "shares a boundary". That is, the left and right portions of the data each contain one identical data point at the division.
}
\examples{
# TODO write examples
}
\references{
Orr, G. W., Green, H. J., Hughson, R. L., & Bennett, G. W. (1982). A computer linear regression model to determine ventilatory anaerobic threshold. Journal of Applied Physiology Respiratory Environmental and Exercise Physiology, 52(5), 1349–1352. https://doi.org/10.1152/jappl.1982.52.5.1349
}
