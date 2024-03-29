% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/averages.R
\name{avg_exercise_test}
\alias{avg_exercise_test}
\title{Average gas exchange data from an exercise test}
\usage{
avg_exercise_test(
  .data,
  method = "breath",
  calc_type = "rolling",
  time_col = "time",
  roll_window = 15,
  bin_w = 15,
  align = "center",
  mos = "mean",
  roll_trim = 0,
  bin_trim = 0,
  cutoff = 0.04,
  fs = 1,
  order = 3
)
}
\arguments{
\item{.data}{Breath-by-breath gas exchange data.}

\item{method}{Choose between \code{breath} averages, \code{time} averages, or \code{digital} filtering.}

\item{calc_type}{Choose \code{rolling}, \code{bin}, or \code{bin-roll}.}

\item{time_col}{The name of the column with time.}

\item{roll_window}{How many seconds or breaths to include if rolling.}

\item{bin_w}{Bin size of breaths or time.}

\item{align}{If using a rolling method, how to align the rolling average. Default is \code{"center"} Other choices include \code{"left"}, and \code{"right"}.}

\item{mos}{'Measure of center'. Choices include \code{"mean"} (default) or \code{"median"}.}

\item{roll_trim}{Indicate if you want a trimmed mean. Roll_trim removes a number of data points equal to \code{roll_trim} This is used to emulate MCG's "mid-5-of-7" averaging method. \code{roll_trim} must be a positive, even integer.}

\item{bin_trim}{See \code{roll_trim}.}

\item{cutoff}{The cutoff frequency in Hz. Only used by digital filter.}

\item{fs}{The sampling frequency in Hz. Only used by digital filter.}

\item{order}{The Butterworth low-pass filter order. Only used by digital filter.}
}
\value{
A data frame.
}
\description{
This function averages first by either breath, time, or digital filtering
If averaging by breath or time averages, it can also perform rolling or bin averages. Furthermore, you can specify if you want a whole or trimmed mean.
}
\details{
If you combine rolling and bin averages with \code{calc_type = bin_roll} it is important to note how \code{roll_window} and \code{bin_w} interact. This first creates bin averages that are evenly divisible by the \code{roll_window}. For example, if your bin_w is 5, the function first computes bin averages every 5 breaths or seconds, depending on the \code{method} parameter. Then, if the roll_window was \code{15}, the rolling average would include 3 points in that average because 15 ÷ 5 = 3.

\code{roll_window} must be evenly divisible by \code{bin_w}
}
\examples{

# TODO write an example

# TODO force people to use this function twice if doing bin-roll or rolling-bin?

}
