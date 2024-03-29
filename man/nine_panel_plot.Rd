% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nine_panel_plot.R
\name{nine_panel_plot}
\alias{nine_panel_plot}
\title{Construct Wasserman's Nine Panel Plot}
\usage{
nine_panel_plot(
  .data,
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
  petco2 = "petco2"
)
}
\arguments{
\item{.data}{Gas exchange data frame or tibble}

\item{vt1_dat}{A data frame or tibble of values at VT1}

\item{vt2_dat}{A data frame or tibble of values at VT2/RC}

\item{vo2}{Name of the `VO2` variable}

\item{vco2}{Name of the `VCO2` variable}

\item{ve}{Name of the `pulmonary ventilation (VE)` variable}

\item{time}{Name of the `time` variable}

\item{hr}{Name of the `heart rate` variable}

\item{speed}{Name of the `speed` variable}

\item{grade}{Name of the `grade` variable}

\item{watts}{Name of the `watts` variable}

\item{vt}{Name of the `tidal volume (Vt)` variable}

\item{rer}{Name of the `respiratory exchange ratio (RER)` variable}

\item{peto2}{Name of the `end tidal oxygen (PetO2)` variable}

\item{petco2}{Name of the `end tidal carbon dioxide (PetCO2)` variable}
}
\value{
a ggplot object
}
\description{
This is a famous representation of a cardiopulmonary exercise test often used in clinical practice.
}
\examples{
# TODO write an example
}
