library(gasExchangeR)
library(tidyverse)
library(janitor)
library(broom)
library(Deriv)
library(Ryacas)
library(readxl)
library(devtools)

# d1 sign change is most similar to Wisén and Wohlfart (2004)
# This function finds the first derivative of the VE/VCO2 vs. time
# and notes VT2 as the LAST time the first derivative crosses above zero
# The same procedure can be used to find VT1, but use VE/VO2 vs. time instead
