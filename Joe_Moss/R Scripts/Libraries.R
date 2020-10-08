################################################################################
## Project: Introduction to PLS
## Script name: Libraries.R
## Script purpose: This script loads all the required R libraries
## Date: September 2020
## Author: Joe Moss (joe.moss@york.ac.uk)
## Organisation: York Health Economics Consortium (YHEC)
################################################################################


# List of required packages ==========================================================================================================

Packages <- c("doSNOW", "foreach", "knitr", "scales", "devtools", "reshape2", "plotly", "tidyverse")

#lapply(Packages, install.packages, character.only = TRUE, dependicies = TRUE)
#devtools::install_github("DARTH-git/dampack")

lapply(c(Packages,"dampack"), library, character.only = TRUE)
