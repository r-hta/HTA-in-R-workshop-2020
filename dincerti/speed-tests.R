# Dependencies and settings ----------------------------------------------------
rm(list = ls())
source("speed-tests-fun.R")
library("hesim")
library("heemod")
library("data.table")

# Run simulations
heemod_1000 <- run_heemod(1000)
hesim_1000 <- run_hesim_indiv(1000)
hesim_1000 <- run_hesim_cohort(1000)