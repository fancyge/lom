##########################################################################################
####### Preparation for ACE modelling - filter species and pathway features
####### set cutoff for prevalence, abundance and variance
##########################################################################################
library(dplyr)
library(stringr)
library(tidyverse)
library(compositions) 

##################################################################
### filter species 
##################################################################
sp_abun <- read.csv("./inputs/comp_metaphlan.csv")
species_abun <- sp_abun[, -c(1:4)]

prevalence_cutoff <- 0.20   # ≥20% of samples must have non-zero
abundance_cutoff  <- 0.001  # ≥0.1% relative abundance
var_cutoff        <- 0.001  # variance threshold

n_samples <- nrow(species_abun)

prev_filter <- colSums(species_abun > 0) / n_samples >= prevalence_cutoff
abun_filter <- colMeans(species_abun) >= abundance_cutoff

log_abun <- log1p(species_abun)  # log(1+x) transform
var_filter <- apply(log_abun, 2, sd) >= var_cutoff

keep_species <- prev_filter & abun_filter & var_filter
filtered_abun <- cbind(sp_abun[, c(1:4)], species_abun[, keep_species, drop = FALSE])


##################################################################
### filter pathways
##################################################################
pwy_abun <- read.delim("./inputs/func_us_425.txt")
pwy_abun <- pwy_abun[,-1]
epsilon <- 1e-6
n_samples <- nrow(pwy_abun)

prev_filter <- colSums(pwy_abun > epsilon) / n_samples >= 0.20
var_cutoff  <- quantile(apply(pwy_abun, 2, sd), 0.25)   # drop lowest 25% SD
var_filter  <- apply(pwy_abun, 2, sd) >= var_cutoff
abun_filter <- colMeans(abs(pwy_abun)) >= 0.001

keep_pwy <- prev_filter & var_filter & abun_filter
filtered_pwy <- cbind(sp_abun[, c(1:4)], pwy_abun[, keep_pwy, drop = FALSE])

