#!/usr/bin/Rscript

# Load libraries
library(tidyverse)
library(glue)
library(RMINC)
library(RcppCNPy)
library(normtools)

# NOTE: if modifying this script to accomodate other datasets, then also modify normalization model (if necessary) to work with covariates in other data

###########################
# Arguments
###########################

# Testing arguments
outdir <- "data/outputs/isocortex_thalamus_left"

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]

# Input dataframe
df <- read_csv("data/outputs/wt154.csv")
num_files <- nrow(df)

# Input files
mask_file_thalamus <- glue("{outdir}/thalamus_mask_highres.mnc")
mask_file_cortex <- glue("{outdir}/cortex_mask_lowres.mnc")

###########################
# Read data
###########################

mask_vol_thalamus <- mincGetVolume(mask_file_thalamus)
mask_vol_cortex <- mincGetVolume(mask_file_cortex)

p_thal <- length(mask_vol_thalamus[mask_vol_thalamus > 0.5])
p_cortex <- length(mask_vol_cortex[mask_vol_cortex > 0.5])

###########################
# THALAMUS
###########################

# Setup data matrix
data_matrix_thalamus <- matrix(NA, nrow = num_files, ncol=p_thal)
dim(data_matrix_thalamus)

# Read data matrix
# Takes about 1 minute
prog <- progress_estimated(num_files)
for (i in 1:num_files) {
  data_matrix_thalamus[i, ] <- mincGetVolume(df$resampled_determinants_CCFv3_abs_50[i])[mask_vol_thalamus > 0.5]
  print(prog$tick())
}

# Normalize data matrix
# Takes about 5 minutes on fievel when using left side thalamus mask at 50um
# - 7-10 minutes if full thalamus mask at 50um
# - may hang if using screen
Sys.time()
normalized_data_matrix_thalamus <- norm_by_model(data_matrix_thalamus, formula = ~ TwoLevel_Group + thalamus, df = df)
Sys.time()

# Save both matrices
npySave(filename = glue("{outdir}/data_matrix_thalamus.npy"), object=data_matrix_thalamus)
npySave(filename = glue("{outdir}/normalized_data_matrix_thalamus.npy"), object=normalized_data_matrix_thalamus)

###########################
# CORTEX
###########################

# Setup data matrix
data_matrix_cortex <- matrix(NA, nrow = num_files, ncol=p_cortex)
dim(data_matrix_cortex)

# Read data matrix
# Takes about 30 seconds
prog <- progress_estimated(num_files)
for (i in 1:num_files) {
  data_matrix_cortex[i, ] <- mincGetVolume(df$resampled_determinants_CCFv3_abs_200[i])[mask_vol_cortex > 0.5]
  print(prog$tick())
}

# Normalize data matrix
# Takes about 1 minute
Sys.time()
normalized_data_matrix_cortex <- norm_by_model(data_matrix_cortex, formula = ~ TwoLevel_Group + isocortex, df = df)
Sys.time()

# Save both matrices
npySave(filename = glue("{outdir}/data_matrix_cortex.npy"), object=data_matrix_cortex)
npySave(filename = glue("{outdir}/normalized_data_matrix_cortex.npy"), object=normalized_data_matrix_cortex)

print(glue("Done!"))