#!/usr/bin/Rscript

# Load libraries
library(tidyverse)
library(glue)
library(RMINC)
library(RcppCNPy)
library(psych)

###########################
# Arguments
###########################

# Testing arguments
# outdir <- "data/outputs/isocortex_thalamus_left"

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]

# Input files
mask_file_thalamus <- glue("{outdir}/thalamus_mask_highres.mnc")
mask_file_cortex <- glue("{outdir}/cortex_mask_lowres.mnc")
data_matrix_file_thalamus <- glue("{outdir}/normalized_data_matrix_thalamus.npy")
data_matrix_file_cortex <- glue("{outdir}/normalized_data_matrix_cortex.npy")

###########################
# Read data
###########################

print(glue("[{Sys.time()}] Reading data"))

mask_vol_thalamus <- mincGetVolume(mask_file_thalamus)
mask_vol_cortex <- mincGetVolume(mask_file_cortex)

data_matrix_thalamus <- npyLoad(data_matrix_file_thalamus)
data_matrix_cortex <- npyLoad(data_matrix_file_cortex)

p_thal <- length(mask_vol_thalamus[mask_vol_thalamus > 0.5])
p_cortex <- length(mask_vol_cortex[mask_vol_cortex > 0.5])

###########################
# Compute correlation
###########################

dim(data_matrix_thalamus)
dim(data_matrix_cortex)



# Takes 2-3 minutes
print(glue("[{Sys.time()}] Computing correlation"))
cormat <- cor(data_matrix_thalamus, data_matrix_cortex)

# Takes about 30 seconds
print(glue("[{Sys.time()}] Computing z-transform"))
zmat <- fisherz(cormat)

# Examine a few things
dim(zmat)
zmat[1:5, 1:5]

print(glue("[{Sys.time()}] Minimum thalamocortical z-transformed correlation is {round(min(zmat), 5)}"))
print(glue("[{Sys.time()}] Minimum thalamocortical z-transformed correlation is {round(max(zmat), 5)}"))

# Write out SC (z) matrix

print(glue("[{Sys.time()}] Writing out data"))

npySave(filename = glue("{outdir}/correlation_z.npy"), object = zmat)

print(glue("[{Sys.time()}] Done!"))

# Test that the SC map is OK
for (k in c(1,2,3,10,100,200,300,1000, 15000, 15001, 15002, 25000, 40000, 50000, 50005, 50050, 50500, 55000, 75000)) {
  print(glue("Working on k = {k}"))
  ovol <- mask_vol_cortex
  ovol[] <- 0
  ovol[mask_vol_cortex > 0.5] <- zmat[k,]
  min(zmat[k,])
  max(zmat[k,])
  thalamus_mask_indices <- which(mask_vol_thalamus > 0.5)
  voxel_coords_test <- mincVectorToVoxelCoordinates(mask_file_thalamus, thalamus_mask_indices[k])
  world_coords_test <- mincConvertVoxelToWorld(mask_file_thalamus, voxel_coords_test[1], voxel_coords_test[2], voxel_coords_test[3])
  mincWriteVolume(ovol, glue("tmp/SC_z_test_{k}_at_world_{world_coords_test[1]}_{world_coords_test[2]}_{world_coords_test[3]}.mnc"))
}
