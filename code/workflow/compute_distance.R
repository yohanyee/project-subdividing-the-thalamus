#!/usr/bin/Rscript

# Load libraries
library(tidyverse)
library(glue)
library(RMINC)
library(RcppCNPy)

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

###########################
# Read data
###########################

print(glue("[{Sys.time()}] Reading data"))

mask_vol_thalamus <- mincGetVolume(mask_file_thalamus)
mask_vol_cortex <- mincGetVolume(mask_file_cortex)

p_thal <- length(mask_vol_thalamus[mask_vol_thalamus > 0.5])
p_cortex <- length(mask_vol_cortex[mask_vol_cortex > 0.5])

###########################
# Get world coords of thalamic voxels
###########################

print(glue("[{Sys.time()}] Getting world coordinates [x,y,z] of thalamic voxels"))

thalamus_mask_indices <- which(mask_vol_thalamus > 0.5)

# Get world coords of thalamic voxels
# mincVectorToVoxelCoordinates returns z, y, x when seen in Display
voxel_coords_thalamus <- matrix(mincVectorToVoxelCoordinates(mask_file_thalamus, thalamus_mask_indices), ncol=3)

# Get world coordinate
# Since mincConvertVoxelToWorld takes in data as z, y, x, and the input is given in this order, the output is in the correct order as seen in Display
# Takes about two minutes
world_coords_thalamus <- t(mincConvertVoxelMatrix(filename = mask_file_thalamus, voxel_matrix = t(voxel_coords_thalamus)[1:3,]))

###########################
# Get world coords of thalamic voxels
###########################

print(glue("[{Sys.time()}] Getting world coordinates [x,y,z] of cortical voxels"))

cortex_mask_indices <- which(mask_vol_cortex > 0.5)

# Get world coords of cortical voxels
# mincVectorToVoxelCoordinates returns z, y, x when seen in Display
voxel_coords_cortex <- matrix(mincVectorToVoxelCoordinates(mask_file_cortex, cortex_mask_indices), ncol=3)

# Get world coordinate
# Since mincConvertVoxelToWorld takes in data as z, y, x, and the input is given in this order, the output is in the correct order as seen in Display
# Takes a few seconds
world_coords_cortex <- t(mincConvertVoxelMatrix(filename = mask_file_cortex, voxel_matrix = t(voxel_coords_cortex)[1:3,]))

###########################
# Compute distance
###########################

print(glue("[{Sys.time()}] Computing pairwise distance"))

distmat <- matrix(NA, nrow=p_thal, ncol=p_cortex)
prog <- txtProgressBar(max=p_thal, style=3)

# Loop over each dimension
# Should take around 2-3 minutes
for (i in 1:p_thal) {
  for (j in 1:p_cortex) {
    distmat[i,j] <- sqrt((world_coords_thalamus[i,1] - world_coords_cortex[j,1])^2 + (world_coords_thalamus[i,2] - world_coords_cortex[j,2])^2 + (world_coords_thalamus[i,3] - world_coords_cortex[j,3])^2)
  }
  setTxtProgressBar(prog, i)
}
close(prog)

dim(distmat)
distmat[1:5, 1:5]

print(glue("[{Sys.time()}] Minimum thalamocortical distance is {min(distmat)} mm"))
print(glue("[{Sys.time()}] Minimum thalamocortical distance is {max(distmat)} mm"))

# Write out distance file

print(glue("[{Sys.time()}] Writing out data"))

npySave(filename = glue("{outdir}/euclidean_distance.npy"), object = distmat)

print(glue("[{Sys.time()}] Done!"))

# Test that the distance map is OK
k <- 15000
ovol <- mask_vol_cortex
ovol[] <- 0
ovol[mask_vol_cortex > 0.5] <- distmat[k,]
min(distmat[k,])
max(distmat[k,])
voxel_coords_test <- mincVectorToVoxelCoordinates(mask_file_thalamus, thalamus_mask_indices[k])
world_coords_test <- mincConvertVoxelToWorld(mask_file_thalamus, voxel_coords_test[1], voxel_coords_test[2], voxel_coords_test[3])
mincWriteVolume(ovol, glue("tmp/distance_test_at_world_{world_coords_test[1]}_{world_coords_test[2]}_{world_coords_test[3]}.mnc"))
