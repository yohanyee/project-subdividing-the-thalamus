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
# decompdir <- "data/outputs/isocortex_thalamus_left/decomposition_z/kmeans/10_clusters"
# n_clusters <- 10

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
decompdir <- args[2]
n_clusters <- as.integer(args[3])

# Input files
mask_file_thalamus <- glue("{outdir}/thalamus_mask_highres.mnc")
mask_file_cortex <- glue("{outdir}/cortex_mask_lowres.mnc")
array_file_labels <- glue("{decompdir}/cluster_labels.npy")
array_file_centers <- glue("{decompdir}/cluster_centers.npy")


###########################
# Read data
###########################

print(glue("[{Sys.time()}] Reading data"))

mask_vol_thalamus <- mincGetVolume(mask_file_thalamus)
mask_vol_cortex <- mincGetVolume(mask_file_cortex)

p_thal <- length(mask_vol_thalamus[mask_vol_thalamus > 0.5])
p_cortex <- length(mask_vol_cortex[mask_vol_cortex > 0.5])

k_labels <- npyLoad(array_file_labels)
k_centers <- npyLoad(array_file_centers)

###########################
# Write out
###########################

# Link averages
file.symlink(normalizePath("data/resources/average_template_50.mnc"), glue("{decompdir}/average_template_50.mnc"))
file.symlink(normalizePath("data/resources/average_template_200.mnc"), glue("{decompdir}/average_template_200.mnc"))
file.symlink(normalizePath("data/resources/AMBA_relabeled_50um_right_side.mnc"), glue("{decompdir}/AMBA_relabeled_50um_right_side.mnc"))
file.symlink(normalizePath("data/resources/AMBA_relabeled_50um_left_side.mnc"), glue("{decompdir}/AMBA_relabeled_50um_left_side.mnc"))
file.symlink(normalizePath("data/resources/AMBA_relabeled_50um_right_side_THonly.mnc"), glue("{decompdir}/AMBA_relabeled_50um_right_side_THonly.mnc"))
file.symlink(normalizePath("data/resources/AMBA_relabeled_50um_left_side_THonly.mnc"), glue("{decompdir}/AMBA_relabeled_50um_left_side_THonly.mnc"))

# Write out cluster labels
ovol <- mask_vol_thalamus
ovol[] <- 0
ovol[mask_vol_thalamus > 0.5] <- k_labels
mincWriteVolume(ovol, glue("{decompdir}/cluster_labels.mnc"))

# Merge cluster labels with AMBA labels
system(cmd <- glue("mincmath -2 -add {decompdir}/cluster_labels.mnc data/resources/AMBA_relabeled_50um_right_side.mnc {decompdir}/cluster_labels_with_AMBA_relabeled_50um.mnc"))

# Merge cluster labels with AMBA labels, TH only
system(cmd <- glue("mincmath -2 -add {decompdir}/cluster_labels.mnc data/resources/AMBA_relabeled_50um_right_side_THonly.mnc {decompdir}/cluster_labels_with_AMBA_relabeled_50um_THonly.mnc"))

# Write out individual masks
for (i in 1:n_clusters) {
  ovol <- mask_vol_thalamus
  ovol[] <- 0
  outvec <- rep(0, length.out=p_thal)
  outvec[which(k_labels==i)] <- 1
  ovol[mask_vol_thalamus > 0.5] <- outvec
  mincWriteVolume(ovol, glue("{decompdir}/cluster_mask_{i}.mnc"))
}

# Write out cluster centers
for (i in 1:n_clusters) {
  ovol <- mask_vol_cortex
  ovol[] <- 0
  ovol[mask_vol_cortex > 0.5] <- k_centers[i,]
  mincWriteVolume(ovol, glue("{decompdir}/cluster_center_{i}.mnc"))
}

# Subdivide cortex based on center max
ovol <- mask_vol_cortex
ovol[] <- 0
ovol[mask_vol_cortex > 0.5] <- apply(k_centers, 2, which.max)
mincWriteVolume(ovol, glue("{decompdir}/cluster_center_max.mnc"))

# Subdivide cortex based on center max, after scaling
ovol <- mask_vol_cortex
ovol[] <- 0
ovol[mask_vol_cortex > 0.5] <- apply(t(scale(t(k_centers))), 2, which.max)
mincWriteVolume(ovol, glue("{decompdir}/cluster_center_scaled_max.mnc"))

# Subdivide cortex based on center min
ovol <- mask_vol_cortex
ovol[] <- 0
ovol[mask_vol_cortex > 0.5] <- apply(k_centers, 2, which.min)
mincWriteVolume(ovol, glue("{decompdir}/cluster_center_min.mnc"))

# Subdivide cortex based on center min, after scaling
ovol <- mask_vol_cortex
ovol[] <- 0
ovol[mask_vol_cortex > 0.5] <- apply(t(scale(t(k_centers))), 2, which.min)
mincWriteVolume(ovol, glue("{decompdir}/cluster_center_scaled_min.mnc"))
