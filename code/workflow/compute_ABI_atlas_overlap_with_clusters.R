#!/usr/bin/Rscript

# NOTES:
#
# Compute overlaps (both with sensitivity and ppv, along with Dice (F1 score) and Jaccard metric) between cluster labels for a given number of clusters and actual atlas
# Allen atlas: assess overlap at each node of thalamus tree, after filtering structures not in atlas

# Load libraries
library(tidyverse)
library(glue)
library(RMINC)
library(grid)
library(MRIcrotome)
library(data.tree)
library(patchwork)

source("code/functions/functions.R")

###########################
# Arguments
###########################

print(glue("[{Sys.time()}] Input arguments"))

# Args
args <- commandArgs(trailingOnly = T)
indir <- args[1]
mode <- args[2]
n_clusters <- as.integer(args[3])

# Testing arguments
# indir <- "data/outputs/isocortex_thalamus_left"
# mode <- "residual"
# n_clusters <- 8

###

###########################
# Define atlas paths and work on plot directory
###########################

thalamus_atlas_file <- "data/resources/AMBA_relabeled_backsampled_50um_split_THonly_left.mnc"
thalamus_tree_file <- "data/resources/AMBA_relabeled_backsampled_25um_split_hanat.RData"
outdir <- glue("{indir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/ABI_atlas_overlaps")

print(glue("[{Sys.time()}] Output directory is: {outdir}"))
# Cleanup directory if it exists
if (dir.exists(outdir)) {
  
  # Remove directory if it exists
  print(glue("[{Sys.time()}] Output directory already exists, therefore removing it and starting fresh..."))
  unlink(outdir, recursive = T)
  
} else {
  
  # Create output plot directory
  print(glue("[{Sys.time()}] Creating output directory"))
  dir.create(outdir, showWarnings = T, recursive = T)
  
}

###########################
# Read data
###########################

print(glue("[{Sys.time()}] Reading data"))

thalamus_atlas_vol <- mincGetVolume(thalamus_atlas_file)
cluster_label_vol <- mincGetVolume(glue("{indir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/cluster_labels.mnc"))
thalamus_mask_indices <- which(thalamus_atlas_vol > 0.5)
full_atlas_vol <- mincGetVolume("data/resources/AMBA_relabeled_backsampled_50um_split.mnc")

thalamus_atlas_vec <- as.integer(round(thalamus_atlas_vol))[thalamus_atlas_vol > 0.5]
cluster_label_vec <- as.integer(round(cluster_label_vol))[thalamus_atlas_vol > 0.5]

# Load tree
print(glue("[{Sys.time()}] Loading ABI tree definitions"))
load(thalamus_tree_file)
thalamus_tree <- tree_original$`Basic cell groups and regions`$`Brain stem`$Interbrain$Thalamus

###########################
# Preprocess data
#
# Remove bad labels around edges that come up due to resampling 25um high resolution atlas to 50um
#
###########################

print(glue("[{Sys.time()}] Preprocessing atlas data"))

# Thalamus labels
thalamus_atlas_labels <- sort(unique(thalamus_atlas_vec))
thalamus_tree_labels <- thalamus_tree$Get("id")

# Which structures (defined in tree) are in the atlas?
names(thalamus_tree_labels[which((thalamus_tree_labels %in% thalamus_atlas_labels))])

# Resampling of thalamus mask at 50um causes a few voxels of neighbouring structures to come into thalamus atlas -- remove these 
# Note: these voxels are also in mask_vol_thalamus, therefore the cluster labels should also be fixed
# Which labels not in thalamus tree?
non_thalamus_labels <- thalamus_atlas_labels[which(!c(thalamus_atlas_labels %in% thalamus_tree_labels))]
non_thalamus_structures <- sort(tree_original$Get("id")[which(tree_original$Get("id") %in% non_thalamus_labels)])

# Remove these voxels from both thalamus atlas and cluster labels
length(which(thalamus_atlas_vec %in% non_thalamus_labels))
thalamus_atlas_vec[which(thalamus_atlas_vec %in% non_thalamus_labels)] <- 0
cluster_label_vec[which(thalamus_atlas_vec %in% non_thalamus_labels)] <- 0

# Get labels vector again without non-thalamus labels
thalamus_atlas_labels <- sort(unique(thalamus_atlas_vec))[-1]

print(glue("[{Sys.time()}] Rewriting processed (fixed) atlases"))

# Rewrite out atlas
print(glue("[{Sys.time()}] Writing thalamus atlas"))
ovol <- thalamus_atlas_vol
ovol[] <- 0
ovol[thalamus_atlas_vol > 0.5] <- thalamus_atlas_vec
mincWriteVolume(ovol, glue("{outdir}/thalamus_atlas_volume_fixed.mnc"))

print(glue("[{Sys.time()}] Writing SC-based cluster labels"))
ovol <- thalamus_atlas_vol
ovol[] <- 0
ovol[thalamus_atlas_vol > 0.5] <- cluster_label_vec
mincWriteVolume(ovol, glue("{outdir}/cluster_labels_fixed.mnc"))

# Write out flipped atlas and merged volume
print(glue("[{Sys.time()}] Choosing side"))
ovol <- mvol <- full_atlas_vol
if (mean(thalamus_atlas_labels) > 20000) {
  # If left side
  print(glue("> Clusters on left, thalamus atlas on right"))
  ovol[!(as.integer(round(ovol[])) %in% (thalamus_atlas_labels - 20000))] <- 0
  mvol[!(as.integer(round(mvol[])) %in% (thalamus_atlas_labels - 20000))] <- 0
} else {
  # Else labels are on right side
  print(glue("> Clusters on right, thalamus atlas on left"))
  ovol[!(as.integer(round(ovol[])) %in% (thalamus_atlas_labels + 20000))] <- 0
  mvol[!(as.integer(round(mvol[])) %in% (thalamus_atlas_labels + 20000))] <- 0
}

print(glue("[{Sys.time()}] Writing thalamus atlas, opposite to clusters"))
mincWriteVolume(ovol, glue("{outdir}/thalamus_atlas_volume_opposite_fixed.mnc"))

# Write out merged atlas
print(glue("[{Sys.time()}] Writing merged cluster labels and thalamus atlas"))
mvol[cluster_label_vol > 0.5] <- cluster_label_vol[cluster_label_vol > 0.5]
mincWriteVolume(mvol, glue("{outdir}/cluster_labels_fixed_with_thalamus_atlas_volume_fixed.mnc"))

###########################
# Define functions
###########################

print(glue("[{Sys.time()}] Defining functions"))

# Get atlas label info
get_overlap_statistics <- function(cluster_mask, atlas_mask) {
  tp <- length(which(cluster_mask==TRUE & atlas_mask==TRUE))
  tn <- length(which(cluster_mask==FALSE & atlas_mask==FALSE))
  fp <- length(which(cluster_mask==TRUE & atlas_mask==FALSE))
  fn <- length(which(cluster_mask==FALSE & atlas_mask==TRUE))
  
  sensitivity <- tp/(tp+fn)
  specificity <- tn/(tn+fp)
  ppv <- tp/(tp+fp)
  npv <- tn/(tn+fn)
  accuracy <- (tp + tn)/(tp + tn + fp + fn)
  f1 <- 2*tp/(2*tp + fp + fn)
  jaccard <- tp/(tp + fn + fp)
  
  return(c(sensitivity=sensitivity,
           specificity=specificity,
           ppv=ppv,
           npv=npv,
           accuracy=accuracy,
           f1=f1,
           jaccard=jaccard))
}

compare_atlas_labels <- function(thalamus_atlas_vec, cluster_label_vec, cluster_label, atlas_label) {
  structure_name <- names(which(thalamus_tree$Get("id")==atlas_label))[1]
  structure_node <- FindNode(thalamus_tree, structure_name)
  if (!structure_node$is_unsplit) {
    structure_node <- structure_node$parent
  }
  structure_name <- structure_node$name
  structure_acro <- structure_node$acronym
  structure_level <- structure_node$level
  structure_labels <- sort(unique(structure_node$Get("id")))
  structure_order <- structure_node$graph_order
  structure_color <- structure_node$color_hex_triplet
  cl_mask <- (cluster_label_vec %in% cluster_label)
  at_mask <- (thalamus_atlas_vec %in% structure_labels)
  
  olstats <- get_overlap_statistics(cluster_mask = cl_mask, atlas_mask = at_mask)
  outdf <- tibble(cluster_label=cluster_label, 
                  structure_level=structure_level, 
                  structure_label=atlas_label, 
                  structure_name=structure_name,
                  structure_acronym=structure_acro, 
                  structure_order=structure_order,
                  structure_color=structure_color,
                  sensitivity=olstats["sensitivity"],
                  specificity=olstats["specificity"],
                  ppv=olstats["ppv"],
                  npv=olstats["npv"],
                  accuracy=olstats["accuracy"],
                  f1=olstats["f1"],
                  jaccard=olstats["jaccard"])
  
  return(outdf)
}

###########################
# Reread data
###########################

print(glue("[{Sys.time()}] Computing overlaps between each cluster and each atlas label"))

# Only compute overlaps between nonrepeated labels
# Labels on left side (leaf, split) are essentially the same as their parent structure (unsplit)
# Therefore no need to compute these leaf labels (id > 20000)
atlas_labels <- sort(unique(thalamus_tree$Get("id")))
r_strucs <- atlas_labels[which(atlas_labels < 20000)]
l_strucs <- atlas_labels[which(atlas_labels > 20000)]-20000
setdiff(r_strucs, l_strucs) # Upstream (nonleaf) structures 
setdiff(l_strucs, r_strucs) # Verify that we're not missing any labels
atlas_labels <- r_strucs

# Loop over cluster labels
# Loop over atlas labels
if (exists("df_overlaps")) {
  rm(df_overlaps)
}
for (cluster_label in 1:n_clusters) {
  print(glue("Working on cluster label = {cluster_label} (of {n_clusters})"))
  for (atlas_label in atlas_labels) {
    odf <- compare_atlas_labels(thalamus_atlas_vec = thalamus_atlas_vec, 
                                cluster_label_vec = cluster_label_vec, 
                                cluster_label = cluster_label, 
                                atlas_label = atlas_label)
    if (exists("df_overlaps")) {
      if (nrow(merge(odf,df_overlaps)) == 0) {
        df_overlaps <- rbind(df_overlaps, odf)
      }
    } else {
      df_overlaps <- odf
    }
  }
}
df_overlaps <- df_overlaps %>% 
  mutate(n_clusters=n_clusters) %>%
  select(n_clusters, everything()) %>%
  arrange(cluster_label, structure_level, structure_label)

###########################
# Write out overlap and symlink
###########################

print(glue("[{Sys.time()}] Writing out overlap data and symlinking cluster data directory"))

# Output overlap data
write_csv(x=df_overlaps, path=glue("{outdir}/overlap_data.csv"))

###################################

print(glue("[{Sys.time()}] Done!"))
