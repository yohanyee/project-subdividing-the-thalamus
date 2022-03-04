#!/usr/bin/Rscript

# NOTES:
#
# Analyze overlaps (both with sensitivity and ppv, along with combined metric) between cluster labels for a given number of clusters and actual atlas
# Allen atlas: assess overlap at each node of thalamus tree, after filtering structures not in atlas

# Load libraries
library(tidyverse)
library(glue)
library(RMINC)
library(data.tree)
library(data.table)
library(foreach)
library(doParallel)


###########################
# Arguments
###########################

print(glue("[{Sys.time()}] Input arguments"))

# Args
args <- commandArgs(trailingOnly = T)
indir <- args[1]
mode <- args[2]
n_clusters <- as.integer(args[3])
n_permutations <- as.integer(args[4])
n_threads <- as.integer(args[5])

# Testing arguments
# outdir <- "data/outputs/isocortex_thalamus_left"
# mode <- "residual"
# n_clusters <- 8
# n_permutations <- 5000
# n_threads <- 4

#
print(glue("> Using data from: {indir}"))
print(glue("> SC type: {mode}"))
print(glue("> Number of clusters: {n_clusters}"))
print(glue("> Number of permutations: {n_permutations}"))
print(glue("> Number of threads: {n_threads}"))

###########################
# Define atlas paths and output directory
###########################

thalamus_tree_file <- "data/resources/AMBA_relabeled_backsampled_25um_split_hanat.RData"
outdir <- glue("{indir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/ABI_atlas_overlaps")

print(glue("> Output directory is: {outdir}"))

###########################
# Read data
###########################

print(glue("[{Sys.time()}] Reading data"))

thalamus_atlas_vol <- mincGetVolume(glue("{outdir}/thalamus_atlas_volume_fixed.mnc"))
cluster_label_vol <- mincGetVolume(glue("{outdir}/cluster_labels_fixed.mnc"))

thalamus_atlas_vec <- as.integer(round(thalamus_atlas_vol))[thalamus_atlas_vol > 0.5]
cluster_label_vec <- as.integer(round(cluster_label_vol))[thalamus_atlas_vol > 0.5]

# Load tree
print(glue("[{Sys.time()}] Loading ABI tree definitions"))
load(thalamus_tree_file)
thalamus_tree <- tree_original$`Basic cell groups and regions`$`Brain stem`$Interbrain$Thalamus

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
n_voxels <- length(thalamus_atlas_vec)

# Loop over cluster labels
# Loop over atlas labels
if (exists("df_overlaps")) {
  rm(df_overlaps)
}

############################

print(glue("[{Sys.time()}] Setting up multithreading"))

cl <- makeCluster(n_threads, outfile="", rscript_args="--vanilla")
registerDoParallel(cl) 

export_objects <- c("thalamus_tree", "FindNode", "Get", "tibble",
                    "compare_atlas_labels", 
                    "thalamus_atlas_vec", "cluster_label_vec", "n_voxels",
                    "%>%", "mutate", "select", "arrange", "everything")

# For testing, uncomment below line
# n_permutations <- 100



#########################################
# Do the permutations
print(glue("[{Sys.time()}] Running loop!"))

prog <- txtProgressBar(max=n_permutations, style = 3)
out <- foreach(i=1:n_permutations, 
               .export = c(export_objects), 
               .verbose = TRUE, 
               .inorder = FALSE) %dopar% {
                 
                 # Permute labels               
                 cluster_label_vec_permuted <- cluster_label_vec[sample(1:n_voxels, size=n_voxels, replace=F)]
                 
                 # Compute permutation overlaps for each cluster
                 if (exists("df_overlaps")) {
                   rm(df_overlaps)
                 }
                 for (cluster_label in 1:n_clusters) {
                   
                   # Compute permutation overlaps for each label
                   for (atlas_label in atlas_labels) {
                     odf <- compare_atlas_labels(thalamus_atlas_vec = thalamus_atlas_vec, 
                                                 cluster_label_vec = cluster_label_vec_permuted, 
                                                 cluster_label = cluster_label, 
                                                 atlas_label = atlas_label)
                     
                     # Combine
                     if (exists("df_overlaps")) {
                       if (nrow(merge(odf,df_overlaps)) == 0) {
                         df_overlaps <- rbind(df_overlaps, odf)
                       }
                     } else {
                       df_overlaps <- odf
                     }
                     
                   }
                   
                 }
                 
                 # Put together permutation data for individual run
                 df_overlaps <- df_overlaps %>% 
                   mutate(n_clusters=n_clusters, 
                          perm=i) %>%
                   select(n_clusters, perm, everything()) %>%
                   arrange(cluster_label, structure_level, structure_label)
                 
                 # Progress and output
                 setTxtProgressBar(prog, i)
                 df_overlaps
               }
close(prog)

print(glue("[{Sys.time()}] Done loop. Stopping cluster"))

registerDoSEQ()
stopCluster(cl)

# Bind data together
print(glue("[{Sys.time()}] Binding data output"))
df_overlaps_perm <- rbindlist(out)

# Clean up 
print(glue("[{Sys.time()}] Cleaning up"))
rm(out)
gc()

# Output data
print(glue("[{Sys.time()}] Writing out data"))
write_csv(x = df_overlaps_perm, path = glue("{outdir}/overlaps_permuted.csv"))

###################################

print(glue("[{Sys.time()}] Done!"))
