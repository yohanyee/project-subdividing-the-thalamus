#!/usr/bin/Rscript

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
# n_clusters <- 6

###########################
# File paths
###########################

print(glue("[{Sys.time()}] Setting file paths"))

cluster_data_dir <- glue("{indir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/")
ABI_overlap_file <- glue("{cluster_data_dir}/ABI_atlas_overlaps/overlap_data.csv")
ABI_overlap_permutation_file <- glue("{cluster_data_dir}/ABI_atlas_overlaps/overlaps_permuted.csv")
outfile <- glue("{cluster_data_dir}/ABI_atlas_overlaps/overlap_data_with_pvalues.csv")

###########################
# Read data
###########################

print(glue("[{Sys.time()}] Setting file paths"))

df_overlaps <- read_csv(ABI_overlap_file)
df_overlaps_permuted <- read_csv(ABI_overlap_permutation_file)

###########################
# Define functions
###########################

print(glue("[{Sys.time()}] Defining function to compute p-value"))

get_p <- function(clust_label, struc_name, col_name) {
  
  # True metric value
  true_val <- df_overlaps %>%
    filter(cluster_label==clust_label,
           structure_name==struc_name) %>%
    pull(col_name)
  
  # Permutation distribution
  perm_dist <- df_overlaps_permuted %>%
    filter(cluster_label==clust_label,
           structure_name==struc_name) %>%
    pull(col_name)  

  # Histogram, uncomment for testing
  # Confirm that across different metrics, permutation distribution is normal
  # hist(perm_dist)
    
  # Compute p-value based on fitting normal distribution to permutation distribution
  # Note: use one sided p-value
  # Otherwise, when significance of cluster *not* overlapping structure is also measured
  # For two-sided pvalue, use: pval <- 2*min(iuc, 1-iuc)
  iuc <- pnorm(q = true_val, mean = mean(perm_dist), sd = sd(perm_dist))
  pval <- 1-iuc # One sided pvalue  
  
  # Return pvalue
  return(pval)
}


###########################
# Compute p-values
###########################

print(glue("[{Sys.time()}] Setting up dataframe"))

df_overlaps_withp <- df_overlaps %>%
  mutate(sensitivity_p=NA, 
         specificity_p=NA, 
         ppv_p=NA, 
         npv_p=NA,
         accuracy_p=NA,
         f1_p=NA,
         jaccard_p=NA)

print(glue("[{Sys.time()}] Computing p-values"))

# Takes about 2 minutes
prog <- txtProgressBar(max=nrow(df_overlaps_withp), style = 3)
for (i in 1:nrow(df_overlaps_withp)) {
  clust_label <- df_overlaps_withp$cluster_label[[i]]
  struc_name <- df_overlaps_withp$structure_name[[i]]
  
  df_overlaps_withp$sensitivity_p[[i]] <- get_p(clust_label, struc_name, "sensitivity")
  df_overlaps_withp$specificity_p[[i]] <- get_p(clust_label, struc_name, "specificity")
  df_overlaps_withp$ppv_p[[i]] <- get_p(clust_label, struc_name, "ppv")
  df_overlaps_withp$npv_p[[i]] <- get_p(clust_label, struc_name, "npv")
  df_overlaps_withp$accuracy_p[[i]] <- get_p(clust_label, struc_name, "accuracy")
  df_overlaps_withp$f1_p[[i]] <- get_p(clust_label, struc_name, "f1")
  df_overlaps_withp$jaccard_p[[i]] <- get_p(clust_label, struc_name, "jaccard")
  
  setTxtProgressBar(prog, i)
}
close(prog)

print(glue("[{Sys.time()}] Correcting p-values"))

df_overlaps_withp <- df_overlaps_withp %>%
  mutate(sensitivity_pfdr=p.adjust(sensitivity_p, method="fdr"), 
         specificity_pfdr=p.adjust(specificity_p, method="fdr"), 
         ppv_pfdr=p.adjust(ppv_p, method="fdr"), 
         npv_pfdr=p.adjust(npv_p, method="fdr"),
         accuracy_pfdr=p.adjust(accuracy_p, method="fdr"),
         f1_pfdr=p.adjust(f1_p, method="fdr"),
         jaccard_pfdr=p.adjust(jaccard_p, method="fdr"),
         sensitivity_pbonferroni=p.adjust(sensitivity_p, method="bonferroni"), 
         specificity_pbonferroni=p.adjust(specificity_p, method="bonferroni"), 
         ppv_pbonferroni=p.adjust(ppv_p, method="bonferroni"), 
         npv_pbonferroni=p.adjust(npv_p, method="bonferroni"),
         accuracy_pbonferroni=p.adjust(accuracy_p, method="bonferroni"),
         f1_pbonferroni=p.adjust(f1_p, method="bonferroni"),
         jaccard_pbonferroni=p.adjust(jaccard_p, method="bonferroni")
  )

###########################
# Write out results
###########################

print(glue("[{Sys.time()}] Writing dataframe"))
write_csv(df_overlaps_withp, path = outfile)

###########################
# Done!

print(glue("[{Sys.time()}] Done!"))
