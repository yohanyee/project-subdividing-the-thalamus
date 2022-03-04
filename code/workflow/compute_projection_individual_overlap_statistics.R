library(tidyverse)
library(glue)
library(RMINC)
library(RcppCNPy)
library(data.table)

# Testing arguments
outdir <- "data/outputs/isocortex_thalamus_left"
mode <- "residual"
n_clusters <- 6

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
mode <- args[2]
n_clusters <- as.integer(args[3])

#
sc_thres <- 0.1
tracer_thres <- 0.1
seed_thres <- 0.9

# SC file
sc_file <- glue("{outdir}/correlation_residual.npy")

# Masks
mask_file_thalamus <- glue("{outdir}/thalamus_mask_highres.mnc")
mask_file_cortex <- glue("{outdir}/cortex_mask_lowres.mnc")

# Read volumes
mask_vol_thalamus <- mincGetVolume(mask_file_thalamus)
mask_vol_cortex <- mincGetVolume(mask_file_cortex)

p_thal <- length(mask_vol_thalamus[mask_vol_thalamus > 0.5])
p_cortex <- length(mask_vol_cortex[mask_vol_cortex > 0.5])

# Read SC data
scdat <- npyLoad(sc_file)
dim(scdat)

# Get seed and target overlapping experiments
df_seed_all <- 1:n_clusters %>%
  map_dfr(function(k) {
    out <- read_csv(glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/projection_overlaps/{k}/seed_overlap_cluster_{k}.csv")) %>%
      mutate(set="seed_overlap",
             n_clusters=n_clusters,
             cluster=k) %>%
      select(set, n_clusters, cluster, everything())
    out
  })

df_tracer_all <- 1:n_clusters %>%
  map_dfr(function(k) {
    out <- read_csv(glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/projection_overlaps/{k}/tracer_overlap_cluster_{k}.csv")) %>%
      mutate(set="seed_overlap",
             n_clusters=n_clusters,
             cluster=k) %>%
      select(set, n_clusters, cluster, everything())
    out
  })

# Filter out experiments in which a seed or target does not fully overlap with the cluster
df_seed <- df_seed_all %>%
  filter(percent_of_seed_covered_by_roi >= 0.5)
table(df_seed$cluster) # Number of experiments with seed in each cluster

df_tracer <- df_tracer_all %>%
  filter(percent_of_roi_covered_by_tracer >= 0.5)
table(df_tracer$cluster) # Number of experiments with tracer projecting to cluster

# Define function
get_overlap_statistics <- function(sc, tracer) {
  tp <- length(which(sc==TRUE & tracer==TRUE))
  tn <- length(which(sc==FALSE & tracer==FALSE))
  fp <- length(which(sc==TRUE & tracer==FALSE))
  fn <- length(which(sc==FALSE & tracer==TRUE))
  
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

##########################
# For each seed (injection site) overlapping with thalamus experiment
seed_olstats_list <- list()
for (k in 1:nrow(df_seed)) {
  
  print(glue("Working on experiment {k} of {nrow(df_seed)}"))
  
  seed_row <- df_seed[k,]
  
  # Injection site
  seed_indices <- which(mincGetVolume(seed_row$seed_file)[mask_vol_thalamus > 0.5] > seed_thres)
  
  # Cortical map associated with injection site voxels
  sc_vec <- colMeans(scdat[seed_indices,])
  
  # Test variables
  sc_pos <- sc_vec >= sc_thres
  sc_neg <- sc_vec <= -sc_thres
  
  # Likefile
  likefile_200um <- "data/resources/average_template_200.mnc"
  
  # Get 200um tracer
  tracer_file_50 <- seed_row$tracer_file
  tracer_file_200 <- tempfile(pattern = "resampled_tracer-", tmpdir = "/dev/shm/", fileext = ".mnc")
  system(cmd <- glue("mincresample -2 -like {likefile_200um} -keep_real_range {tracer_file_50} {tracer_file_200}"))
  
  # Get tracer vector and cleanup
  tvec <- mincGetVolume(tracer_file_200)[mask_vol_cortex > 0.5]
  tpos <- tvec >= tracer_thres
  file.remove(tracer_file_200)
  
  # Get cortical tracer (ground truth) overlap with positive SC region (test)
  seed_olstats_pos <- get_overlap_statistics(sc_pos, tpos) %>% 
    t %>%
    as_tibble() %>%
    mutate(sc_direction="positive") %>%
    select(sc_direction, everything()) %>%
    cbind(seed_row, .)
  
  # Get cortical tracer (ground truth) overlap with negative SC region (test)
  seed_olstats_neg <- get_overlap_statistics(sc_neg, tpos) %>% 
    t %>%
    as_tibble() %>%
    mutate(sc_direction="negative") %>%
    select(sc_direction, everything()) %>%
    cbind(seed_row, .)
  
  # Bind information together
  seed_olstats <- rbind(seed_olstats_pos, seed_olstats_neg)
  seed_olstats_list[[k]] <- seed_olstats
}
seed_olstats <- rbindlist(seed_olstats_list)

##########################
# For each tracer overlapping with thalamus experiment
tracer_olstats_list <- list()
for (k in 1:nrow(df_tracer)) {
  
  print(glue("Working on experiment {k} of {nrow(df_tracer)}"))
  
  tracer_row <- df_tracer[k,]
  
  # Injection site
  tracer_indices <- which(mincGetVolume(tracer_row$tracer_file)[mask_vol_thalamus > 0.5] > tracer_thres)
  
  # Cortical map associated with injection site voxels
  sc_vec <- colMeans(scdat[tracer_indices,])
  
  # Test variables
  sc_pos <- sc_vec >= sc_thres
  sc_neg <- sc_vec <= -sc_thres
  
  # Likefile
  likefile_200um <- "data/resources/average_template_200.mnc"
  
  # Get 200um injection site
  seed_file_50 <- tracer_row$seed_file
  seed_file_200 <- tempfile(pattern = "resampled_seed-", tmpdir = "/dev/shm/", fileext = ".mnc")
  system(cmd <- glue("mincresample -2 -like {likefile_200um} -keep_real_range {seed_file_50} {seed_file_200}"))
  
  # Get tracer vector and cleanup
  svec <- mincGetVolume(seed_file_200)[mask_vol_cortex > 0.5]
  spos <- svec >= seed_thres
  file.remove(seed_file_200)
  
  # Get cortical injection site (ground truth) overlap with positive SC region (test)
  tracer_olstats_pos <- get_overlap_statistics(sc_pos, spos) %>% 
    t %>%
    as_tibble() %>%
    mutate(sc_direction="positive") %>%
    select(sc_direction, everything()) %>%
    cbind(tracer_row, .)
  
  # Get cortical injection site (ground truth) overlap with negative SC region (test)
  tracer_olstats_neg <- get_overlap_statistics(sc_neg, spos) %>% 
    t %>%
    as_tibble() %>%
    mutate(sc_direction="negative") %>%
    select(sc_direction, everything()) %>%
    cbind(tracer_row, .)
  
  # Bind information together
  tracer_olstats <- rbind(tracer_olstats_pos, tracer_olstats_neg)
  tracer_olstats_list[[k]] <- tracer_olstats
}
tracer_olstats <- rbindlist(tracer_olstats_list)

# Write out
write_csv(seed_olstats, glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/projection_overlaps/individual_seed_overlap_statistics.csv"))
write_csv(tracer_olstats, glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/projection_overlaps/individual_tracer_overlap_statistics.csv"))


############################
