library(tidyverse)
library(glue)
library(RMINC)
library(data.table)
library(foreach)
library(doParallel)

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
n_perms <- 5000
n_threads <- 20

# Masks
mask_file_thalamus <- glue("{outdir}/thalamus_mask_highres.mnc")
mask_file_cortex <- glue("{outdir}/cortex_mask_lowres.mnc")

# Read volumes
mask_vol_thalamus <- mincGetVolume(mask_file_thalamus)
mask_vol_cortex <- mincGetVolume(mask_file_cortex)

p_thal <- length(mask_vol_thalamus[mask_vol_thalamus > 0.5])
p_cortex <- length(mask_vol_cortex[mask_vol_cortex > 0.5])

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

# For each cluster
cluster_olstats_list <- list()
for (k in 1:n_clusters) {
  
  # Input files
  sc_file <- glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/cluster_center_{k}.mnc")
  sc_vol <- mincGetVolume(sc_file)
  sc_vec <- sc_vol[mask_vol_cortex > 0.5]
  
  # Test variables
  sc_pos <- sc_vec >= sc_thres
  sc_neg <- sc_vec <= -sc_thres
  
  # Likefile
  likefile_200um <- "data/resources/average_template_200.mnc"
  
  # Create tracer (seed overlap) merged for each cluster
  #merged_tracer_file <- tempfile(pattern = "merged_tracer-", tmpdir = "/dev/shm/", fileext = ".mnc")
  seed_overlap_merged_tracer_file_50um <- glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/projection_overlaps/{k}/seed_overlap_merged_max_tracer.mnc")
  seed_overlap_merged_seed_file_50um <- glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/projection_overlaps/{k}/seed_overlap_merged_max_seed.mnc")
  seed_overlap_merged_tracer_file_200um <- glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/projection_overlaps/{k}/seed_overlap_merged_max_tracer_200um.mnc")
  seed_overlap_merged_seed_file_200um <- glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/projection_overlaps/{k}/seed_overlap_merged_max_seed_200um.mnc")

  seed_overlap_tracer_files <- df_seed %>%
    filter(cluster==k) %>%
    pull(tracer_file)
  seed_overlap_seed_files <- df_seed %>%
    filter(cluster==k) %>%
    pull(seed_file)
  system(cmd <- glue("mincmath -2 -max {paste(seed_overlap_tracer_files, collapse = ' ')} {seed_overlap_merged_tracer_file_50um} -clobber"))
  system(cmd <- glue("mincmath -2 -max {paste(seed_overlap_seed_files, collapse = ' ')} {seed_overlap_merged_seed_file_50um} -clobber"))
  system(cmd <- glue("mincresample -2 -like {likefile_200um} -keep_real_range {seed_overlap_merged_tracer_file_50um} {seed_overlap_merged_tracer_file_200um} -clobber"))
  system(cmd <- glue("mincresample -2 -like {likefile_200um} -keep_real_range {seed_overlap_merged_seed_file_50um} {seed_overlap_merged_seed_file_200um} -clobber"))
  
  # Create tracer (tracer overlap) merged for each cluster
  tracer_overlap_merged_tracer_file_50um <- glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/projection_overlaps/{k}/tracer_overlap_merged_max_tracer.mnc")
  tracer_overlap_merged_seed_file_50um <- glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/projection_overlaps/{k}/tracer_overlap_merged_max_seed.mnc")
  tracer_overlap_merged_tracer_file_200um <- glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/projection_overlaps/{k}/tracer_overlap_merged_max_tracer_200um.mnc")
  tracer_overlap_merged_seed_file_200um <- glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/projection_overlaps/{k}/tracer_overlap_merged_max_seed_200um.mnc")
  
  tracer_overlap_tracer_files <- df_tracer %>%
    filter(cluster==k) %>%
    pull(tracer_file)
  tracer_overlap_seed_files <- df_tracer %>%
    filter(cluster==k) %>%
    pull(seed_file)
  system(cmd <- glue("mincmath -2 -max {paste(tracer_overlap_tracer_files, collapse = ' ')} {tracer_overlap_merged_tracer_file_50um} -clobber"))
  system(cmd <- glue("mincmath -2 -max {paste(tracer_overlap_seed_files, collapse = ' ')} {tracer_overlap_merged_seed_file_50um} -clobber"))
  system(cmd <- glue("mincresample -2 -like {likefile_200um} -keep_real_range {tracer_overlap_merged_tracer_file_50um} {tracer_overlap_merged_tracer_file_200um} -clobber"))
  system(cmd <- glue("mincresample -2 -like {likefile_200um} -keep_real_range {tracer_overlap_merged_seed_file_50um} {tracer_overlap_merged_seed_file_200um} -clobber"))
  
  # Compute overlap statistics for seed-cluster overlap (tracer in cortex vs SC in cortex)
  seed_overlap_vol <- mincGetVolume(seed_overlap_merged_tracer_file_200um)
  seed_overlap_vec <- seed_overlap_vol[mask_vol_cortex > 0.5]
  seed_overlap_pos <- seed_overlap_vec >= tracer_thres
  
  seed_olstats_pos <- get_overlap_statistics(sc_pos, seed_overlap_pos) %>% 
    t %>%
    as_tibble() %>%
    mutate(comparison="seed_overlap",
           sc_direction="positive",
           n_clusters=n_clusters,
           cluster=k) %>%
    select(comparison, n_clusters, cluster, sc_direction, everything())
  seed_olstats_neg <- get_overlap_statistics(sc_neg, seed_overlap_pos) %>% 
    t %>%
    as_tibble() %>%
    mutate(comparison="seed_overlap",
           sc_direction="negative",
           n_clusters=n_clusters,
           cluster=k) %>%
    select(comparison, n_clusters, cluster, sc_direction, everything())
  
  # Compute overlap statistics for tracer-cluster overlap (seed in cortex vs SC in cortex)
  tracer_overlap_vol <- mincGetVolume(tracer_overlap_merged_seed_file_200um)
  tracer_overlap_vec <- tracer_overlap_vol[mask_vol_cortex > 0.5]
  tracer_overlap_pos <- tracer_overlap_vec >= seed_thres
  
  tracer_olstats_pos <- get_overlap_statistics(sc_pos, tracer_overlap_pos) %>% 
    t %>%
    as_tibble() %>%
    mutate(comparison="tracer_overlap",
           sc_direction="positive",
           n_clusters=n_clusters,
           cluster=k) %>%
    select(comparison, n_clusters, cluster, sc_direction, everything())
  tracer_olstats_neg <- get_overlap_statistics(sc_neg, tracer_overlap_pos) %>% 
    t %>%
    as_tibble() %>%
    mutate(comparison="tracer_overlap",
           sc_direction="negative",
           n_clusters=n_clusters,
           cluster=k) %>%
    select(comparison, n_clusters, cluster, sc_direction, everything())
  
  # Bind information together
  cluster_olstats <- rbind(seed_olstats_pos, seed_olstats_neg, tracer_olstats_pos, tracer_olstats_neg)
  cluster_olstats_list[[k]] <- cluster_olstats
}
olstats <- rbindlist(cluster_olstats_list)

# Write out
write_csv(olstats, glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/projection_overlaps/merged_overlap_statistics.csv"))


############################

print(glue("[{Sys.time()}] Setting up multithreading"))
start_time <- Sys.time()

cl <- makeCluster(n_threads, outfile="", rscript_args="--vanilla")
registerDoParallel(cl) 

export_objects <- c("n_perms", "n_clusters", "outdir", "mode", "sc_thres", "seed_thres", "tracer_thres", "mask_vol_cortex",
                    "df_seed", "df_tracer", "mincGetVolume", "get_overlap_statistics", "glue", "rbindlist",
                    "%>%", "as_tibble", "mutate", "select", "filter", "nrow", "sample_n", "pull", "everything",
                    "start_time")

# For testing, uncomment below line
# n_perms <- 100

print(glue("[{Sys.time()}] Running loop!"))

# Repeat, generating permutation distribution
# For each permutation
# perm_list <- list()
#prog <- txtProgressBar(max=n_perms, style = 3)
perm_list <- 
  foreach(i=1:n_perms, 
          .export = c(export_objects), 
          .verbose = TRUE, 
          .inorder = FALSE) %dopar% {
            
            # For each cluster
            cluster_olperm_list <- list()
            for (k in 1:n_clusters) {
              
              # Input files
              print(glue("[{start_time} -> {Sys.time()} | {i}/{n_perms} | {k}/{n_clusters}] Read SC map"))
              sc_file <- glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/cluster_center_{k}.mnc")
              sc_vol <- mincGetVolume(sc_file)
              sc_vec <- sc_vol[mask_vol_cortex > 0.5]
              
              # Test variables
              sc_pos <- sc_vec >= sc_thres
              sc_neg <- sc_vec <= -sc_thres
              
              # Likefile
              likefile_200um <- "data/resources/average_template_200.mnc"
              
              # Create tracer (seed overlap) merged for each cluster
              print(glue("[{start_time} -> {Sys.time()} | {i}/{n_perms} | {k}/{n_clusters}] seed_overlap: Create merged tracer map"))
              #merged_tracer_file <- tempfile(pattern = "merged_tracer-", tmpdir = "/dev/shm/", fileext = ".mnc")
              seed_overlap_merged_tracer_file_50um <- tempfile(pattern = "merged_tracer-", tmpdir = "/dev/shm/", fileext = "_50um.mnc")
              seed_overlap_merged_tracer_file_200um <- gsub('_50um', '_200um', seed_overlap_merged_tracer_file_50um)
              
              n_seed <- df_seed %>%
                filter(cluster==k) %>%
                nrow()
              seed_overlap_tracer_files <- df_seed %>%
                sample_n(size = n_seed, replace = F) %>%
                pull(tracer_file)
              seed_overlap_seed_files <- df_seed %>%
                sample_n(size = n_seed, replace = F) %>%
                pull(seed_file)
              system(cmd <- glue("mincmath -2 -max {paste(seed_overlap_tracer_files, collapse = ' ')} {seed_overlap_merged_tracer_file_50um} -clobber"))
              system(cmd <- glue("mincresample -2 -like {likefile_200um} -keep_real_range {seed_overlap_merged_tracer_file_50um} {seed_overlap_merged_tracer_file_200um}"))
              
              Sys.sleep(1)
              
              # Create tracer (tracer overlap) merged for each cluster
              print(glue("[{start_time} -> {Sys.time()} | {i}/{n_perms} | {k}/{n_clusters}] tracer_overlap: Create merged seed map"))
              tracer_overlap_merged_seed_file_50um <- tempfile(pattern = "merged_seed-", tmpdir = "/dev/shm/", fileext = "_50um.mnc")
              tracer_overlap_merged_seed_file_200um <-gsub('_50um', '_200um', tracer_overlap_merged_seed_file_50um)
              
              n_tracer <- df_tracer %>%
                filter(cluster==k) %>%
                nrow()
              tracer_overlap_tracer_files <- df_tracer %>%
                sample_n(size = n_tracer, replace = F) %>%
                pull(tracer_file)
              tracer_overlap_seed_files <- df_tracer %>%
                sample_n(size = n_tracer, replace = F) %>%
                pull(seed_file)
              system(cmd <- glue("mincmath -2 -max {paste(tracer_overlap_seed_files, collapse = ' ')} {tracer_overlap_merged_seed_file_50um} -clobber"))
              system(cmd <- glue("mincresample -2 -like {likefile_200um} -keep_real_range {tracer_overlap_merged_seed_file_50um} {tracer_overlap_merged_seed_file_200um}"))
              
              Sys.sleep(1)
              
              #########################
              
              # Compute overlap statistics for seed-cluster overlap (tracer in cortex vs SC in cortex)
              print(glue("[{start_time} -> {Sys.time()} | {i}/{n_perms} | {k}/{n_clusters}] seed_overlap: Get vector"))
              seed_overlap_vol <- mincGetVolume(seed_overlap_merged_tracer_file_200um)
              seed_overlap_vec <- seed_overlap_vol[mask_vol_cortex > 0.5]
              seed_overlap_pos <- seed_overlap_vec >= tracer_thres
              
              # Cleanup
              print(glue("[{start_time} -> {Sys.time()} | {i}/{n_perms} | {k}/{n_clusters}] seed_overlap: Cleanup tempfiles"))
              file.remove(c(seed_overlap_merged_tracer_file_50um, seed_overlap_merged_tracer_file_200um))
              
              print(glue("[{start_time} -> {Sys.time()} | {i}/{n_perms} | {k}/{n_clusters}] seed_overlap: Compute stats"))
              seed_olperm_pos <- get_overlap_statistics(sc_pos, seed_overlap_pos) %>% 
                t %>%
                as_tibble() %>%
                mutate(comparison="seed_overlap",
                       sc_direction="positive",
                       n_clusters=n_clusters,
                       cluster=k,
                       run=i) %>%
                select(run, comparison, n_clusters, cluster, sc_direction, everything())
              seed_olperm_neg <- get_overlap_statistics(sc_neg, seed_overlap_pos) %>% 
                t %>%
                as_tibble() %>%
                mutate(comparison="seed_overlap",
                       sc_direction="negative",
                       n_clusters=n_clusters,
                       cluster=k,
                       run=i) %>%
                select(run, comparison, n_clusters, cluster, sc_direction, everything())
              
              # Compute overlap statistics for tracer-cluster overlap (seed in cortex vs SC in cortex)
              print(glue("[{start_time} -> {Sys.time()} | {i}/{n_perms} | {k}/{n_clusters}] tracer_overlap: Get vector"))
              tracer_overlap_vol <- mincGetVolume(tracer_overlap_merged_seed_file_200um)
              tracer_overlap_vec <- tracer_overlap_vol[mask_vol_cortex > 0.5]
              tracer_overlap_pos <- tracer_overlap_vec >= seed_thres
              
              # Cleanup
              print(glue("[{start_time} -> {Sys.time()} | {i}/{n_perms} | {k}/{n_clusters}] tracer_overlap: Cleanup tempfiles"))
              file.remove(c(tracer_overlap_merged_seed_file_50um, tracer_overlap_merged_seed_file_200um))
              
              print(glue("[{start_time} -> {Sys.time()} | {i}/{n_perms} | {k}/{n_clusters}] tracer_overlap: Compute stats"))
              tracer_olperm_pos <- get_overlap_statistics(sc_pos, tracer_overlap_pos) %>% 
                t %>%
                as_tibble() %>%
                mutate(comparison="tracer_overlap",
                       sc_direction="positive",
                       n_clusters=n_clusters,
                       cluster=k,
                       run=i) %>%
                select(run, comparison, n_clusters, cluster, sc_direction, everything())
              tracer_olperm_neg <- get_overlap_statistics(sc_neg, tracer_overlap_pos) %>% 
                t %>%
                as_tibble() %>%
                mutate(comparison="tracer_overlap",
                       sc_direction="negative",
                       n_clusters=n_clusters,
                       cluster=k,
                       run=i) %>%
                select(run, comparison, n_clusters, cluster, sc_direction, everything())
              
              # Bind information together
              print(glue("[{start_time} -> {Sys.time()} | {i}/{n_perms} | {k}/{n_clusters}] Bind together"))
              cluster_olperm <- rbind(seed_olperm_pos, seed_olperm_neg, tracer_olperm_pos, tracer_olperm_neg)
              cluster_olperm_list[[k]] <- cluster_olperm
            }
            
            print(glue("[{start_time} -> {Sys.time()} | {i}/{n_perms}] Binding everything together"))
            olperm <- rbindlist(cluster_olperm_list)
            # perm_list[[i]] <- olperm
            
            #setTxtProgressBar(prog, i)
            olperm
          }

print(glue("[{Sys.time()}] Stopping cluster"))

#close(prog)
registerDoSEQ()
stopCluster(cl)

print(glue("[{Sys.time()}] Binding"))
olperm_all <- rbindlist(perm_list)
rm(perm_list)

print(glue("[{Sys.time()}] Writing out"))
write_csv(olperm_all, glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/projection_overlaps/merged_overlap_permutations.csv"))