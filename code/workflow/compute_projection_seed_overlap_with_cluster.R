# Takes about 1 hour mins to run

library(tidyverse)
library(RMINC)
library(glue)

query <- read_csv("data/resources/projection/query.csv")
experiment_dir <- "data/resources/projection/GridData_50um"

# Arguments
print(glue("[{Sys.time()}] Input arguments"))

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
mode <- args[2]
n_clusters <- as.integer(args[3])
clust <- as.integer(args[4])

#indir <- "data/outputs/isocortex_thalamus_left"
#mode <- "residual"
#n_clusters <- 6
#clust <- 1

# Create output directories
print(glue("[{Sys.time()}] Creating output directories"))
outdir <- glue("{indir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/projection_overlaps/{clust}")
dir.create(outdir, showWarnings = T, recursive = T)

# Get cluster definition
print(glue("[{Sys.time()}] Getting cluster definitions"))
roifile <- glue("{indir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/cluster_mask_{clust}.mnc")
roivol <- mincGetVolume(roifile)
roisize <- sum(roivol)

# Get queries of interest
print(glue("[{Sys.time()}] Getting ABI experiments of interest"))
qids <- query %>%
  filter(`product-id`==5)

# Output overlaps file
print(glue("[{Sys.time()}] Computing overlaps"))
overlaps_file <- glue("{outdir}/seed_overlap_cluster_{clust}.csv")

# Compute overlaps if the file doesn't exist
if (!file.exists(overlaps_file)) {
  out_list <- list()
  k <- 1
  prog <- txtProgressBar(max=nrow(qids), style=3)
  for (i in 1:nrow(qids)) {
    id <- qids$id[i]
    sfile<- glue("{experiment_dir}/{id}/masks/injection_fraction_bin0.9.mnc")
    sfile_flipped <- glue("{experiment_dir}/{id}/masks/injection_fraction_bin0.9_flipped.mnc")
    
    if (file.exists(sfile)) {
      svol <- mincGetVolume(sfile)
      ssize <- sum(svol)
      
      if (ssize >= 1) {
        voxels_seed_covered_by_roi <- sum(roivol[svol > 0.5])
        voxels_roi_covered_by_seed <- sum(svol[roivol > 0.5])
        percent_of_seed_covered_by_roi <- voxels_seed_covered_by_roi/ssize
        percent_of_roi_covered_by_seed <- voxels_roi_covered_by_seed/roisize
        out_list[[k]] <- tibble(experiment=id,
                                hemisphere="right",
                                seed_size=ssize,
                                roi_size=roisize,
                                voxels_seed_covered_by_roi=voxels_seed_covered_by_roi,
                                voxels_roi_covered_by_seed=voxels_roi_covered_by_seed,
                                percent_of_seed_covered_by_roi=percent_of_seed_covered_by_roi,
                                percent_of_roi_covered_by_seed=percent_of_roi_covered_by_seed
        )
        k <- k+1
      }
    }
    
    if (file.exists(sfile_flipped)) {
      svol <- mincGetVolume(sfile_flipped)
      ssize <- sum(svol)
      
      if (ssize >= 1) {
        voxels_seed_covered_by_roi <- sum(roivol[svol > 0.5])
        voxels_roi_covered_by_seed <- sum(svol[roivol > 0.5])
        percent_of_seed_covered_by_roi <- voxels_seed_covered_by_roi/ssize
        percent_of_roi_covered_by_seed <- voxels_roi_covered_by_seed/roisize
        out_list[[k]] <- tibble(experiment=id,
                                hemisphere="left",
                                seed_size=ssize,
                                roi_size=roisize,
                                voxels_seed_covered_by_roi=voxels_seed_covered_by_roi,
                                voxels_roi_covered_by_seed=voxels_roi_covered_by_seed,
                                percent_of_seed_covered_by_roi=percent_of_seed_covered_by_roi,
                                percent_of_roi_covered_by_seed=percent_of_roi_covered_by_seed
        )
        k <- k + 1
      }
    }
    
    setTxtProgressBar(prog, i)
  }
  overlaps_df <- bind_rows(out_list)
  write_csv(overlaps_df, overlaps_file)
} else {
  overlaps_df <- read_csv(overlaps_file)
}
close(prog)

# # Visualize distribution of overlaps
# print(glue("[{Sys.time()}] Visualizing overlap distribution (above 0.01 overlap)"))
# ol_hist <- overlaps_df %>%
#   filter(percent_of_seed_covered_by_roi >= 0.01) %>%
#   ggplot(aes(x=percent_of_seed_covered_by_roi)) + 
#   geom_histogram(binwidth = 0.05, color='grey20', fill='grey50') +
#   xlab('Percent of injection seed covered by cluster ROI') +
#   ylab('Number of projection images\n(filtered above 0.01)') +
#   theme_bw()
# 
# # Save visualization
# print(glue("[{Sys.time()}] Saving histogram"))
# png(glue("{outdir}/coverage_of_seed.png"), width = 2000, height=1000, res=300)
# print(ol_hist)
# dev.off()

# Filter away low overlaps
print(glue("[{Sys.time()}] Processing dataframe"))
odf <- overlaps_df %>%
  filter(percent_of_seed_covered_by_roi >= 0.01) %>%
  mutate(tracer_file=case_when(hemisphere=="right" ~ glue("{experiment_dir}/{experiment}/projection_density.mnc"),
                               hemisphere=="left" ~ glue("{experiment_dir}/{experiment}/projection_density_flipped.mnc")),
         seed_file=case_when(hemisphere=="right" ~ glue("{experiment_dir}/{experiment}/injection_fraction.mnc"),
                             hemisphere=="left" ~ glue("{experiment_dir}/{experiment}/injection_fraction_flipped.mnc"))
  )

# Merge with query data
odf <- odf %>%
  left_join(query, by=c("experiment"="id"))

# Write out csv
print(glue("[{Sys.time()}] Writing overlap (above 0.01) data"))
odf %>% 
  write_csv(overlaps_file)

##################################################
# Work with usable tracer experiments
# Usable -> specific to cluster, i.e. atleast 90% overlap

print(glue("[{Sys.time()}] Working with usable data (overlap above 0.9)"))

odf_usable <- odf %>%
  filter(percent_of_seed_covered_by_roi >= 0.9)

# Symlink individual projection seed and targets
print(glue("[{Sys.time()}] Creating overlap directory"))
linkdir <- glue("{outdir}/seed_overlap_experiments")
dir.create(linkdir)

if (nrow(odf_usable) >= 1) {
  print(glue("[{Sys.time()}] Symlinking overlapping projection data"))
  prog <- txtProgressBar(max=nrow(odf_usable), style=3)
  for (i in 1:nrow(odf_usable)) {
    file.symlink(normalizePath(odf_usable$seed_file[[i]]), glue("{linkdir}/{odf_usable$experiment[[i]]}_{odf_usable$hemisphere[[i]]}_seed.mnc"))
    file.symlink(normalizePath(odf_usable$tracer_file[[i]]), glue("{linkdir}/{odf_usable$experiment[[i]]}_{odf_usable$hemisphere[[i]]}_tracer.mnc"))
    setTxtProgressBar(prog, i)
  }
  close(prog)
  
  # Save merged projection seed and target
  print(glue("[{Sys.time()}] Max merging projection data"))
  system(cmd <- glue("mincmath -2 -max {linkdir}/*_seed.mnc {outdir}/seed_overlap_merged_max_seed.mnc"))
  system(cmd <- glue("mincmath -2 -max {linkdir}/*_tracer.mnc {outdir}/seed_overlap_merged_max_tracer.mnc"))
} else {
  print("No experiments in which seed overlaps with cluster ROI more than 90%!")
}

###################################################
# Done!