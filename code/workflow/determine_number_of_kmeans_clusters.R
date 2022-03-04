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
# mode <- "residual"

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
mode <- args[2]

# Input files
mask_file_thalamus <- glue("{outdir}/thalamus_mask_highres.mnc")
mask_file_cortex <- glue("{outdir}/cortex_mask_lowres.mnc")
sc_file <- glue("{outdir}/correlation_{mode}.npy")

# Output csv
outcsv <- glue("{outdir}/decomposition_{mode}/kmeans/scree.csv")
outpng_wss <- glue("{outdir}/decomposition_{mode}/kmeans/scree_wss.png")
outpng_wssoverbss <- glue("{outdir}/decomposition_{mode}/kmeans/scree_wssoverbss.png")

###########################
# Read data
###########################

print(glue("[{Sys.time()}] Reading mask data"))

mask_vol_thalamus <- mincGetVolume(mask_file_thalamus)
mask_vol_cortex <- mincGetVolume(mask_file_cortex)

p_thal <- length(mask_vol_thalamus[mask_vol_thalamus > 0.5])
p_cortex <- length(mask_vol_cortex[mask_vol_cortex > 0.5])

# Takes a 3-5 minutes to read in all the data
print(glue("[{Sys.time()}] Reading SC matrix"))
scmat <- npyLoad(sc_file)


###########################
# Compute cluster distances
###########################

# Compute total sum of squares
print(glue("[{Sys.time()}] Computing total sum of squares"))
tss <- sum(scale(scmat, scale=F)^2)

# Loop over clusters
print(glue("[{Sys.time()}] Computing within sum of squares"))
n_clust_list <- seq(from=2, to=42, by=2)

df <- tibble(n_clusters=integer(length(n_clust_list)),
             wss=numeric(length(n_clust_list)),
             tss=tss,
             bss=numeric(length(n_clust_list)))

for (k in 1:length(n_clust_list)) {
  n_clusters <- n_clust_list[[k]]
  print(glue("Working on number of clusters: {n_clusters}"))
  
  labels_file <- glue("{outdir}/decomposition_{mode}/kmeans/{n_clusters}_clusters/cluster_labels.npy")
  labels <- npyLoad(labels_file)
  
  # Compute within some of squares
  prog <- progress_estimated(n_clusters)
  wss <- 1:n_clusters %>% 
    map_dbl(function(i) {
      out <- sum(scale(scmat[labels==i,], scale=F)^2)
      print(prog$tick())
      out
    }) %>% 
    sum()
  
  # Compute between sum of squares
  bss <- tss - wss
  
  # Output
  df$n_clusters[[k]] <- n_clusters
  df$wss[[k]] <- wss
  df$bss[[k]] <- tss - wss
  
  write_csv(df, outcsv)
}

# Writing out results
print(glue("[{Sys.time()}] Writing dataframe"))
write_csv(df, outcsv)

print(glue("[{Sys.time()}] Outputting plots"))

png(outpng_wss, width = 1024, height = 768, res = 300)
print(ggplot(df, aes(x=n_clusters, y=(wss))) +
        geom_point() + 
        geom_path() +
        xlab("Number of clusters") +
        ylab("Total within sum of squares") +
        theme_bw()
)
dev.off()

png(outpng_wssoverbss, width = 1024, height = 768, res = 300)
print(ggplot(df, aes(x=n_clusters, y=(wss/bss))) +
        geom_point() + 
        geom_path() +
        xlab("Number of clusters") +
        ylab("Total within sum of squares /\n between sum of squares") +
        theme_bw()
)
dev.off()

# Done
print(glue("[{Sys.time()}] Done!"))



