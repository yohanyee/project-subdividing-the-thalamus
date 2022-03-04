#!/usr/bin/Rscript
#
# Submit with qbatch: echo "Rscript code/workflow/compute_transcriptomic_similarity [outdir] [n threads]" | qbatch --mem [MEM]G -w 5:58:59 --ppj [N_THREADS] -N thalamus_ts -
# Where N_THREADS is the number of requested threads
# and MEM = 10G (for the storage and saving of TS matrix) + 5 GB (for the base R) + N_THREADS * 4 GB
#
# So for example, a high performance computing submission would be:
#
# cd $BASEDIR
# module load mice-env/1.0.8
# echo "/usr/bin/time -v Rscript code/workflow/compute_transcriptomic_similarity data/outputs/isocortex_thalamus_left 20" | qbatch --mem 95G -w 5:58:59 --ppj 20 -N thalamus_ts --logdir code/workflow/logs -
#

# Load libraries
library(tidyverse)
library(glue)
library(RMINC)
library(RcppCNPy)
library(foreach)
library(doParallel)
library(data.table)

###########################
# Arguments
###########################

# Testing arguments
# outdir <- "data/outputs/isocortex_thalamus_left"
# n_threads <- 4

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
n_threads <- as.integer(args[2])

# Input files
mask_file_thalamus <- glue("{outdir}/thalamus_mask_highres.mnc")
mask_file_cortex <- glue("{outdir}/cortex_mask_lowres.mnc")
mask_file_coverage <- glue("{outdir}/filtered_normalized_expression_coverage.mnc")

gexp_file <- glue("{outdir}/filtered_normalized_expression.npz")
gene_file <- glue("{outdir}/filtered_normalized_expression_gene_list.csv")

###########################
# Read data
###########################

print(glue("[{Sys.time()}] Reading data"))

mask_vol_thalamus <- mincGetVolume(mask_file_thalamus)
mask_vol_cortex <- mincGetVolume(mask_file_cortex)
mask_vol_coverage <- mincGetVolume(mask_file_coverage)

p_thal <- length(mask_vol_thalamus[mask_vol_thalamus > 0.5])
p_cortex <- length(mask_vol_cortex[mask_vol_cortex > 0.5])

gexp <- npyLoad(gexp_file)
genes <- read_csv(gene_file) %>% select(symbol=msg.genes.original_symbol, name=msg.genes.original_name, entrez=msg.genes.entrez_id)

dim(genes)

###########################

thalamus_mask_indices <- which(mask_vol_thalamus > 0.5)
coverage_mask_indices <- which(mask_vol_coverage > 0.5)
gexp_cortical_indices <- which(mask_vol_cortex[mask_vol_coverage > 0.5] > 0.5)

gexp_cortex <- gexp[,gexp_cortical_indices]

dim(gexp)
dim(gexp_cortex)

############################

print(glue("[{Sys.time()}] Setting up multithreading"))

cl <- makeCluster(n_threads, outfile="", rscript_args="--vanilla")
registerDoParallel(cl) 

export_objects <- c("mincVectorToVoxelCoordinates", "mincConvertVoxelToWorld", "mincConvertWorldToVoxel", "minc.dimensions.sizes")

# Loop over thalamic voxels
n_iterations <- length(thalamus_mask_indices)

# For testing, uncomment below line
# n_iterations <- 100

print(glue("[{Sys.time()}] Running loop!"))

prog <- txtProgressBar(max=n_iterations, style = 3)
out <- foreach(i=1:n_iterations, 
               .export = c(export_objects), 
               .verbose = TRUE, 
               .inorder = TRUE) %dopar% {
                 
                 # Get voxel coords at highres
                 # mincVectorToVoxelCoordinates returns z, y, x when seen in Display
                 voxel_coords_highres <- mincVectorToVoxelCoordinates(mask_file_thalamus, thalamus_mask_indices[i])
                 
                 # Get world coordinate
                 # Since mincConvertVoxelToWorld takes in data as z, y, x, and the input voxel_coords here is given in this order, the output is in the correct order
                 world_coords <- mincConvertVoxelToWorld(mask_file_thalamus, voxel_coords_highres[1], voxel_coords_highres[2], voxel_coords_highres[3])
                 
                 # Get voxel coordinate at 200um
                 # Input coordinates as x, y, z; returns voxel coords as z, y, x
                 voxel_coords <- mincConvertWorldToVoxel(mask_file_coverage, world_coords[1], world_coords[2], world_coords[3])
                 
                 # Also as z, y, x
                 sizes <- minc.dimensions.sizes(mask_file_coverage)
                 vector_coords <- (sizes[3]*sizes[2])*voxel_coords[1] + (sizes[3]*voxel_coords[2]) + voxel_coords[3] + 1
                 
                 # Confirm the coordinates match
                 # all(mincVectorToVoxelCoordinates(mask_file_coverage, vectorCoord = vector_coords)==voxel_coords)
                 expression_matrix_index <- which(coverage_mask_indices==vector_coords)
                 
                 seed_exp_vec <- gexp[,expression_matrix_index]
                 coexpression <- apply(gexp_cortex, MARGIN = 2, FUN=cor, seed_exp_vec, use="pairwise.complete.obs")
                 
                 # Test that coexpression mapping is correct
                 # coexpression <- apply(gexp, MARGIN = 2, FUN=cor, seed_exp_vec, use="pairwise.complete.obs")
                 # ovol <- mask_vol_coverage
                 # ovol[] <- 0
                 # ovol[mask_vol_coverage > 0.5] <- coexpression
                 # mincWriteVolume(ovol, glue("tmp/coexpression_test_at_world_{world_coords[1]}_{world_coords[2]}_{world_coords[3]}.mnc"))
                 
                 setTxtProgressBar(prog, i)
                 
                 coexpression
               }

print(glue("[{Sys.time()}] Stopping cluster"))

close(prog)
registerDoSEQ()
stopCluster(cl)

# Bind data together
print(glue("[{Sys.time()}] Binding data output"))
outdf <- matrix(NA, nrow=p_thal, ncol=p_cortex)
prog <- txtProgressBar(max=p_thal, style=3)
for (l in 1:p_thal) {
  outdf[l, ] <- out[[l]]
  setTxtProgressBar(prog, l)
}
close(prog)

# Garbage collection
print(glue("[{Sys.time()}] Garbage collection"))
rm(out)
gc()

# Write out coexpression file
print(glue("[{Sys.time()}] Writing out data"))

npySave(filename = glue("{outdir}/coexpression.npy"), object = outdf)

print(glue("[{Sys.time()}] Done!"))