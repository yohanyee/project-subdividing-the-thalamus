#!/bin/bash
# Run this to redo the analysis
# Run each line manually and examine the results, as opposed to executing this file

#############################################
# Setup
#############################################

# Working directory
export BASEDIR="/path/to/this/repository"
cd $BASEDIR

# Mouse Imaging Centre specifics
module load mice-env/1.0.8 # Require module for numpy / itk_convert used by transform_space.py
module load python-modules/20200317 # For scikit sklearn

#############################################
# PREPROCESSING
#############################################

# Get masks
#   Download from ABI server and process at both 50um and 200um
#   The cortex laplacian is downloaded (only available) at 10um, transforming this to MINC requires ~30GB, takes 5 minutes or so
/usr/bin/time -v /bin/bash code/workflow/get_masks.sh

# Get input files
#   Read scanbase csv (for 154 mice) and select appropriate columns, point to input files,
#   point to Jacobian determinants that are already resampled, and get overall volumes
/usr/bin/time -v Rscript code/workflow/get_input_files.R

#############################################
# MAIN COMPUTATION
# 
# For a given thalamus (50um) mask and cortex (50um) mask
# Set up an output directory, symlink the masks to that specific directory
# Idea is to keep the output modular, can replace a different mask if needed, 
#   e.g. the opposite hemisphere, or cerebral cortex instead of isocortex
#
# Therefore, all programs below should take single project directory containing masks and work with that
#############################################

# First, setup this project output directory
MASKDIR="data/outputs/masks"
OUTDIR="data/outputs/isocortex_thalamus_left"
mkdir -p $OUTDIR

# Choose masks and transcriptomic data to use
ln -s `realpath $MASKDIR/TH_coronal_coverage_50um_left.mnc` $OUTDIR/thalamus_mask_highres.mnc
ln -s `realpath $MASKDIR/TH_coronal_coverage_200um_left.mnc` $OUTDIR/thalamus_mask_lowres.mnc
ln -s `realpath $MASKDIR/Isocortex_coronal_coverage_200um_left.mnc` $OUTDIR/cortex_mask_lowres.mnc
ln -s `realpath $MASKDIR/Isocortex_coronal_coverage_50um_left.mnc` $OUTDIR/cortex_mask_highres.mnc
ln -s `realpath $MASKDIR/laplacian_50um_under_Isocortex_coronal_coverage_50um_left.mnc` $OUTDIR/cortex_laplacian_highres.mnc
ln -s `realpath $MASKDIR/laplacian_200um_under_Isocortex_coronal_coverage_200um_left.mnc` $OUTDIR/cortex_laplacian_lowres.mnc
ln -s `realpath data/resources/average_template_50.mnc` $OUTDIR/average_template_50.mnc
ln -s `realpath data/resources/average_template_200.mnc` $OUTDIR/average_template_200.mnc
ln -s `realpath data/resources/coronal_200um_gene_mtx_normalized_filtered.npz` $OUTDIR/filtered_normalized_expression.npz
ln -s `realpath data/resources/coronal_200um_gene_table_normalized_filtered.csv` $OUTDIR/filtered_normalized_expression_gene_list.csv
ln -s `realpath data/resources/coronal_200um_coverage_bin0.8.mnc` $OUTDIR/filtered_normalized_expression_coverage.mnc

# Get WT154 matrices
#   Return n x p_thal matrix for thalamus at 50um, 
#   n x p_cortex matrix for cortex at 200um 
#   consisting of absolute Jacobian determinant values.
#   Also return normalized versions of above data matrices, regressing out group and structure volume
#   10 minutes to run 
/usr/bin/time -v Rscript code/workflow/get_data_matrices_wt154.R $OUTDIR

# Compute and store p_thal x p_cortex (Fisher-transformed) correlation matrix from normalized data
#   Should take around 5 minutes to run
/usr/bin/time -v Rscript code/workflow/compute_sc.R $OUTDIR

#############################################
# OPTIONAL COMPUTATION
#
# If you want to regress out distance and transcriptomic, also run the following
#############################################

# Get distance matrix (p_thal x p_matric)
#   Should take 5-10 mins to run locally
/usr/bin/time -v Rscript code/workflow/compute_distance.R $OUTDIR

# Get transcriptomic similarity matrix (p_thal x p_cortex)
#   Run on compute cluster using 20 cores using command below, should take about 3-4 hours to run
#   Or run locally: /usr/bin/time -v Rscript code/workflow/compute_transcriptomic_similarity.R $OUTDIR 4
#   Resulting matrix should be 5 GB for single hemisphere at 50um thalamic voxels
echo "/usr/bin/time -v Rscript code/workflow/compute_transcriptomic_similarity.R data/outputs/isocortex_thalamus_left 20" | qbatch --mem 95G -w 5:58:59 --ppj 20 -N thalamus_ts --logdir code/workflow/logs -

# Compute and store residual correlation matrix along with the following parameters for each thalamic voxel:
#   betas, t-stat, p-value for: TS, distance; and also R-squared, AIC, and BIC
#   Should take around 20 mins
/usr/bin/time -v Rscript code/workflow/compute_sc_residuals.R $OUTDIR

#############################################
# MAIN CLUSTERING / THALAMUS PARCELLATION BASED ON STRUCTURAL COVARIANCE
#
# k-means cluster thalamic voxels and for k=2..40 and store cluster labels and mean cortical values 
# Note: AMBA has 79(!) unique labels within the thalamus (mnc file)
# Levels (starting at brain=1, gray=2, ...) are: 
#   level 5: 1 structure (thalamus)
#         6: 3 structures (sensory vs polymodal)
#         7: 15 structures
#         8: 45 structures
# Select optimal k based on both scree
#############################################

# k-means cluster: 20 sets of clusters starting at 2 and increment by 2
for k in {2..42..2}
do
    # Initialization for 5 jobs can take an hour, so expect script to take up to 5 hours total, for 5 threads (n_init=10, max_iter=500) (assuming each iter takes around 10 seconds)
    # Cluster thalamic data in SC matrix
    echo "/usr/bin/time -v python3 code/workflow/decompose_kmeans.py $OUTDIR ${OUTDIR}/decomposition_z/kmeans/${k}_clusters correlation_z.npy ${k}" | qbatch --mem 90G -w 5:58:59 --ppj 5 -N km_z_${k} --logdir code/workflow/logs -
    # Optional: cluster thalamic data in SC residual (distance and transcriptomic similarity regressed out) matrix
    echo "/usr/bin/time -v python3 code/workflow/decompose_kmeans.py $OUTDIR ${OUTDIR}/decomposition_residual/kmeans/${k}_clusters correlation_residual.npy ${k}" | qbatch --mem 90G -w 5:58:59 --ppj 5 -N km_r_${k} --logdir code/workflow/logs -
done

# Determine optimal number of components based on scree plot
#   Estimated to use 22GB, 10-20 mins for each script
# Optimal number of components in SC matrix
/usr/bin/time -v Rscript code/workflow/determine_number_of_kmeans_clusters.R ${OUTDIR} z
# Optional: Optimal number of components in SC residual (distance and transcriptomic similarity regressed out) matrix
/usr/bin/time -v Rscript code/workflow/determine_number_of_kmeans_clusters.R ${OUTDIR} residual

#############################################
# ANALYSIS: OVERLAP WITH ABI ATLAS
#############################################

# Compute overlap data locally
#   Takes about 10 minutes to run, ~ 2GB
for k in {2..42..2}
do
    # 2 clusters: 12 seconds, 1.25G
    # 16 clusters: 27 seconds, 1.25G
    # 20 clusters: 34 seconds, 1.25G
    # 42 clusters: 1:14, 1.26G
    echo "Working on n_clusters = ${k}"
    /usr/bin/time -v Rscript code/workflow/compute_ABI_atlas_overlap_with_clusters.R $OUTDIR z ${k} 2>&1 | tee code/workflow/logs/ol_ABI_${k}_z.log 
    /usr/bin/time -v Rscript code/workflow/compute_ABI_atlas_overlap_with_clusters.R $OUTDIR residual ${k} 2>&1 | tee code/workflow/logs/ol_ABI_${k}.log #Optional
done

# Compute overlap permutation data on compute cluster
#   NOTE: the non-permuted overlap data must be first computed, since this script below relies on reading in minc files (*_fixed.mnc) that 
#   is written out by the above script
for k in {2..16..2}
do
    # Submit this on compute cluster
    # 2 clusters: approximately 6000 permutations / hour using 20 threads 
    # 16 clusters: approximately 2000 permutations / hour using 20 threads
    # 42 clusters: approximately 800 permutations / hour using 20 threads
    # Writes out a csv file with permutation distribution that is 130MB*(${k}) in size
    echo "/usr/bin/time -v Rscript code/workflow/compute_ABI_atlas_overlap_with_clusters_permutation.R $OUTDIR z ${k} 10000 20" | qbatch --mem 16G -w 23:58:59 --ppj 20 -N ol_ABI_perm_${k}_z --logdir code/workflow/logs -
    echo "/usr/bin/time -v Rscript code/workflow/compute_ABI_atlas_overlap_with_clusters_permutation.R $OUTDIR residual ${k} 10000 20" | qbatch --mem 16G -w 23:58:59 --ppj 20 -N ol_ABI_perm_${k} --logdir code/workflow/logs - #Optional
done

# Compute p-values
#   Generally quick, script runs in 5-10 minutes
#   Specifically, it takes about 0.5-1.5 minutes * ${k}
#   16 clusters: 21 minutes
#   Requires enough memory to read in permutation file 
#   Script shouldn't use more than 8 GB even when number of clusters is maximal (2GB + 130*${k})
#   Can only be run when the above permutation data has been generated
for k in {2..16..2}
do
    /usr/bin/time -v Rscript code/workflow/compute_ABI_atlas_overlap_with_clusters_pvalues.R $OUTDIR z ${k} 2>&1 | tee code/workflow/logs/ol_ABI_p_${k}_z.log
    /usr/bin/time -v Rscript code/workflow/compute_ABI_atlas_overlap_with_clusters_pvalues.R $OUTDIR residual ${k} 2>&1 | tee code/workflow/logs/ol_ABI_p_${k}.log #Optional
done

#############################################
# ANALYSIS: OVERLAP WITH PROJECTION (MOUSE CONNECTIVITY) DATA
# Compute projection data overlap
# 
# Overlaps measured in two different ways, corresponding to each end of the tracer data (seed and target):
# 
# For seed overlaps (percent of seed covered by cluster)
# - seed (injection fraction) binarized at 0.9 ("confident that seed is where tracer was injected")
# - seed overlaps csv filtered at 0.01
# - seed usable data thresholded at 0.9 ("seed is specific to cluster")
#
# For tracer overlaps (percent of cluster covered by tracer, since tracer may project elsewhere too)
# - tracer (projection density) binarized at 0.1 ("there is a measureable tracer signal")
# - tracer overlaps csv filtered at 0.01 
# - tracer usable data thresholded at 0.1 ("enough of the cluster is covered by the tracer to count it as an overlap")
#
# If running on compute cluster, you may run into graphics issues, in which case you should keep the plotting part commented

# Compute overlaps
for k in {2..16..2}
do
    for l in `seq 1 $k`
    do
        echo "Working on number of clusters: ${k}, cluster number: ${l}"
        # Compute overlaps with clusters built from SC matrix
        echo "/usr/bin/time -v Rscript code/workflow/compute_projection_seed_overlap_with_cluster.R $OUTDIR z ${k} ${l}" | qbatch --mem 8G -w 1:58:59 --ppj 1 -N ol_seed_${k}_${l}_z --logdir code/workflow/logs -
        echo "/usr/bin/time -v Rscript code/workflow/compute_projection_tracer_overlap_with_cluster.R $OUTDIR z ${k} ${l}" | qbatch --mem 8G -w 1:58:59 --ppj 1 -N ol_tracer_${k}_${l}_z --logdir code/workflow/logs -
        # Optional: Compute overlaps with clusters built from SC residual (distance and transcriptomic similarity regressed out) matrix
        echo "/usr/bin/time -v Rscript code/workflow/compute_projection_seed_overlap_with_cluster.R $OUTDIR residual ${k} ${l}" | qbatch --mem 8G -w 1:58:59 --ppj 1 -N ol_seed_${k}_${l} --logdir code/workflow/logs -
        echo "/usr/bin/time -v Rscript code/workflow/compute_projection_tracer_overlap_with_cluster.R $OUTDIR residual ${k} ${l}" | qbatch --mem 8G -w 1:58:59 --ppj 1 -N ol_tracer_${k}_${l} --logdir code/workflow/logs -
    done
done

# Compute overlap statistics for merged tracers, merging experiments specific to each cluster
#   Relies on seed_overlap_cluster_{k}.csv and tracer_overlap_cluster_{k}.csv files in projection_overlaps/${k} folder being generated from scripts above
#   Run on compute cluster in parallel (for generating permutations):
echo "/usr/bin/time -v Rscript code/workflow/compute_projection_merged_overlap_statistics_and_permutations.R $OUTDIR z 6" | qbatch --mem 24G -w 12:58:59 --ppj 20 -N ol_merged_perm_6_z --logdir code/workflow/logs -
echo "/usr/bin/time -v Rscript code/workflow/compute_projection_merged_overlap_statistics_and_permutations.R $OUTDIR residual 6" | qbatch --mem 24G -w 12:58:59 --ppj 20 -N ol_merged_perm_6 --logdir code/workflow/logs - # Optional


# Compute individual overlap statistics
#   Also relies on seed_overlap_cluster_{k}.csv and tracer_overlap_cluster_{k}.csv files in projection_overlaps/${k} folder being generated from scripts above
/usr/bin/time -v Rscript code/workflow/compute_projection_individual_overlap_statistics.R $OUTDIR z 6 2>&1 | tee code/workflow/logs/ol_individual_projectionstats_6_z.log
/usr/bin/time -v Rscript code/workflow/compute_projection_individual_overlap_statistics.R $OUTDIR residual 6 2>&1 | tee code/workflow/logs/ol_individual_projectionstats_6.log # Optional

