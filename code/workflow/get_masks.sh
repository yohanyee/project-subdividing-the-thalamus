#!/bin/bash

# Make output directory for masks
MASKDIR=$BASEDIR/data/outputs/masks
mkdir -p $MASKDIR

#######################################################################################################
# Structure mask processing

# Download ABI structure masks at 50 um
# Definitions: http://api.brain-map.org/api/v2/structure_graph_download/1.json
# Download server: http://download.alleninstitute.org/informatics-archive/current-release/
wget -O $MASKDIR/structure_315.nrrd http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/structure_masks/structure_masks_50/structure_315.nrrd # Isocortex, Isocortex / 70FF71 
wget -O $MASKDIR/structure_549.nrrd http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/structure_masks/structure_masks_50/structure_549.nrrd # Thalamus, TH / FF7080
wget -O $MASKDIR/structure_688.nrrd http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/structure_masks/structure_masks_50/structure_688.nrrd # Cerebral cortex, CTX / B0FFB8 (includes hippocampus)

# Download ABI structure masks at 25 um, for high res plots
wget -O $MASKDIR/structure_315_25um.nrrd http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/structure_masks/structure_masks_25/structure_315.nrrd # Isocortex, Isocortex / 70FF71 
wget -O $MASKDIR/structure_549_25um.nrrd http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/structure_masks/structure_masks_25/structure_549.nrrd # Thalamus, TH / FF7080
wget -O $MASKDIR/structure_688_25um.nrrd http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/structure_masks/structure_masks_25/structure_688.nrrd # Cerebral cortex, CTX / B0FFB8 (includes hippocampus)

# Convert to 50um data to MINC
python3 code/scripts/transform_space.py $MASKDIR/structure_315.nrrd $MASKDIR/Isocortex_50um.mnc
python3 code/scripts/transform_space.py $MASKDIR/structure_549.nrrd $MASKDIR/TH_50um.mnc
python3 code/scripts/transform_space.py $MASKDIR/structure_688.nrrd $MASKDIR/CTX_50um.mnc

# Downsample 50um data
mincresample -2 -like data/resources/average_template_200.mnc -nearest_neighbour -keep_real_range $MASKDIR/Isocortex_50um.mnc $MASKDIR/Isocortex_200um.mnc
mincresample -2 -like data/resources/average_template_200.mnc -nearest_neighbour -keep_real_range $MASKDIR/TH_50um.mnc $MASKDIR/TH_200um.mnc
mincresample -2 -like data/resources/average_template_200.mnc -nearest_neighbour -keep_real_range $MASKDIR/CTX_50um.mnc $MASKDIR/CTX_200um.mnc

# Convert to 25um data to MINC
python3 code/scripts/transform_space.py $MASKDIR/structure_315_25um.nrrd $MASKDIR/Isocortex_25um.mnc
python3 code/scripts/transform_space.py $MASKDIR/structure_549_25um.nrrd $MASKDIR/TH_25um.mnc
python3 code/scripts/transform_space.py $MASKDIR/structure_688_25um.nrrd $MASKDIR/CTX_25um.mnc

# Intersect each mask with both coronal and sagittal coverage maps to produce 6 masks for potential use
# (only 2 of these will be used, one for thalamus and one for cortex, related to either coronal or sagittal coverage
minccalc -2 -expression "A[0]>0.5 && A[1]>0.5" data/resources/sagittal_200um_coverage_bin0.8_resampled_50um.mnc $MASKDIR/Isocortex_50um.mnc $MASKDIR/Isocortex_sagittal_coverage_50um.mnc
minccalc -2 -expression "A[0]>0.5 && A[1]>0.5" data/resources/sagittal_200um_coverage_bin0.8_resampled_50um.mnc $MASKDIR/TH_50um.mnc $MASKDIR/TH_sagittal_coverage_50um.mnc
minccalc -2 -expression "A[0]>0.5 && A[1]>0.5" data/resources/sagittal_200um_coverage_bin0.8_resampled_50um.mnc $MASKDIR/CTX_50um.mnc $MASKDIR/CTX_sagittal_coverage_50um.mnc

minccalc -2 -expression "A[0]>0.5 && A[1]>0.5" data/resources/coronal_200um_coverage_bin0.8_resampled_50um.mnc $MASKDIR/Isocortex_50um.mnc $MASKDIR/Isocortex_coronal_coverage_50um.mnc
minccalc -2 -expression "A[0]>0.5 && A[1]>0.5" data/resources/coronal_200um_coverage_bin0.8_resampled_50um.mnc $MASKDIR/TH_50um.mnc $MASKDIR/TH_coronal_coverage_50um.mnc
minccalc -2 -expression "A[0]>0.5 && A[1]>0.5" data/resources/coronal_200um_coverage_bin0.8_resampled_50um.mnc $MASKDIR/CTX_50um.mnc $MASKDIR/CTX_coronal_coverage_50um.mnc

# Also do at 200um
minccalc -2 -expression "A[0]>0.5 && A[1]>0.5" data/resources/sagittal_200um_coverage_bin0.8.mnc $MASKDIR/Isocortex_200um.mnc $MASKDIR/Isocortex_sagittal_coverage_200um.mnc
minccalc -2 -expression "A[0]>0.5 && A[1]>0.5" data/resources/sagittal_200um_coverage_bin0.8.mnc $MASKDIR/TH_200um.mnc $MASKDIR/TH_sagittal_coverage_200um.mnc
minccalc -2 -expression "A[0]>0.5 && A[1]>0.5" data/resources/sagittal_200um_coverage_bin0.8.mnc $MASKDIR/CTX_200um.mnc $MASKDIR/CTX_sagittal_coverage_200um.mnc

minccalc -2 -expression "A[0]>0.5 && A[1]>0.5" data/resources/coronal_200um_coverage_bin0.8.mnc $MASKDIR/Isocortex_200um.mnc $MASKDIR/Isocortex_coronal_coverage_200um.mnc
minccalc -2 -expression "A[0]>0.5 && A[1]>0.5" data/resources/coronal_200um_coverage_bin0.8.mnc $MASKDIR/TH_200um.mnc $MASKDIR/TH_coronal_coverage_200um.mnc
minccalc -2 -expression "A[0]>0.5 && A[1]>0.5" data/resources/coronal_200um_coverage_bin0.8.mnc $MASKDIR/CTX_200um.mnc $MASKDIR/CTX_coronal_coverage_200um.mnc

# Create a cortical-thalamic label set
minccalc -2 -expression "A[0] + A[1] + A[2] + 3*A[3]" data/resources/average_template_50_mask.mnc $MASKDIR/CTX_50um.mnc $MASKDIR/Isocortex_50um.mnc $MASKDIR/TH_50um.mnc $MASKDIR/combined_label_set_50um.mnc


#######################################################################################################
# Split into hemispheres

# Split 50um masks into left hemisphere only
left_mask_file_50="data/resources/left_side_50um.mnc"
for f in $MASKDIR/*_50um.mnc
do
  outfile_50=`echo $f | sed "s/.mnc/_left.mnc/g"`
  echo "Working on file: ${f} -> ${outfile_50}"
  mincmask ${f} ${left_mask_file_50} ${outfile_50} -clobber
done

# Split 200um masks into left hemisphere only
left_mask_file_200="data/resources/left_side_200um.mnc"
for f in $MASKDIR/*_200um.mnc
do
  outfile_200=`echo $f | sed "s/.mnc/_left.mnc/g"`
  echo "Working on file: ${f} -> ${outfile_200}"
  mincmask ${f} ${left_mask_file_200} ${outfile_200} -clobber
done

# Split 25um masks into left hemisphere only
left_mask_file_25="data/resources/left_side_25um.mnc"
for f in $MASKDIR/*_25um.mnc
do
  outfile_25=`echo $f | sed "s/.mnc/_left.mnc/g"`
  echo "Working on file: ${f} -> ${outfile_25}"
  mincmask ${f} ${left_mask_file_25} ${outfile_25} -clobber
done



# Split 50um masks into right hemisphere only
right_mask_file_50="data/resources/right_side_50um.mnc"
for f in $MASKDIR/*_50um.mnc
do
  outfile_50=`echo $f | sed "s/.mnc/_right.mnc/g"`
  echo "Working on file: ${f} -> ${outfile_50}"
  mincmask ${f} ${right_mask_file_50} ${outfile_50} -clobber
done

# Split 200um masks into right hemisphere only
right_mask_file_200="data/resources/right_side_200um.mnc"
for f in $MASKDIR/*_200um.mnc
do
  outfile_200=`echo $f | sed "s/.mnc/_right.mnc/g"`
  echo "Working on file: ${f} -> ${outfile_200}"
  mincmask ${f} ${right_mask_file_200} ${outfile_200} -clobber
done

# Split 25um masks into right hemisphere only
right_mask_file_25="data/resources/right_side_25um.mnc"
for f in $MASKDIR/*_25um.mnc
do
  outfile_25=`echo $f | sed "s/.mnc/_right.mnc/g"`
  echo "Working on file: ${f} -> ${outfile_25}"
  mincmask ${f} ${right_mask_file_25} ${outfile_25} -clobber
done

#######################################################################################################
# Create definitions
cat << EOF > $MASKDIR/combined_label_set_defs.csv
Structure,right.label,left.label
rest_of_brain,1,1
rest_of_cerebral_cortex,2,2
isocortex,3,3
thalamus,4,4
EOF


#######################################################################################################
# Laplacian processing

# Get laplacian, convert to MINC, downsample to 50um and 200um, and delete nrrd / 10um MINC files
# Transforming space requires up to 30GB
wget -O $MASKDIR/laplacian_10.nrrd http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/laplacian_10.nrrd
python3 code/scripts/transform_space.py $MASKDIR/laplacian_10.nrrd $MASKDIR/laplacian_10.mnc
mincresample -2 -like data/resources/average_template_50.mnc -keep_real_range $MASKDIR/laplacian_10.mnc $MASKDIR/laplacian_50um.mnc
mincresample -2 -like data/resources/average_template_200.mnc -keep_real_range $MASKDIR/laplacian_10.mnc $MASKDIR/laplacian_200um.mnc
rm $MASKDIR/laplacian_10.mnc

# Mask (at 50um) within Isocortex, and also coronal and sagittal coverage masks 
mincmask $MASKDIR/laplacian_50um.mnc $MASKDIR/Isocortex_50um.mnc $MASKDIR/laplacian_50um_under_Isocortex_50um.mnc
mincmask $MASKDIR/laplacian_50um.mnc $MASKDIR/Isocortex_coronal_coverage_50um.mnc $MASKDIR/laplacian_50um_under_Isocortex_coronal_coverage_50um.mnc
mincmask $MASKDIR/laplacian_50um.mnc $MASKDIR/Isocortex_sagittal_coverage_50um.mnc $MASKDIR/laplacian_50um_under_Isocortex_sagittal_coverage_50um.mnc

# Mask (at 200um) within Isocortex, and also coronal and sagittal coverage masks 
mincmask $MASKDIR/laplacian_200um.mnc $MASKDIR/Isocortex_200um.mnc $MASKDIR/laplacian_200um_under_Isocortex_200um.mnc
mincmask $MASKDIR/laplacian_200um.mnc $MASKDIR/Isocortex_coronal_coverage_200um.mnc $MASKDIR/laplacian_200um_under_Isocortex_coronal_coverage_200um.mnc
mincmask $MASKDIR/laplacian_200um.mnc $MASKDIR/Isocortex_sagittal_coverage_200um.mnc $MASKDIR/laplacian_200um_under_Isocortex_sagittal_coverage_200um.mnc

# Mask left, 50um
mincmask $MASKDIR/laplacian_50um.mnc ${left_mask_file_50} $MASKDIR/laplacian_50um_left.mnc
mincmask $MASKDIR/laplacian_50um_under_Isocortex_50um.mnc ${left_mask_file_50} $MASKDIR/laplacian_50um_under_Isocortex_50um_left.mnc
mincmask $MASKDIR/laplacian_50um_under_Isocortex_coronal_coverage_50um.mnc ${left_mask_file_50} $MASKDIR/laplacian_50um_under_Isocortex_coronal_coverage_50um_left.mnc
mincmask $MASKDIR/laplacian_50um_under_Isocortex_sagittal_coverage_50um.mnc ${left_mask_file_50} $MASKDIR/laplacian_50um_under_Isocortex_sagittal_coverage_50um_left.mnc

# Mask right, 50um
mincmask $MASKDIR/laplacian_50um.mnc ${right_mask_file_50} $MASKDIR/laplacian_50um_right.mnc
mincmask $MASKDIR/laplacian_50um_under_Isocortex_50um.mnc ${right_mask_file_50} $MASKDIR/laplacian_50um_under_Isocortex_50um_right.mnc
mincmask $MASKDIR/laplacian_50um_under_Isocortex_coronal_coverage_50um.mnc ${right_mask_file_50} $MASKDIR/laplacian_50um_under_Isocortex_coronal_coverage_50um_right.mnc
mincmask $MASKDIR/laplacian_50um_under_Isocortex_sagittal_coverage_50um.mnc ${right_mask_file_50} $MASKDIR/laplacian_50um_under_Isocortex_sagittal_coverage_50um_right.mnc

# Mask left, 200um
mincmask $MASKDIR/laplacian_200um.mnc ${left_mask_file_200} $MASKDIR/laplacian_200um_left.mnc
mincmask $MASKDIR/laplacian_200um_under_Isocortex_200um.mnc ${left_mask_file_200} $MASKDIR/laplacian_200um_under_Isocortex_200um_left.mnc
mincmask $MASKDIR/laplacian_200um_under_Isocortex_coronal_coverage_200um.mnc ${left_mask_file_200} $MASKDIR/laplacian_200um_under_Isocortex_coronal_coverage_200um_left.mnc
mincmask $MASKDIR/laplacian_200um_under_Isocortex_sagittal_coverage_200um.mnc ${left_mask_file_200} $MASKDIR/laplacian_200um_under_Isocortex_sagittal_coverage_200um_left.mnc

# Mask right, 200um
mincmask $MASKDIR/laplacian_200um.mnc ${right_mask_file_200} $MASKDIR/laplacian_200um_right.mnc
mincmask $MASKDIR/laplacian_200um_under_Isocortex_200um.mnc ${right_mask_file_200} $MASKDIR/laplacian_200um_under_Isocortex_200um_right.mnc
mincmask $MASKDIR/laplacian_200um_under_Isocortex_coronal_coverage_200um.mnc ${right_mask_file_200} $MASKDIR/laplacian_200um_under_Isocortex_coronal_coverage_200um_right.mnc
mincmask $MASKDIR/laplacian_200um_under_Isocortex_sagittal_coverage_200um.mnc ${right_mask_file_200} $MASKDIR/laplacian_200um_under_Isocortex_sagittal_coverage_200um_right.mnc
