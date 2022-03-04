#!/usr/bin/Rscript

# Load libraries
library(tidyverse)
library(glue)
library(RMINC)
library(RcppCNPy)
library(psych)
library(data.table)

###########################
# Arguments
###########################

# Testing arguments
# outdir <- "data/outputs/isocortex_thalamus_left"

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]

# Input files
mask_file_thalamus <- glue("{outdir}/thalamus_mask_highres.mnc")
mask_file_cortex <- glue("{outdir}/cortex_mask_lowres.mnc")
coexpression_file <- glue("{outdir}/coexpression.npy")
distance_file <- glue("{outdir}/euclidean_distance.npy")
sc_file <- glue("{outdir}/correlation_z.npy")

###########################
# Read data
###########################

print(glue("[{Sys.time()}] Reading mask data"))

mask_vol_thalamus <- mincGetVolume(mask_file_thalamus)
mask_vol_cortex <- mincGetVolume(mask_file_cortex)

p_thal <- length(mask_vol_thalamus[mask_vol_thalamus > 0.5])
p_cortex <- length(mask_vol_cortex[mask_vol_cortex > 0.5])

# Takes a 3-5 minutes to read in all the data
print(glue("[{Sys.time()}] Reading SC (z) matrix"))
zmat <- npyLoad(sc_file)
print(glue("[{Sys.time()}] Reading transcriptomic similarity matrix"))
cmat <- npyLoad(coexpression_file)
print(glue("[{Sys.time()}] Reading Euclidean distance matrix"))
dmat <- npyLoad(distance_file)

###########################
# Compute residuals
###########################

print(glue("[{Sys.time()}] Computing residuals"))
# Should take 10-15 minutes in total

resmat <- matrix(NA, nrow=p_thal, ncol=p_cortex)
df_intercept <- tibble(term="Intercept", i=numeric(p_thal), estimate=numeric(p_thal), statistic=numeric(p_thal), pvalue=numeric(p_thal))
df_dist <- tibble(term="Distance", i=numeric(p_thal), estimate=numeric(p_thal), statistic=numeric(p_thal), pvalue=numeric(p_thal))
df_coexp <- tibble(term="Coexpression", i=numeric(p_thal), estimate=numeric(p_thal), statistic=numeric(p_thal), pvalue=numeric(p_thal))
df_model <- tibble(term="Model", i=numeric(p_thal), adjusted_rsq=numeric(p_thal), aic=numeric(p_thal), bic=numeric(p_thal))

# Clean up before starting
gc()

# Start
prog <- txtProgressBar(max=p_thal, style=3)
for (i in 1:p_thal) {
  
  df <- tibble(z=zmat[i,], coexp=fisherz(cmat[i,]), d=dmat[i,])
  fit <- lm(z ~ d + coexp, df)
  sc_res <- residuals(fit)
  
  # Testing
  #summary(fit)
  #hist(sc_res)
  #qplot(sc_res, df$z) + geom_abline(slope=1, intercept = 0, lty=3, color='red') + theme_bw()
  #plot(fit)
  #ggplot(df, aes(x=coexp, y=z)) + geom_point(aes(color=d)) + geom_smooth(method=lm) + theme_bw()
  #ggplot(df, aes(x=d, y=z)) + geom_point(aes(color=coexp)) + geom_smooth(method=lm)  + theme_bw()

  resmat[i,] <- sc_res
  
  sfit <- summary(fit)
  
  df_intercept$i[[i]] <- i
  df_intercept$estimate[[i]] <- sfit$coefficients["(Intercept)", "Estimate"]
  df_intercept$statistic[[i]] <- sfit$coefficients["(Intercept)", "t value"]
  df_intercept$pvalue[[i]] <- sfit$coefficients["(Intercept)", "Pr(>|t|)"]
  
  df_dist$i[[i]] <- i
  df_dist$estimate[[i]] <- sfit$coefficients["d", "Estimate"]
  df_dist$statistic[[i]] <- sfit$coefficients["d", "t value"]
  df_dist$pvalue[[i]] <- sfit$coefficients["d", "Pr(>|t|)"]
  
  df_coexp$i[[i]] <- i
  df_coexp$estimate[[i]] <- sfit$coefficients["coexp", "Estimate"]
  df_coexp$statistic[[i]] <- sfit$coefficients["coexp", "t value"]
  df_coexp$pvalue[[i]] <- sfit$coefficients["coexp", "Pr(>|t|)"]
  
  df_model$i[[i]] <- i
  df_model$adjusted_rsq[[i]] <- sfit$adj.r.squared
  df_model$aic[[i]] <- AIC(fit)
  df_model$bic[[i]] <- BIC(fit)
  
  setTxtProgressBar(prog, i)
}
close(prog)

# Garbage collection
print(glue("[{Sys.time()}] Removing unneeded objects"))
rm(list = c("zmat", "cmat", "dmat"))
gc()

# Examine a few things
print(glue("[{Sys.time()}] Some stats..."))
dim(resmat)
resmat[1:5, 1:5]

print(glue("[{Sys.time()}] Minimum thalamocortical z-transformed residual correlation is {round(min(resmat), 5)}"))
print(glue("[{Sys.time()}] Minimum thalamocortical z-transformed residual correlation is {round(max(resmat), 5)}"))

# Bind dataframes together
print(glue("[{Sys.time()}] Binding dataframes together"))
df_params <- rbindlist(l = list(df_intercept, df_dist, df_coexp))

###########################
# Write out data
###########################

# Write out SC (residual) matrix
print(glue("[{Sys.time()}] Writing out data"))
npySave(filename = glue("{outdir}/correlation_residual.npy"), object = resmat)

# Write out model dataframes
print(glue("[{Sys.time()}] Writing out dataframes"))
write_csv(df_params, path = glue("{outdir}/correlation_regression_model_parameters.csv"))
write_csv(df_model, path = glue("{outdir}/correlation_regression_model_fit.csv"))

# Function
write_minc <- function(dat, outfile, clobber=T) {
  buffer <- mask_vol_thalamus
  buffer[] <- 0
  buffer[mask_vol_thalamus > 0.5] <- dat
  mincWriteVolume(buffer, output.filename = outfile, like.filename = mask_file_thalamus, clobber=clobber)
}

# Output model estimates (intercept, dist, coexp) as minc files
write_minc(df_intercept$estimate, outfile = glue("{outdir}/correlation_regression_model_beta_intercept.mnc"))
write_minc(df_dist$estimate, outfile = glue("{outdir}/correlation_regression_model_beta_distance.mnc"))
write_minc(df_coexp$estimate, outfile = glue("{outdir}/correlation_regression_model_beta_coexp.mnc"))

# Output model fit as minc file (adj R squared)
write_minc(df_model$adjusted_rsq, outfile = glue("{outdir}/correlation_regression_model_adjusted_rsq.mnc"))

print(glue("[{Sys.time()}] Done!"))

# Test that the residual map is OK
for (k in c(1,2,3,10,100,200,300,1000, 15000, 15001, 15002, 25000, 40000, 50000, 50005, 50050, 50500, 55000, 75000)) {
  print(glue("Working on k = {k}"))
  ovol <- mask_vol_cortex
  ovol[] <- 0
  ovol[mask_vol_cortex > 0.5] <- resmat[k,]
  min(resmat[k,])
  max(resmat[k,])
  thalamus_mask_indices <- which(mask_vol_thalamus > 0.5)
  voxel_coords_test <- mincVectorToVoxelCoordinates(mask_file_thalamus, thalamus_mask_indices[k])
  world_coords_test <- mincConvertVoxelToWorld(mask_file_thalamus, voxel_coords_test[1], voxel_coords_test[2], voxel_coords_test[3])
  mincWriteVolume(ovol, glue("tmp/SC_res_test_{k}_at_world_{world_coords_test[1]}_{world_coords_test[2]}_{world_coords_test[3]}.mnc"))
}
