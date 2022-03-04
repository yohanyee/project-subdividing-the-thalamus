# Functions
write_minc <- function(dat, outfile, like.file, like.vol, clobber=T) {
  buffer <- like.vol
  buffer[] <- 0
  buffer[like.vol > 0.5] <- dat
  mincWriteVolume(buffer, output.filename = outfile, like.filename = like.file, clobber=clobber)
}

write_minc_cortex <- function(dat, outfile, clobber=T) {
  write_minc(dat, outfile, like.file = mask_file_cortex, like.vol = mask_vol_cortex, clobber)
}

write_minc_thalamus <- function(dat, outfile, clobber=T) {
  write_minc(dat, outfile, like.file = mask_file_thalamus, like.vol = mask_vol_thalamus, clobber)
}