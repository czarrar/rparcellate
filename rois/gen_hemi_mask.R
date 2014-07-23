#!/usr/bin/env Rscript

# This script will create a mask for the left and right hemisphere, separately.
suppressMessages(library(niftir))

# Load the data
cat("Loading Data\n")
stdfile <- file.path(Sys.getenv("FSLDIR"), "data", "standard", "MNI152_T1_2mm.nii.gz")
hdr <- read.nifti.header(stdfile)
img <- read.nifti.image(stdfile)
img[,,] <- 0

# Get the left-hemisphere
cat("Left Hemisphere\n")
lh <- img[,,]
lh[47:dim(lh)[1],,] <- 1
write.nifti(lh, hdr, outfile="lh_mask_2mm.nii.gz")

# Get the right-hemisphere
cat("Right Hemisphere\n")
rh <- img[,,]
rh[1:45,,] <- 1
write.nifti(rh, hdr, outfile="rh_mask_2mm.nii.gz")
