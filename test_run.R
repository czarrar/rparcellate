source("functions.R")
library(biganalytics)

# We'll run our parcellation on the localizer data

# Masks for each hemisphere
greyfile    <- "rois/standard_grey_25pc_2mm.nii.gz"
grey        <- load_nifti(greyfile)
grey_mask   <- grey$img
hdr         <- grey$hdr

# 
preprocdir  <- "/mnt/nfs/psych/faceMemoryMRI/analysis/preprocessing/tb9226/Localizer/run01.feat"
funcfile    <- file.path(preprocdir, "reg_standard/filtered_func_data.nii.gz")
maskfile    <- file.path(preprocdir, "reg_standard/mask.nii.gz")

func.all    <- read.big.nifti(funcfile) # list$img, list$hdr
mask.sd     <- colsd(func.all) > 0
mask.func   <- read.mask(maskfile)

mask        <- (grey_mask!=0) & mask.func & mask.sd
mask.hemi   <- split_hemispheres(mask, hdr)

# scale and mask
func.lh     <- scale(func.all[,mask.hemi$lh])
func.rh     <- scale(func.all[,mask.hemi$rh])
rm(func.all)


## 2. REHO (just for lh) ######################

func            <- func.lh
tfunc           <- t(func)
mask            <- mask.hemi$lh
mask3d          <- mask
dim(mask3d)     <- hdr$dim

mask_file       <- "step00_mask.nii.gz"
write.nifti(as.numeric(mask3d), hdr, outfile=mask_file)

#func.rmse       <- reho.rmse(func, mask3d)
#write.nifti(func.rmse, hdr, mask, outfile=reho_file)
#mask_file       <- "step01_reho_mask.nii.gz"
#write.nifti(func.rmse!=0, hdr, mask, outfile=mask_file)

reho_file       <- "step01_reho.nii.gz"
raw_cmd         <- "3dReHo -mask %s -inset %s -prefix %s"
cmd             <- sprintf(raw_cmd, mask_file, funcfile, reho_file)
system(cmd)

# Smooth Results (by 2mm)
sm_file     <- "step01_reho_smooth02mm.nii.gz"
raw_cmd     <- "3dBlurInMask -input %s -FWHM 2 -mask %s -prefix %s"
cmd         <- sprintf(raw_cmd, reho_file, mask_file, sm_file)
cat(cmd, "\n")
system(cmd)
reho_sm     <- read.nifti.image(sm_file)

# Regenerate the mask
new_mask    <- reho_sm!=0
mask_file   <- "step01_reho_mask.nii.gz"
write.nifti(as.numeric(new_mask), hdr, outfile=mask_file)

old_mask3d  <- mask3d
old_mask    <- as.vector(old_mask3d)
mask3d      <- new_mask
mask        <- as.vector(new_mask)

drop_inds   <- which((old_mask[old_mask] - mask[old_mask]) > 0)
if (length(drop_inds) > 0) {
  func <- func[,-drop_inds]
  tfunc <- tfunc[-drop_inds,]
}

#' ### Peak Detection
#' I have used a minimum distance of 8mm or 4 voxels between peaks.
#+ seeds-peak
peak_file1  <- "step01_reho_smooth02mm_peaks.nii.gz"
peak_file2  <- "step01_reho_smooth02mm_peaks.txt"
raw_cmd     <- "3dExtrema -maxima -volume -closure -mask_file %s -output %s %s 1> %s"
cmd         <- sprintf(raw_cmd, mask_file, peak_file1, sm_file, peak_file2)
cat(cmd, "\n")
system(cmd)

#' Grab those peaks
peaks_img   <- read.nifti.image(peak_file1)
peak_inds   <- which(peaks_img == 1)
peak_mask_inds <- which(peaks_img[mask] == 1)


## 3. K-Means #############

# Use the identified peaks to create cluster centers
# which are the average between the center node and
# the neighoring voxels
centers <- searchlight(rowMeans, func, mask3d, process_nodes=peak_inds, percent_neighbors=0)
centers <- tfunc[peak_mask_inds,]

# Now we can compute the k-means clustering
km1_file <- "step02_kmeans.nii.gz"
system.time(km1 <- kmeans(tfunc, centers, nstart=20, iter.max=200)) # don't think i need nstart
write.nifti(km1$cluster, hdr, mask, outfile=km1_file, overwrite=T)


## 4. Spatially Constrain #############

cl_file     <- "step02_kmeans_spatclust.nii.gz"
raw_cmd     <- "3dclust -savemask %s -1noneg -isovalue -dxyz=1 0 -1 %s &> /dev/null"
cmd         <- sprintf(raw_cmd, cl_file, km1_file)
cat(cmd, "\n")
system(cmd)

# find the clusters with size less than 10
# then assign them to the nearest cluster
# now with those new clusters defined, try the k-means again started from that point with a spatial constraint
km1.clust <- read.nifti.image(cl_file)[mask]
tab <- table(km1.clust)
