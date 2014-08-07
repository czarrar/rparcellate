source("region_growing.R")

applywarp <- function(infile, reffile, outfile, warpfile, interp="trilinear", verbose=TRUE) {
    interp_opts <- c("nn", "trilinear", "sinc", "spline")
    if (!(interp %in% interp_opts)) stop("interp '", interp, "' not found")
    
    raw_cmd     <- "applywarp -i %s -r %s -o %s -w %s --interp=%s"
    cmd         <- sprintf(raw_cmd, 
                            infile, 
                            reffile, 
                            outfile, 
                            warpfile, 
                            interp)
    if (verbose) cat(cmd, "\n")
    system(cmd)
}

# Masks for each hemisphere
priorfile   <- "anatomical_priors/fsl_maxprob25.nii.gz"

basedir     <- "/mnt/nfs/psych/faceMemoryMRI/analysis/preprocessing"
subjects    <- list.files(basedir)

for (sub in subjects) {
    preprocdir <- file.path(basedir, sub, "Localizer/run01.feat")
    funcfile    <- file.path(preprocdir, "reg_standard/filtered_func_data.nii.gz")
    maskfile    <- file.path(preprocdir, "reg_standard/mask.nii.gz")

    if (!file.exists(funcfile)) {
        regdir      <- file.path(preprocdir, "reg")
        applywarp(
            infile=file.path(preprocdir, "filtered_func_data.nii.gz"), 
            reffile=file.path(regdir, "standard.nii.gz"), 
            outfile=funcfile, 
            warpfile=file.path(regdir, "example_func2standard_warp.nii.gz"), 
            interp="spline"
        )
    }

    if (!file.exists(maskfile)) {
        regdir      <- file.path(preprocdir, "reg")
        applywarp(
            infile=file.path(preprocdir, "mask.nii.gz"), 
            reffile=file.path(regdir, "standard.nii.gz"), 
            outfile=maskfile, 
            warpfile=file.path(regdir, "example_func2standard_warp.nii.gz"), 
            interp="nn"
        )
    }
}

featdirs  <- file.path(basedir, subjects, "Localizer/run01.feat")
subjects  <- subjects[file.exists(featdirs)]
featdirs  <- featdirs[file.exists(featdirs)]
funcfiles <- file.path(featdirs, "reg_standard/filtered_func_data.nii.gz")
maskfiles <- file.path(featdirs, "reg_standard/mask.nii.gz")

outdir <- "ztest5"; dir.create(outdir)
mask_file <- file.path(outdir, "concatenate_mask.nii.gz")
func_file <- file.path(outdir, "concatenate_func.nii.gz")
temporally_concatenate(funcfiles, maskfiles, func_file, mask_file)

outdir2 <- file.path("ztest5/region_growing")
parcels <- region_growing_wrapper(func_file, mask_file, priorfile, outdir=outdir2, roi.scale=10000)


###
# Subject-Specific
###

# utilities
hdr      <- read.nifti.header(file.path(outdir2, "mask.nii.gz"))
mask     <- read.mask(file.path(outdir2, "mask.nii.gz"))
mparcels <- parcels[mask]
uparcels <- sort(unique(mparcels))
nparcels <- length(uparcels)
rparcels <- vector("integer", length(mparcels))
for (i in 1:nparcels) rparcels[mparcels == uparcels[i]] <- i

# now let's get subject specific data (here one subject)
func    <- read.big.nifti4d(funcfiles[1])
func    <- do.mask(func, mask)
func    <- scale_fast(func, to.copy=F)
func    <- func[,]/sqrt(nrow(func)-1)

# compute the mean ts of each parcel for this subject
# note this is nparcels x ntpts, whereas func is ntpts x nvoxs
parcel_ts <- laply(1:nparcels, function(i) {
    rowMeans(func[,mparcels==uparcels[i]])
}, .progress="text")

# now compute spatial map associated with each mean parcel ts
parcel_maps <- parcel_ts %*% func

# with these maps, we want to know what where each cluster belongs
# calculate the closest parcellation for each voxel
parcel_max_cor <- aaply(parcel_maps, 2, which.max, .progress="text")
# then see voxels that are the same as in the prior parcellation
parcel_same    <- parcel_max_cor == rparcels

# for voxels that are not the same, then add them to a region if in the same cluster (spatially contigious)
# basically for all 
# thinking if i actually do some cluster
write.nifti(parcel_max_cor, hdr, mask, outfile="tmp_max.nii.gz")
cmd <- "3dclust -1noneg -isovalue -prefix tmp_clust.nii.gz -dxyz=1 0 1 tmp_max.nii.gz 1> tmp_max.txt"
cat(cmd, "\n"); system(cmd)
parcel_max_cor_clust <- read.nifti.image("tmp_clust.nii.gz")[mask]