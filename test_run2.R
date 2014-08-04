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

preprocdir  <- "/mnt/nfs/psych/faceMemoryMRI/analysis/preprocessing/tb9226/Localizer/run01.feat"
funcfile    <- file.path(preprocdir, "reg_standard/filtered_func_data.nii.gz")
maskfile    <- file.path(preprocdir, "reg_standard/mask.nii.gz")

if (!file.exists(funcfile)) {
    regdir      <- file.path(preprocdir, "reg")
    applywarp(
        infile=file.path(preprocdir, "filtered_func_data.nii.gz"), 
        reffile=file.path(regdir, "standard.nii.gz"), 
        outfile=funcfile, 
        warpfile=file.path(regdir, "example_func2standard_warp.nii.gz"), 
        interp="trilinear"
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

parcels <- region_growing_wrapper(funcfile, maskfile, priorfile, outdir="ztest")
