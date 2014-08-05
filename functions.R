suppressMessages(library(niftir))
source("searchlight_funs.R")

roidir <- "~/Dropbox/Research/yale/rparcellate"

load_nifti <- function(fn) {
	img <- read.nifti.image(fn)
	hdr <- read.nifti.header(fn)
	return(list(img=img, hdr=hdr))
}

matrix_like <- function(ref, val=NA, ...) matrix(val, nrow(ref), ncol(ref), ...)
array_like <- function(ref, val=NA, ...) array(val, dim(ref), ...)

split_hemispheres <- function(roi, hdr) {
	# Get the center coordinate
	center.coord <- as.numeric(hdr$qto.ijk %*% c(0,0,0,1))[-4]
	cat("xyz: 0 0 0 => ijk:", center.coord, "\n")
    
  # Note the above is based on the first element being 0
  # need to convert 0 => 1
  cx <- center.coord[1] + 1 # center of x-axis
  ex <- dim(roi)[1] # end of x-axis

    # Split the hemispheres
	lh_roi <- roi[,,]
	rh_roi <- roi[,,]
	lh_roi[1:cx,,] <- 0	
	rh_roi[cx:ex,,] <- 0

	return(list(lh=as.vector(lh_roi)!=0, rh=as.vector(rh_roi)!=0))
}

# Temporally concatenate data from multiple subjects / runs
temporally_concatenate <- function(in_func_files, in_mask_files, out_func_file, out_mask_file) {
	require(plyr)

	n <- length(in_func_files)
	if (n != length(in_mask_files)) stop("number of func files doesn't equal number of mask files")

	# Mask
	cat("combine masks\n")
	masks <- sapply(in_mask_files, read.mask)
	mask <- rowSums(masks) == n

	# Functionals
	cat("combine functionals\n")
	funcs <- laply(in_func_files, function(in_func_file) {
		func <- read.big.nifti(func_file)
		func <- do.mask(func, mask)
		as.matrix(func)
	})
	# TODO: merge some of the dimensions here
	dim(funcs) <- c(dim(funcs)...)

	# Save
	cat("save\n")
	hdr <- read.nifti.header(in_func_files[[1]])
	write.nifti(mask, hdr, outfile=out_mask_file, overwrite=T)
	write.nifti(funcs, hdr, mask, outfile=out_func_file, overwrite=T)

	invisible(NULL)
}
