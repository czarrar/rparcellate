suppressMessages(library(niftir))
source("/mnt/nfs/psych/rparcellate/command/lib/searchlight_funs.R")

load_nifti <- function(fn) {
	img <- read.nifti.image(fn)
	hdr <- read.nifti.header(fn)
	return(list(img=img, hdr=hdr))
}

matrix_like <- function(ref, val=NA, ...) matrix(val, nrow(ref), ncol(ref), ...)
array_like <- function(ref, val=NA, ...) array(val, dim(ref), ...)

run <- function(raw_cmd, ...) {
	cmd <- sprintf(raw_cmd, ...)
	cat(cmd, "\n")
	system(cmd)
}

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
	suppressMessages(require(plyr))
	suppressMessages(require(niftir))
	suppressMessages(require(biganalytics))

	n <- length(in_func_files)
	if (n != length(in_mask_files)) stop("number of func files doesn't equal number of mask files")

	# Functionals
	cat("combining", n, "functionals\n")
	## collect the dimensions to pre-make the combined functional
	dims        <- sapply(in_func_files, function(fn) read.nifti.header(fn)$dim)
	nvols       <- dims[4,]
	funcs       <- matrix(0, prod(dims[1:3,1]), sum(nvols))
	mask_vars   <- matrix(0, prod(dims[1:3,1]), n)
	for (i in 1:n) {
		cat(i, "...", sep="")
		# read in data
		func            <- read.big.nifti4d(in_func_files[i])
        # generate mask from var
        mask_var        <- colvar(func) > 0
        mask_vars[,i]   <- mask_var
        # mask dataset
        func            <- do.mask(func, mask_var)
		# scale (mean = 0, sd = 1)
        scale_fast(func, to.copy=F)
		# save
		adj             <- sum(nvols[0:(i-1)])
		tinds           <- (1:nvols[i]) + adj
		funcs[mask_var,tinds]   <- t(func[,])/sqrt(nvols[i]-1) # variance normalize
        rm(func); gc(F,T)
	}
	cat("\n")

	# Mask
	cat("combine masks\n")
	masks.main  <- sapply(in_mask_files, read.mask)
	mask.main   <- (rowSums(masks.main) == n)
	mask.var    <- rowSums(mask_vars) == n
	mask        <- mask.main & mask.var
    
    # Mask and Resize
    cat("mask and resize\n")
    funcs[!mask,]   <- 0
    dim(funcs)  <- c(dims[1:3,1], sum(nvols))
    
	# Save
	cat("save mask\n")
	hdr <- read.nifti.header(in_mask_files[[1]])
	write.nifti(mask, hdr, outfile=out_mask_file, overwrite=T)
    
    cat("save func\n")
	hdr$dim <- c(hdr$dim, ncol(funcs))
	hdr$pixdim <- c(hdr$pixdim, 1)
	write.nifti(funcs, hdr, outfile=out_func_file, overwrite=T)

	invisible(NULL)
}

# use afni or fsl
temporally_concatenate_afni <- function(in_func_files, in_mask_files, out_func_file, out_mask_file) {
	hdr <- read.nifti.header(in_mask_files[1])

	# Center and optionally scale each functional image
	# Compute the standard deviation of each image for the mask as well
	cat("computing standard deviation\n")
	mask_sd_files <- sapply(in_func_files, function(infile) {
		outfile <- tempfile(pattern="mask_sd_", fileext=".nii.gz")
		run("3dTstat -stdev -prefix %s %s", outfile, infile)
		# TODO: also call the -mean option and save that output
		# TODO: subtract the mean and divide by the standard deviation, save this as well
		outfile
	})

	# Read in the mask files and combine
	cat("combine masks\n")
	masks <- sapply(c(in_mask_files, mask_sd_files), read.mask)
	mask <- rowSums(masks) == ncol(masks)
	write.nifti(mask, hdr, outfile=out_mask_file, overwrite=T)

	# TODO: only combine the normalized dataset

	# Combine the functionals
	cat("combining the functionals\n")
	run("3dTcat -verb -output %s %s", out_func_file, paste(in_func_files, collapse=" "))

	# Removing temporary files
	cat("removing temporary files\n")
	file.remove(mask_sd_files)

	invisible(NULL)
}
