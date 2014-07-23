roidir <- "~/Dropbox/Research/yale/rparcellate"

load_roi <- function(roifile) {
	img <- read.nifti.image(roifile)
	hdr <- read.nifti.header(roifile)
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

	return(list(lh=lh_roi, rh=rh_roi))
}

# reho function
# - does the computation for a given node
# - that gets the neighbors and loops through the nodes

# unconstrained k-means

# enforce spatial contiguity
# and reassign lost nodes

# final k-means with spatial enforcement
