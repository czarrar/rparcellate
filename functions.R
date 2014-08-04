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

# reho function
# - does the computation for a given node
# - that gets the neighbors and loops through the nodes

#' Root mean square of time-series
#'
#' Returns the RMSE between a set of time-series
#'
#' @param ts.mat matrix of time-series (ntpts x nregions)
#' 
#' @export
#' 
#' @return vector
#'
#' @examples
#' ts.mat <- matrix(rnorm(500), 100, 5)
#' rmse(ts.mat)
rmse <- function(ts.mat) {
	mean.ts <- rowMeans(ts.mat)
	diff.sq <- sweep(ts.mat, 1, mean.ts)^2
	mse	    <- rowMeans(diff.sq)
	sqrt(mean(mse))
}

#' Compute rmse-based reho
reho.rmse <- function(ts.mat, mask, ...) {
    searchlight(rmse, ts.mat, mask, ...) 
}


# Next step is to do hierarchical clustering
# We would want to ward linkage but to actually determine link with maximal similarity?
