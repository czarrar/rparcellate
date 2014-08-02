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


# region growing
# - start with some regions (assignment and time-courses)
# - compute correlation between each region and all the neighbors
# - join if the correlation is more than 90% of the maximum correlation btw region + neis
# - iterate
region_growing <- function(func, start_nodes, mask3d) {
	## Setup ##
    mask        <- as.vector(mask3d)
	nnodes      <- ncol(func)
	ntpts		<- nrow(func)
	nregions    <- length(start_nodes)
	regions	    <- vector("numeric", nnodes)

	## Initialize region time-course ##
	# Average within a radius of 1.5 times the voxel size
	start_neis 	<- find_neighbors_masked(mask3d, start_nodes, 
										 nei=1, nei.dist=1.5, 
										 verbose=FALSE) # add this!
	for (ri in 1:nregions) {
		neis <- start_neis[[ri]]
        regions[neis] <- ri
	}
	
	# Iterate until no voxel is unassigned
	nleft <- sum(regions == 0)
	cat("# left", nleft, "\n")
	while (nleft > 0) {
		new_regions <- regions[]
        
        # calculate region time-series
        region_ts <- sapply(1:nregions, function(ri) {
            rowMeans(func[,new_regions==ri,drop=F])
        })
	    
		# Loop through each region
		for (ri in 1:nregions) {
			# Get region ts
			rts <- region_ts[,ri]
			
			# Get neighborhood ts
			region_nodes	<- which(regions == ri)
			lst_neis 		<- find_neighbors_masked(mask3d, region_nodes, 
										include.self=FALSE, verbose=FALSE)
			nei_nodes		<- unique(unlist(lst_neis))
			nts				<- func[,nei_nodes,drop=F]
			
			# Compute correlation btw region and neighbors
			rn_cors			<- cor(rts, nts)[1,]
			
			# Select those with 90% of maximum correlation
			max_cor			<- 0.9 * max(rn_cors)
			neis_join  		<- nei_nodes[rn_cors > max_cor]
			
			# Assign nodes to region
			# but if already assigned see if our correlation is greater
			for (nei in neis_join) {
				if (new_regions[nei] == 0) {
					new_regions[nei] <- ri
				} else {
					prev_reg    <- new_regions[nei]
					prev_cor    <- cor(region_ts[,prev_reg], func[,nei])
					cur_cor 	<- rn_cors[nei_nodes == nei]
					if (prev_cor > cur_cor) new_regions[nei] <- ri
				}
			}
		}
		        
		regions <- new_regions
		nleft <- sum(regions == 0)
		cat("# left", nleft, "\n")
	}
	
	return(regions)
}


# Next step is to do hierarchical clustering
# We would want to ward linkage but to actually determine link with maximal similarity?
