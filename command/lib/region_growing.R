# Playing with converting to C code
library(inline)
library(RcppArmadillo)
suppressMessages(library(biganalytics))
script.dir <- dirname(sys.frame(1)$ofile)
source(file.path(script.dir, "functions.R"))

cat("compiling\n")

code <- '
	using namespace arma;

	/**
	** Setup User Vars
	**/

    Rcpp::List list_neis(R_list_neis);
    
    Rcpp::NumericVector r_regions(R_regions);
    vec regions(r_regions.begin(), r_regions.size(), false);
    vec new_regions = regions;

    double nregions =  Rcpp::as<double>(R_nregions);
    double nnodes = (double)regions.n_elem;
    
    Rcpp::NumericMatrix r_node_ts(R_node_ts);
    mat node_ts(r_node_ts.begin(), r_node_ts.nrow(), r_node_ts.ncol(), false);

   	/**
	** Setup Other Vars
	**/

    // Degrees of freedom for when we calculate the correlation
    //double tdf = (double)node_ts.n_rows - 1;
    
    // Vector of neighbors with a given region
    vec region_neis = zeros<vec>(nnodes);
    
    // For each node, highest correlation value to neighboring region
    vec best_cors = zeros<vec>(nnodes);
    
    double max_cor; // Stores maximum correlation of region to neighbors
    arma::uword nei_ind; // Region neighbor index
    
    /**
	** Calculate
	**/

	// Iterate until no voxel is unassigned
	uword nleft = 0;
	for (uword ni=0; ni < regions.n_elem; ni++) {
		if (regions(ni) == 0) {
			nleft = nleft + 1;
		}
	}
	printf("left -> %i\\n", nleft);

	mat region_ts(node_ts.n_rows, nregions);
	double nchanges = 1; // Number of changes to node assignment made

	while (nleft > 0 && nchanges > 0) {
		nchanges = 0; // Number of changes to node assignment made

        // calculate region time-series
        //region_ts.zeros();
        for (double ri=0; ri < nregions; ri++) {
        	// row mean across region nodes
        	region_ts.col(ri) = arma::mean(node_ts.cols(find(regions == (ri+1))), 1);
        }
        region_ts = arma::normalise(region_ts);

	    for (double ri=0; ri < nregions; ri++) {
	        
	        // nodes/voxels for the given region
	        uvec region_nodes = find(regions == (ri+1));
	        
	        // all possible neighbors to nodes within region
	        region_neis.zeros();
			for (uword rni=0; rni<region_nodes.n_elem; rni++) {
	            // for each loop, get the neighbors for one particular node
				SEXP ln = list_neis[region_nodes(rni)];
				uvec node_neis = Rcpp::as<uvec>(ln) - 1;
				// only add node neis that do not belong to another region
				for (uword nni=0; nni<node_neis.n_elem; nni++) {
					if (regions[node_neis(nni)] == 0) {
						region_neis(node_neis(nni)) = 1;
					}
				}
			}
	        // get indices of region neighbors
	        uvec nei_inds = find(region_neis == 1);
	        
	        // if no neighbors left, skip
	        if (nei_inds.n_elem == 0) continue;
	        
	        // since the time-series are normalized, we do a simpler computation 
	        // to get the correlation between this region and its neighbors
	        rowvec rn_cors = region_ts.col(ri).t() * node_ts.cols(nei_inds);
	        
	        // determine correlation threshold for joining
	        max_cor = 0.9 * rn_cors.max();
	        
	        // allow a neighboring voxel to join if
	        // 1. it has a correlation greater than the threshold (`max_cor`)
	        // 2. the correlation is greater than with any other neighboring region
	        for (uword ni=0; ni<nei_inds.n_elem; ni++) {
	            nei_ind = nei_inds(ni);
	            if ((rn_cors(ni) > max_cor) && (new_regions(nei_ind) == 0 || rn_cors(ni) > best_cors(nei_ind))) {
	                new_regions(nei_ind) = ri + 1;
	                best_cors(nei_ind) = rn_cors(ni);
	                nchanges = nchanges + 1;
	            }
	        }
	    }
    	
    	nleft = 0;
    	for (uword ni=0; ni < regions.n_elem; ni++) {
    		regions(ni) = new_regions(ni);
    		if (regions(ni) == 0) {
    			nleft = nleft + 1;
    		}
    	}
		
		//printf("nleft -> %i\\n", nleft);
    	//printf("changes -> %.0f\\n", nchanges);
    }

    return Rcpp::wrap(new_regions);
'

start_region_growing_worker <- cxxfunction(
    signature(R_node_ts="numeric", R_list_neis="List", R_regions="numeric", R_nregions="numeric"), 
    code, plugin="RcppArmadillo"
)

code <- '
    using namespace arma;

    /*** INPUTS ***/

    // Functional data: voxels x time-points
    Rcpp::NumericMatrix r_node_ts(R_node_ts);
    mat node_ts(r_node_ts.begin(), r_node_ts.nrow(), r_node_ts.ncol(), false);
        
    // Current cluster/parcel assignments
    uvec regions = Rcpp::as<uvec>(R_regions);
    uvec new_regions = regions;
    
    // List of neighbors for each voxel
    Rcpp::List list_voxel_neis(R_list_voxel_neis);
    
    // Total number of voxels
    uword ntpts = node_ts.n_rows;
    uword nvoxs = regions.n_elem;
    uword nregions =  Rcpp::as<uword>(R_nregions);
    
    int maxiter = 200;
    
    
    /*** CLUSTER ***/
    
    // Stop when there are no changes in cluster assignment
    int nchanges = 1; int niter = 0;
    mat region_ts(ntpts, nregions);
    uword best_index; double best_dist; uword best_region;
    uword curr_region;
    while (nchanges > 0 && niter < maxiter) {
        nchanges = 0;
        
        // calculate region time-series
        for (uword ri=0; ri < nregions; ri++) {
        	// row mean across region nodes
        	region_ts.col(ri) = arma::mean(node_ts.cols(find(regions == (ri+1))), 1);
        }
        region_ts = arma::normalise(region_ts);
        
        // Loop through each node and check if node should change clusters
        for (uword vi=0; vi<nvoxs; vi++) {
        
            // Get the neigboring regions for given node
            SEXP ll = list_voxel_neis[vi];
            uvec neis = Rcpp::as<uvec>(ll) - 1;
            uvec nei_regions = arma::unique(regions(neis)); // unique neighboring regions
            
            // Skip if neighbors are all part of the same cluster
            if (nei_regions.n_elem == 1) {
                continue;
            }
            
            // Distance to neighboring cluster centers
            // and choose the best cluster with smallest distance to voxel
            // note: current nodes clust should be in the mix here
	        
            // since the time-series are normalized, we do a simpler computation 
	        // to get the correlation between this region and its neighbors
	        rowvec rn_cors = node_ts.col(vi).t() * region_ts.cols(nei_regions - 1);
            
            // pick the region with the largest correlation
            best_dist = rn_cors.max(best_index);
            best_region = nei_regions(best_index);
            
            // If the best cluster is the same, move on
            curr_region = regions(vi);
            if (best_region == curr_region) {
                continue;
            }
            
            // Otherwise, make the change
            nchanges = nchanges + 1;
            new_regions(vi) = best_region;
        }
        
        // copy
    	for (uword vi=0; vi < nvoxs; vi++) {
    		regions(vi) = new_regions(vi);
    	}
        
        niter = niter + 1;
        //printf("%03i: changes -> %i\\n", niter, nchanges);
    }
    
    if (niter >= maxiter) {
        printf("warning: did not converge. %i changes in the last iteration.", nchanges);
    }
    
    return Rcpp::wrap(regions);
'

refine_region_growing_worker <- cxxfunction(
    signature(R_node_ts="numeric", R_list_voxel_neis="List", R_regions="integer", R_nregions="integer"), 
    code, plugin="RcppArmadillo"
)

cat("done compiling\n")


#' func should be scaled
region_growing <- function(func, starting_nodes, mask, hdr, normalize=T) {
    cat("Setup\n")
    mask3d      <- mask
    dim(mask3d) <- hdr$dim
	nnodes      <- ncol(func)
	ntpts		<- nrow(func)
	nregions    <- as.double(length(starting_nodes))
	regions	    <- vector("numeric", nnodes)
    if (normalize) func <- func/sqrt(ntpts-1)
    
    
    cat("Initialize regions\n")
	# Nodes with surrounding nodes within a radius of 1.5 times the voxel size
	start_neis 	<- find_neighbors_masked(mask3d, starting_nodes, 
										 nei=1, nei.dist=1.5, 
										 verbose=FALSE)
	for (ri in 1:nregions) {
		neis <- start_neis[[ri]]
        regions[neis] <- ri
	}
    regions <- as.double(regions)
    
    
    cat("Compile neighbors\n")
    # for each node get nearest neighbors (exclude self)
	lst_neis <- find_neighbors_masked(mask3d, include.self=F, verbose=F)
    
    
    cat("Grow regions\n")
    system.time(regions2 <- start_region_growing_worker(func, lst_neis, regions, nregions))
    
    
    # can have some isolates...remove those
    isolates <- regions2==0
    if (any(isolates)) {
        cat("Removing isolates:", sum(isolates), "\n")
        regions2 <- regions2[!isolates]
        mask3d[mask][isolates] <- F
        func <- func[,!isolates]
        lst_neis <- find_neighbors_masked(mask3d, include.self=F, verbose=F)
    }
    
    
    cat("Refine regions\n")
    regions2 <- as.integer(regions2)
    nregions <- as.integer(length(unique(regions2)))
    system.time(regions3 <- refine_region_growing_worker(func, lst_neis, regions2, nregions))
    
    # DONE!
    list(regions=as.vector(regions3), mask=as.vector(mask3d), regions0=as.vector(regions2), 
    	nregions=nregions)
}

reho_peak_detection <- function(func, mask, hdr, fwhm=2, outprefix=NULL) {
    #  outprefix NULL then create a temp directory
    
    mask3d          <- mask
    dim(mask3d)     <- hdr$dim
    
    if (is.null(outprefix)) {
    	cat("creating temporary directory\n")
    	tmpdir <- tempdir()
    	outprefix <- file.path(tmpdir, "reho")
    } else {
    	tmpdir <- NULL
    }

    cat("running rmse reho\n")
    reho_file       <- paste(outprefix, "map.nii.gz", sep="_")
    func.rmse       <- reho.rmse(func, mask3d, percent_neighbors=1)
    write.nifti(func.rmse, hdr, mask, outfile=reho_file) # don't have this allow anything

    # new mask since not every voxel will have an output
    rmask_file      <- paste(outprefix, "mask.nii.gz", sep="_")
    write.nifti(func.rmse!=0, hdr, mask, outfile=rmask_file)
    
    # Smooth Results (by 2mm)
    cat("smooth results\n")
    sm_file     <- paste(outprefix, "smooth02mm.nii.gz", sep="_")
    raw_cmd     <- "3dBlurInMask -input %s -FWHM 2 -mask %s -prefix %s"
    cmd         <- sprintf(raw_cmd, reho_file, rmask_file, sm_file)
    cat(cmd, "\n")
    system(cmd)
    reho_sm     <- read.nifti.image(sm_file)
    
    #' ### Peak Detection
    cat("peak detection\n")
    #' I have used a minimum distance of 8mm or 4 voxels between peaks.
    #+ seeds-peak
    peak_file1  <- paste(outprefix, "smooth02mm_peaks.nii.gz", sep="_")
    peak_file2  <- paste(outprefix, "smooth02mm_peaks.txt", sep="_")
    raw_cmd     <- "3dExtrema -minima -volume -closure -mask_file %s -output %s %s 1> %s"
    cmd         <- sprintf(raw_cmd, rmask_file, peak_file1, sm_file, peak_file2)
    cat(cmd, "\n")
    system(cmd)

    #' Grab those peaks
    peaks_img   <- read.nifti.image(peak_file1)
    peak_inds   <- which(peaks_img == 1)
    peak_mask_inds <- which(peaks_img[mask] == 1)
    
    cat("Found", length(peak_mask_inds), "peaks based on REHO-like analysis\n")
    
    if (!is.null(tmpdir)) {
    	file.remove(list.files(tmpdir))
    	file.remove(tmpdir)
    }

    list(img=peaks_img, all.inds=peak_inds, mask.inds=peak_mask_inds)
}

region_growing_wrapper <- function(func_file, mask_file, roi_file, outdir, roi.scale=10000) {
	outfile <- file.path(outdir, "parcels.nii.gz")
    if (file.exists(outfile)) stop("output outdir/parcels.nii.gz cannot exist")

	cat("Read in data\n")
    hdr         <- read.nifti.header(mask_file)
    ## assume that mask can be multiple different ROIs
    mask.main   <- read.mask(mask_file)
    rois.all    <- as.vector(read.nifti.image(roi_file))
    rois        <- unique(rois.all[rois.all!=0])
    nrois       <- length(rois)
    nvoxs		<- length(mask.main)
    ## read in functional image
    func.all    <- read.big.nifti(func_file)
    mask.sd     <- colsd(func.all) > 0

    # Get regions in each mask/roi
    cat("Will process", nrois, "mask-rois\n")
    parcels.rois <- sapply(1:nrois, function(i) {
        cat("...", i, "\n", sep="")
        
        mask    <- mask.main & mask.sd & (rois.all == rois[i])
        func    <- scale(func.all[,mask])
        
        outprefix <- ifelse(is.null(outdir), NULL, sprintf("%s/roi%02i_reho", outdir, i))
        peaks   <- reho_peak_detection(func, mask, hdr, fwhm=2, outprefix=outprefix)
        ret     <- region_growing(func, peaks$mask.inds, mask, hdr)
        
        if (ret$nregions > roi.scale) warning("nregions > scale: ", ret$nregions, " vs ", roi.scale)
        parcels <- vector("integer", nvoxs)
        parcels[ret$mask] <- ret$regions + i*roi.scale
        
    	cat("====\n")
        parcels
    })
    
    # Combine
    cat("Combine\n")
    parcels     <- rowSums(parcels.rois)

    cat("Relabel\n")
    regions     <- sort(unique(parcels[parcels!=0]))
    nregions    <- length(regions)
    reparcel	<- vector("numeric", nvoxs)
    for (i in 1:nregions) {
    	inds <- parcels == regions[i]
    	reparcel[inds] <- i
    }

    # if outdir is given, then save the parcels as an image as well
    if (!is.null(outdir)) {
    	cat("Save output\n")
    	mask <- mask.main & mask.sd & (rois.all>0)
    	mask3d <- mask; dim(mask3d) <- hdr$dim
    	write.nifti(mask3d*1, hdr, outfile=file.path(outdir, "mask.nii.gz"))
    	write.nifti(parcels, hdr, outfile=file.path(outdir, "parcels.nii.gz"))
    	write.nifti(reparcel, hdr, outfile=file.path(outdir, "parcels_relabel.nii.gz"))
    }

    invisible(parcels)
}

region_growing_group <- function(func_files, mask_files, roi_file, outdir) {
	cat("Setup group output\n")
	if (file.exists(outdir)) stop("output directory cannot exist")
	dir.create(outdir)

	cat("Temporally concatenate\n")
	mask_file <- file.path(outdir, "concatenate_mask.nii.gz")
	func_file <- file.path(outdir, "concatenate_func.nii.gz")
	temporally_concatenate(func_files, mask_files, func_file, mask_file)

	cat("Group-level parcellation\n")
	region_growing_outdir <- file.path(outdir, "region_growing")
	grp.parcels <- region_growing_wrapper(func_file, mask_file, roi_file, region_growing_outdir)

	cat("Subject-level parcellation\n")
	# for each subject, get mean time-series for each region
	# compute spatial correlation of this time-series with the rest of the brain
	# keep only those voxels that are spatially connected with main cluster
	# and have max correlation greater than any other parcellation
}

# one approach to subject-level is dual regression style
# for each partition, get mean time-series and then compute correlation with everything
# pick max?

