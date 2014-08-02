# Playing with converting to C code
library(inline)
library(RcppArmadillo)

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
		
		printf("left -> %i\\n", nleft);
    	printf("changes -> %.0f\\n", nchanges);
    }

    return Rcpp::wrap(new_regions);
'

test_quick_region_grow_worker <- cxxfunction(
    signature(R_node_ts="numeric", R_list_neis="List", R_regions="numeric", R_nregions="numeric"), 
    code, plugin="RcppArmadillo"
)

test_quick_region_grow <- function(func, start_nodes, mask3d) {
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
    regions <- as.double(regions)
        
    # calculate neighbors for each node
	lst_neis <- find_neighbors_masked(mask3d, include.self=F, verbose=F)
    
    # loop through each region and grow
    print(system.time(ret <- test_quick_region_grow_worker(func/sqrt(ntpts-1), lst_neis, regions, as.double(nregions))))
    
    return(ret)
}

system.time(comp <- test_quick_region_grow(func, peak_mask_inds, mask3d))


test_slow_region_grow <- function(func, start_nodes, mask3d, to.test=FALSE) {
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
    regions <- as.double(regions)
    new_regions <- regions[]
    
    # calculate neighbors for each node
	lst_neis <- find_neighbors_masked(mask3d, include.self=F, verbose=F)
    
    # Iterate until no voxel is unassigned
    nchanges <- 1
	nleft <- sum(regions == 0)
	cat("left =>", nleft, "\n")
	while (nleft > 0 && nchanges > 0) {
       	nchanges <- 0; nchanges1 <- 0; nchanges2 <- 0
	    best_cors <- vector("numeric", nnodes)

        # calculate region time-series
        region_ts <- sapply(1:nregions, function(ri) {
            rowMeans(func[,regions==ri,drop=F])
        })
        region_ts <- scale(region_ts)
	    
        # Loop through each region
		for (ri in 1:nregions) {		
			# Get neighborhood
			region_nodes	<- which(regions == ri)
			nei_nodes		<- unique(unlist(lst_neis[region_nodes]))
			## Exclude nodes that are already assigned to a region
	        #nei_nodes       <- nei_nodes[!(nei_nodes %in% region_nodes)]
	        nei_nodes		<- nei_nodes[regions[nei_nodes] == 0]
	        ## Exclude other
	        if (length(nei_nodes) == 0) next
	        
	        # Get ts
	        rts             <- region_ts[,ri]
			nts				<- func[,nei_nodes,drop=F]
			
			# Compute correlation btw region and neighbors
			#rn_cors			<- cor(rts, nts)[1,]
	        rn_cors			<- (t(rts) %*% nts)[1,]/(nrow(func)-1)
	        
			# Select those with 90% of maximum correlation
			max_cor			<- 0.9 * max(rn_cors)
			neis_join  		<- nei_nodes[rn_cors > max_cor]
	        
			# Assign nodes to region
			# but if already assigned see if our correlation is greater
			for (nei in neis_join) {
				cur_cor 	<- rn_cors[nei_nodes == nei]
				if (new_regions[nei] == 0) {
					new_regions[nei] <- ri
	               	nleft <- nleft - 1
	               	nchanges1 <- nchanges1 + 1
	                best_cors[nei] <- cur_cor
				} else if (cur_cor > best_cors[nei]) {
                    new_regions[nei] <- ri
                    nchanges2 <- nchanges2 + 1
                    best_cors[nei] <- cur_cor
				}
			}
		}

		regions <- new_regions
		nchanges <- nchanges1 + nchanges2
		cat(sprintf("changes => %i (%i,%i)\n", nchanges, nchanges1, nchanges2))
		cat("left =>", nleft, "\n")

		if (to.test) break
    }

    
    return(regions)
}


system.time(comp <- test_quick_region_grow(func, peak_mask_inds, mask3d))
system.time(ref  <- test_slow_region_grow(func, peak_mask_inds, mask3d))

# find the mismatches
tmp0 <- sapply(1:max(ref), function(i) all.equal(as.vector(comp)==i, ref==i))
## i might then focus on each to figure out the deal

system.time(comp2 <- test_quick_region_grow(func, peak_mask_inds, mask3d))
system.time(ref2  <- test_slow_region_grow(func, peak_mask_inds, mask3d))

system.time(ref3  <- test_slow_region_grow(func, peak_mask_inds, mask3d))






system.time(comp <- test_quick_region_grow(func, peak_mask_inds, mask3d))
system.time(ref  <- test_quick_region_ts(func, peak_mask_inds, mask3d))


system.time(comp <- test_quick_region_grow(func, peak_mask_inds, mask3d))


system.time(ref  <- test_slow_region_grow(func, peak_mask_inds, mask3d))
write.nifti(ref, hdr, mask, outfile="step02_region_growing_r.nii.gz", overwrite=T)



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
            nei_nodes       <- nei_nodes[!(nei_nodes %in% region_nodes)]
			nts				<- func[,nei_nodes,drop=F]
			
			# Compute correlation btw region and neighbors
			rn_cors			<- cor(rts, nts)[1,]
			#rn_cors			<- (t(rts) %*% nts)[1,]/(nrow(func)-1)
			
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




# Let's see if we can modify and return with reference
code <- '
    arma::vec regions = Rcpp::as<arma::vec>(R_regions);
    regions(0) = 100;
    return Rcpp::wrap(regions);
'
test_modify_vector <- cxxfunction(
    signature(R_regions="numeric"), 
    code, plugin="RcppArmadillo"
)
regions <- vector("numeric", 10)
cp_regions <- test_modify_vector(as.double(regions))


code <- '
    Rcpp::NumericVector r_regions(R_regions);
    arma::vec regions(r_regions.begin(), r_regions.size(), false);
    regions(0) = 100;
    arma::vec new_regions = regions;
    new_regions(0) = 200;
    return Rcpp::wrap(new_regions);
'
test_modify_vector <- cxxfunction(
    signature(R_regions="numeric"), 
    code, plugin="RcppArmadillo"
)
regions <- vector("numeric", 10)
cp_regions <- test_modify_vector(regions)
all.equal(as.numeric(regions), as.numeric(cp_regions))










# This isn't really region growing
# Instead it is a spatially constrained k-means algorithm that takes in as input
# the previous region growing results
code <- '
    Rcpp::List list_voxel_neis(R_list_voxel_neis);
    Rcpp::IntegerVector clust_ns(R_clust_ns);
    
    Rcpp::IntegerVector r_clusts(R_clusts);
    arma::Col<int> clusts(r_clusts.begin(), r_clusts.size(), false);
    int nvoxs = clusts.n_elem;
    
    Rcpp::NumericMatrix r_tfunc(R_tfunc);
    arma::mat tfunc(r_tfunc.begin(), r_tfunc.nrow(), r_tfunc.ncol(), false);

    Rcpp::NumericMatrix r_centers(R_centers);
    arma::mat centers(r_centers.begin(), r_centers.nrow(), r_centers.ncol(), false);
    
    int nchanges = 0;
    for (int vi=0; vi<nvoxs; vi++) {
        
        SEXP ll = list_voxel_neis[vi];
        arma::uvec neis = Rcpp::as<arma::uvec>(ll) - 1;
        arma::Col<int> neis_clusts = clusts(neis);
        arma::Col<int> uniq_clusts = arma::unique(neis_clusts);
        
        // Neighbors are all part of the same cluster
        if (uniq_clusts.n_elem == 1) {
            continue;
        }
                
        // Distance to neighboring cluster centers
        // and choose the best (cluster with smallest distance to voxel)
        int n = uniq_clusts.n_elem;
        double nei_dist; double best_dist = 0; int best_clust = 1; int uc;
        for (int ci=0; ci<n; ci++) {
            uc = uniq_clusts(ci);
            nei_dist = arma::sum(arma::square(tfunc.row(vi) - centers.row(uc-1)));
            if (ci == 0 || nei_dist < best_dist) {
                best_dist  = nei_dist;
                best_clust = uc;
            }
        }
        
        // If the best cluster is the same, move on
        int curr_clust = clusts[vi];
        if (best_clust == curr_clust) {
            continue;
        }
        
        // Otherwise, make the change
        nchanges = nchanges + 1;
        double old_n = (double)clust_ns[curr_clust];
        double new_n = (double)clust_ns[best_clust];
        // Remove from old center
        centers.row(curr_clust-1) = (centers.row(curr_clust-1)*old_n - tfunc.row(vi)) / (old_n - 1);
        clust_ns[curr_clust-1] = clust_ns[curr_clust-1] - 1;
        // Add to new center
        centers.row(best_clust-1) = (centers.row(best_clust-1)*new_n + tfunc.row(vi)) / (new_n + 1);
        clust_ns[best_clust-1] = clust_ns[best_clust-1] + 1;
        // Assign to new cluster
        clusts(vi) = best_clust;
    }
        
    printf("changes -> %i\\n", nchanges);
    
    return Rcpp::wrap(clusts);
'

spatial_kmeans_worker <- cxxfunction(
    signature(R_clusts="integer", R_list_voxel_neis="List", R_tfunc="numeric", R_centers="numeric", R_clust_ns="integer"), 
    code, plugin="RcppArmadillo"
)

# calculate neighbors for each node
regions <- comp[comp!=0]
mask3d[mask][comp==0] <- F
tfunc <- tfunc[comp!=0,]
nregions <- length(unique(regions))

lst_neis <- find_neighbors_masked(mask3d, include.self=F, verbose=F)
# calculate region time-series
region_ts <- sapply(1:nregions, function(ri) {
    colMeans(tfunc[regions==ri,,drop=F])
})
region_ts <- scale(region_ts)
tregion_ts <- t(region_ts)
region_ns <- as.integer(sapply(1:nregions, function(i) sum(regions==i)))
copy.clusts1 <- as.integer(regions)

system.time(comp2 <- as.numeric(spatial_kmeans_worker(copy.clusts1, lst_neis, tfunc, tregion_ts, region_ns)))

copy.clusts2 <- as.integer(comp2)
system.time(comp3 <- as.numeric(spatial_kmeans_worker(copy.clusts2, lst_neis, tfunc, tregion_ts, region_ns)))
write.nifti(comp3, hdr, as.vector(mask3d), outfile="step03_region_growing_kmeans_c.nii.gz", overwrite=T)
