
# Playing with converting to C code
		new_regions <- regions[]
        
        # calculate region time-series
        region_ts <- sapply(1:nregions, function(ri) {
            rowMeans(func[,new_regions==ri,drop=F])
        })
	    
		## Loop through each region
		
		# Other thought here is to get neighbors for all the nodes
		# at once. In fact this would be the most efficient approach!
		lst_neis <- find_neighbors_masked(mask3d, region_nodes, 
										include.self=FALSE, verbose=FALSE)
		# want to unwrap?
		
		arma::uvec neis = Rcpp::as<arma::uvec>(ll) - 1;
		
		# for (ri in 1:nregions) {
		
		# could loop through lst_neis and save to new variable
		# or could have lst_neis has matrix and null neis are set to 0
		# i think the first option is easier. would be nice if i could 
		# check if unique each time added. otherwise could just have the
		# output vector be a list as long as the nnodes and assign values to it
		# then i can subset the func data frame
			region_nodes = arma::find(regions == (ri+1));
			neis_vec.zeros()
			for (int rni=0; rni<region_nodes.n_elem; rni++) {
				SEXP ln = lst_neis[region_nodes(rni) - 1];
    			arma::uvec nei_inds = Rcpp::as<arma::uvec>(ln) - 1;
				neis_vec.elem(nei_inds).ones();
			}
		# so the neis here will hold all the neighbor hood indices
			nei_inds = arma::find(neis_vec);
			arma::colvec rn_cors = region_ts.col(ri) * func_ts.cols(nei_inds);
		
		# get the correlation threshold
			max_cor = 0.9 .* rn_cors.max();
		
		# grow only to neighborhood vertices greater than thresh
			neis_to_join = nei_inds[arma::find(rn_cors > max_cor)]
		
		# also check if another region did not grow into the node
		   # loop through this neis_to_join
			for (int ni=0; ni<neis_to_join.n_elem; ni++) {
				int nei_to_join = neis_to_join[ni];
				if (rn_cors[nei_to_join] > nei_cors[nei_to_join]) { # rn_cors and nei_cors have different indices
					nei_cors[nei_to_join] = rn_cors[nei_to_join]
					regions[nei_to_join] = ri + 1;
				}				
			}
		
			nei_cors[neis_to_join[ni]] = rn_cors;
		
		# i can either compute the correlation or the dot-product if it is already scaled
		
		# OH check if region_ts are scaled
		
		# }
		
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


