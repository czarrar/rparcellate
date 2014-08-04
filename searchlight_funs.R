# Functions that do a searchlight on a brain image

#' Neighbors for a node
#'
#' Returns the relative indices of neighbors for a node
#'
#' Given a node at the origin, gives the (vector) coordinates of the 
#' neigbouring nodes. Distance of the neighbors is specified by the user.
#'
#' @params dims dimension of array
#' @params nei maximum value for dimensions in x,y,z axis
#' @params nei.dist maximum euclidean distance of neighbor
node_neighbors <- function(dims, nei=1, nei.dist=2, include.self=TRUE) {
	# What are the neighbors for a given node?
    # moffsets gives these neighbors in ijk form
    moffsets    <- expand.grid(list(i=-nei:nei, j=-nei:nei, k=-nei:nei))
    dist        <- sqrt(rowSums(moffsets^2))
    moffsets    <- moffsets[dist<=nei.dist,]
    # offsets gives these neighbors in vector index form
    offsets     <- moffsets$k*dims[1]*dims[2] + moffsets$j*dims[1] + moffsets$i
    if (!include.self) 
	    offsets     <- offsets[offsets!=0]
	return(offsets)
}


#' Find neighbors for a set of nodes
#' 
#' Returns a list of neighbors for each node within a mask
#'
#' Primary purpose of this function is to function as a searchlight. For each 
#' node, a list of neighbors are returned. These neighbors can be used to apply
#' some function to the searchlight using another function.
#'
#' @param mask logical array indicating nodes to consider as neighbors
#' @param nodes indices of nodes to find neighbors for
#' @param include.self whether to include the given node in the neighbor list
#' the self node is included as the first one in the neighboring node list
#' @param nei radius of sphere/box to look for neighbors in
#' @param nei.dist maximum euclidean distance from center node
#' @param thr.nei minimum proportion (0-1) of neighbors
#' @param verbose extent of print (default: TRUE)
#' 
#' @export
#' 
#' @return list
#' 
#' @examples
#' # TODO
find_neighbors <- function(mask, nodes=which(mask), include.self=TRUE, 
                           nei=1, nei.dist=2, thr.nei=0, verbose=TRUE) 
{
    if (is.vector(mask)) stop("please provide the mask as an array")
	progress	<- ifelse(verbose, "text", "none")
    
    dims        <- dim(mask)
    mask        <- as.vector(mask)
    
	offsets     <- node_neighbors(dims, nei, nei.dist, 
											  include.self)
    # Minimum number of neigbors with node required
	min.nei	<- floor(length(offsets)*thr.nei)
	 
    # Let's start
    # for a given node (peak)
    # we get the neighboring indices
    # and ensure that the indices aren't outside the box
    # and ensure that the indices are within the brain mask
    max_ind     <- length(mask) + 1
    nei_by_node <- llply(1:length(nodes), function(i) {
        node_ind <- nodes[i]
        nei_inds <- offsets + node_ind
        nei_inds <- nei_inds[nei_inds > 0 & nei_inds < max_ind]
        nei_inds <- nei_inds[mask[nei_inds]]
		  if (length(nei_inds) < min.nei) nei_inds <- c()
		  return(nei_inds)
    }, .progress=progress, .parallel=FALSE)
    names(nei_by_node) <- nodes
    #attr(nei_by_node, "offsets") <- offsets
	 
	 if (verbose) {
		  nzero <- sum(sapply(nei_by_node, length) == 0)
		  if (nzero > 0) cat(nzero, "nil node neighbors\n")
	 }
	 
    return(nei_by_node)
}


#' Converts nodes in a list to indexing elements withing mask
#' 
#' Returns a new list of neighbors for each node where the neighbors are indices in a masked array
#' 
#' @params neis_by_node0 list of neighbors for each node
#' @params mask mask of array
#'
#' @export
#' 
#' @return list
neighbors_array2mask <- function(neis_by_node0, mask) {
    #mat2arr_inds    		<- which(mask)
    arr2mat_inds    		<- mask*1
    arr2mat_inds[mask] <- 1:sum(mask)
    neis_by_node    		<- lapply(neis_by_node0, function(x) arr2mat_inds[x])
	return(neis_by_node)
}


nodes_mask2array <- function(nodes, mask) {
    mat2arr_inds    <- which(mask)
    nodes2          <- mat2arr_inds[nodes]
	return(nodes2)
}

#' Find neighbors for a set of nodes
#' 
#' Returns a list of neighbors for each node within a mask
#' 
#' Neighbors are indices in a masked array. This function combines 
#' `find_neighbors` with `neighbors_array2mask`.
#' 
#' @params mask mask of array
#' @params ... other options to pass onto `find_neighbors`
#' 
#' @export
#' 
#' @return list
find_neighbors_masked <- function(mask, nodes=1:sum(mask), ...) {
    nodes         <- nodes_mask2array(nodes, as.vector(mask))
	neis_by_node0 <- find_neighbors(mask, nodes, ...)
	neis_by_node  <- neighbors_array2mask(neis_by_node0, as.vector(mask))
	return(neis_by_node)
}


#' Perform searchlight analysis
#'
#' Returns vector of searchlight results
#'
#' Applies some function to each node and it's surrounding neighbors
#'
#' @params fun function to apply to each searchlight
#' @params data 2D functional data (ntpts x nregions)
#' @params mask brain mask or something else to constrain nodes/neighbors examined
#' @params process_nodes nodes within the whole image to examine (assumed that this is a subset of mask)
#' @params include_self for each searchlight, do you include the principal or center node
#' @params percent_neighbors the proportion of neighboring nodes that must exist (default=0.5)
#' @params progress what type of progress information to display (default=text)
#' @params parallel apply operations in parallel (default=False)
#'
#' @export
#'
#' @return vector
searchlight <- function(fun, data, mask, process_nodes=which(mask), 
                        include_self=TRUE, percent_neighbors=0.25, 
                        progress="text", parallel=FALSE, ...) 
{
    # Find the neighbors for each node within your mask
    neis_by_node0   <- find_neighbors(mask, nodes=process_nodes, 
                                      include.self=include_self, 
                                      thr.nei=percent_neighbors)
    neis_by_node    <- neighbors_array2mask(neis_by_node0, as.vector(mask))
    
    # Apply function
    sl_res <- laply(neis_by_node, function(neis) {
        if (length(neis) == 0) {
            return(0)
        } else {
            return(fun(data[,neis], ...))
        }
    }, .progress=progress, .parallel=parallel)
    
    return(sl_res)
}
