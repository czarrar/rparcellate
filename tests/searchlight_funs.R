# Test that 27 nearest neighbors case works
test_that("given a 3x3x3 array, node neighbors works", {
	dims <- c(3, 3, 3)
	
	arr <- array(1:prod(dims), dims)
	ref.offsets <- as.vector(arr) - arr[2,2,2]
	
	comp.offsets <- node_neighbors(dims)
	
	expect_that(ref.offsets, equals(comp.offsets))
})

test_that("given a 91x109x91 array, node neighbors works", {
	dims <- c(91, 109, 91) # standard 2mm brain
	
	arr <- array(1:prod(dims), dims)[1:3,1:3,1:3]
	ref.offsets <- as.vector(arr) - arr[2,2,2]
	
	comp.offsets <- node_neighbors(dims)
	
	expect_that(ref.offsets, equals(comp.offsets))
})


test_that("given a 10x10x10 array, find neighbors works", {
	dims <- c(10, 10, 10)
	nei.dist <- 2
	
	mask <- array(rnorm(prod(dims)) > 0, dims)
	vmask <- as.vector(mask)
	nodes <- which(vmask)
	
	crds <- expand.grid(list(x=1:dims[1], y=1:dims[2], z=1:dims[3]))
	d <- as.matrix(dist(crds))
	ref.nei_by_node <- sapply(nodes, function(i) which(d[i,] < nei.dist & vmask))
	
	comp.nei_by_node <- find_neighbors(mask, nodes, nei.dist=nei.dist)
	
	expect_that(ref.nei_by_node, equals(comp.nei_by_node))
})
