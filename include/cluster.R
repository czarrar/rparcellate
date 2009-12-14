# ---------------------
# CLUSTERING ALGORITHMS
# ---------------------

cluster.kmeans = function(x, centers, iter.max=200, nstart=20, algorithm="Hartigan-Wong") {
    library(cluster)
    return(kmeans(x, centers, iter.max, nstart, algorithm))
}

# Function to calculate spectral clustering solution to distance matrix
## Arguments
## - S: similarity matrix
## - ks: vector of different solutions to calculate
## Results
## - list with members: clusters, eigval, eigvec
##   - clusters: matrix of voxels/regions (rows) x solutions (columns)
##   - eigval: vector of eigen values
##   - eigvec: matrix of eigen vectors
cluster.spectral = function(S, ks) {
    x = as.matrix(S)
    m = nrow(x)
    
    # 1. Local-Scaling
    s = rep(0,nrow(x))
    dota = rowSums(x*x)/2
    dis = crossprod(t(x))
    for (i in 1:m)
        dis[i,] = 2*(-dis[i,] + dota + rep(dota[i],m))  
    dis[dis < 0] = 0    ## fix
    for (i in 1:m)
        s[i] = median(sort(sqrt(dis[i,]))[1:5])
    km = exp(-dis / s%*%t(s))
    
    # 2. Normalize based on row sums
    d = 1/sqrt(rowSums(km))
    l = d * km %*% diag(d)
    
    # 3. Eigen decomposition
    eig = eigen(l, T)
    
    # 4. Loop through and compute cluster solution for each k
    specc.solutions = matrix(NA, m, length(ks))
    for (kk in 1:length(ks)) {
        k = ks[kk]
        xi = eig$vectors[,1:k]
        yi = xi/sqrt(rowSums(xi^2)) # Normalize based on row sums
        specc.solutions[,kk] = cluster.kmeans(yi, k)$cluster
    }
    
    names(dim(specc.solutions)) = c("voxels", "solutions")
    colnames(specc.solutions) = ks
    rownames(specc.solutions) = rownames(S)
    
    return(list(clusters=specc.solutions, eigval=eig$values, eigvec=eig$vectors))
}

# ---------------------
# CLUSTERING VALIDATION
# ---------------------

cluster.valid.r2 = function(dist.mat, clusters.mat) {
    library(vegan)
    
    num.regions = nrow(clusters.mat)
    num.solutions = ncol(clusters.mat)
    
    vec.r2 = c()
    
    for (kk in 1:num.solutions) {
        tmp.adonis = adonis(as.dist(dist.mat) ~ as.factor(clusters.mat[,kk]), permutations=1)
        vec.r2[kk] = tmp.adonis$aov.tab$R2[1]
    }
    
    return(vec.r2)
}

# ---------------------------------------
# CLUSTERING DISTANCE/CONSISTENCY METRICS
# ---------------------------------------

# Function to use with VI distance
logg = function(mmm) {
    logm = log(mmm)
    logm[logm == -Inf] = 0
    return(logm)
}

cluster.dist.vi = function(mem1, mem2, norm=F) {
    nn = length(mem1)
    if (length(mem2) != nn) stop('Need same number of vertices in each
    partition.')
    c1 = length(unique(mem1))
    c2 = length(unique(mem2))

    p1 = matrix(NA, nn, c1)
    for (i in 1:c1) p1[ , i] = (mem1 == i)
    p2 = matrix(NA, nn, c2)
    for (i in 1:c2) p2[ , i] = (mem2 == i)

    p12 = crossprod(p1, p2) / nn
    meila = sum(p12 * (log(apply(p12, 1, sum) %*% t(apply(p12, 2, sum))) - 2
    * logg(p12)))
    varinf = if (norm) - meila / sum(p12 * logg(p12)) else meila
    varinf
}

cluster.dist.cramersv = function(vec1, vec2) {
    library(vcd)
    tab = table(vec1, vec2)
    return(assocstats(tab)$cramer)
}


# ----------------------
# CLUSTERING OTHER STUFF
# ----------------------


###
# Match 2 different cluster memberships
###
# INPUT:
## template.m - ideal vector of cluster memberships
## comparison.m - vector of cluster memberships to match to ideal
## match.method - method in matchClasses (rowmax, greedy, exact)
cluster.match.members = function(template.m, comparison.m, match.method="exact", iter=10) {
	library(e1071)
	num.regions = length(comparison.m)
	new.m = vector(mode="numeric", length=num.regions)

	if(length(unique(template.m)) != length(unique(comparison.m))) {
	
		return(matchUnequalClasses(template.m, comparison.m))
	} else {
		tab = table(template.m, comparison.m)
		conversion = matchClasses(tab, method=match.method, iter=iter)

		for (template.val in 1:num.regions) {
			comparison.val = conversion[template.val]
			new.m[comparison.m==comparison.val] = template.val
		}
	
		return(new.m)
	}
}

matchUnequalClasses = function(template.m, comparison.m) {
	tab = table(comparison.m, template.m)
	
	# Checks
	if(nrow(tab)<ncol(tab))
		stop("template.m must have less partitions than comparison.m")
	if(min(unique(template.m))!=1 | min(unique(comparison.m))!=1)
		stop("partitions template.m and comparison.m must start with 1")
		
	num.col = ncol(tab)
	
	for(i in 1:num.col) {
		which.leading = order(tab[,i], decreasing=TRUE)[1]
		if(which.leading!=i & which.leading>i) {
			tmp.m = comparison.m
			tmp.m[comparison.m==i] = which.leading
			tmp.m[comparison.m==which.leading] = i
			comparison.m = tmp.m
			tab = table(comparison.m, template.m)
		}
	}
	
	return(comparison.m)
}
