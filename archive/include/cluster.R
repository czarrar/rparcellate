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

# Compute all the cluster validations
cluster.valid.all = function(dist.mat, sim.mat, cluster.mat) {
    cluster.validation = list(q=c(), sil=c(), r2=c(), ks=colnames(cluster.mat))
    
    for(kk in 1:ncol(cluster.mat)) {
        ## Compute modularity
        cluster.validation$q[kk] = cluster.valid.modularity(sim.mat, cluster.mat[,kk])

        ## Compute modified silhouette width
        cluster.validation$sil[[kk]] = cluster.valid.modsilhouette(dist.mat, cluster.mat[,kk])

        ## Compute adonis for each cluster solution?
        cluster.validation$r2[kk] = cluster.valid.r2(dist.mat, cluster.mat[,kk])
    }
    
    return(cluster.validation)
}

cluster.valid.r2 = function(dist.mat, vec.k) {
    library(vegan)
    dist.mat = as.dist(dist.mat)
    vec.k = as.factor(vec.k)
    tmp.adonis = adonis(dist.mat ~ vec.k, permutations=1)
    return(tmp.adonis$aov.tab$R2[1])
}

# Computes a modified silhouette width for each cluster solution
# (based on code from Clare Kelly at med nyu edu)
##
## S(i) = (min(AVGD_BETWEEN(i,k)) - AVGD_WITHIN(i)) / max(AVGD_WITHIN(i), min(AVGD_BETWEEN(i,k)))
## where AVGD_WITHIN(i) is the average distance from the i-th point to the
## other points in its own cluster, and AVGD_BETWEEN(i,k) is the average
## distance from the i-th point to points in another cluster k.
##
## Arguments
## - dist.mat: distance matrix
## - vec.k: vector of k partitions of dist.mat (k>1 and is integer)
## Return
## - vector containing silhoette width values (one for each cluster)
cluster.valid.modsilhouette = function(dist.mat, vec.k) {
    dist.mat = as.matrix(dist.mat)
    num.regions = rownames(dist.mat)
    ks = unique(vec.k)
    num.k = length(ks)
    
    # Keep only lower half of distance matrix (for taking accurate means later)
    diag(dist.mat) = NA
    dist.mat[upper.tri(dist.mat)] = NA
    dist.mat[dist.mat==0] = NA
    
    # Vector to hold silhoette values for each cluster partition
    msilw = c()
    
    # Loop through each partition
    for (kk in 1:num.k) {
        k = ks[kk]
        
        ## Compute average within cluster distance
        if (sum(vec.k==k)==1) {
            within.dist = 0
        } else {
            within.dist = mean(dist.mat[vec.k==k,vec.k==k], na.rm=TRUE)
        }
        
        ## Compute all possible between cluster distance
        ## only use the smallest between-cluster dissimilarity
        between.dist = c()
        for (bk in ks[ks!=k]) {
            tmp.btw = mean(c(dist.mat[vec.k==k,vec.k==bk], dist.mat[vec.k==bk,vec.k==k]), na.rm=TRUE)
            between.dist = min(between.dist, tmp.btw)
        }
        
        ## Compute silhouette width for this partition
        msilw[kk] = (between.dist-within.dist)/max(between.dist,within.dist)
    }
    
    return(msilw)
}

# Computes the modularity (Q)
## code is based on brain connectivity toolbox (http://www.brain-connectivity-toolbox.net/)
##
## Arguments
## - A: weighted undirected network (adjacency or weights)...i.e. similarity matrix
## - p: partition labels
## Return
## - Q: modularity score
##
## References: Newman (2006) -- Phys Rev E 74:036104; PNAS 23:8577-8582.
cluster.valid.modularity = function(A, p) {
    K = apply(A, 2, sum)
    N = length(p)
    m = sum(K)
    B = A - (K %*% t(K))/m
    s = matrix(rep(p, times=N), N, N)
    Q = ((!(s-t(s)))*1) * B/m
    Q = sum(Q)
    return(Q)
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

#cluster.kNN <- function(dist.mat, kNN, mutual=F) {
#    dist.mat = as.matrix(dist.mat)
#    num.regions = ncol(dist.mat)
#    
#    # Loop through columns
#    for (i in 1:num.regions) {
#        ## Set specific element that not in top kNN to 0
#        tmp.which = order(dist.mat[,i])[-1]
#        tmp.which = tmp.which[(kNN+1):length(tmp.which)]
#        dist.mat[tmp.which,i] = 0
#    }
#    
#    if (mutual) {
#        mask <- t(dist.mat)!=0
#        dist.mat <- dist.mat * mask
#    }
#    
#    return(dist.mat)
#}

# Creates a sparse matrix by only keep the k smallest connections for each node
cluster.kNN <- function(mat, kNN, isSimMat=F, mutual=F) {
    mat = as.matrix(mat)
        
    mask <- c(T, rep(F, kNN), rep(T, nrow(mat)-kNN-1))
    
    mat <- apply(mat, 2, function(x) {
        ind <- order(x, decreasing=isSimMat)[mask]
        x[ind] <- 0
        return(x)
    })
    
    if (mutual)
        mat <- mat * ((t(mat)!=0)*1)
    
    return(mat)
}

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
