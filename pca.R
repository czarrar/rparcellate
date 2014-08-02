migp <- function(Y, scale=TRUE) {
  # require(rARPACK)
  # note: eigs function in R requires k to not be as large as data dim - 1
  
  s <- length(Y)                # number of subjects
  t <- nrow(Y[[1]])             # number of time-points 
  v <- ncol(Y[[1]])             # number of voxels
  n <- t                        # number of components to output
  
  if (scale) {
    for (i in 1:s) {
      Y[[i]] <- scale(Y[[i]])
    }
  }
  
  o <- sample(s)                # randomise the order subjects are processed
  W <- Y[[o[1]]]                # copy 1st (randomly chosen) subject Y into W
  for (i in 2:s) {              # main loop over all other subjects
    W <- rbind(W, Y[[o[i]]])    # concatenate W with the next subject
    WtW <- tcrossprod(W)        # efficiently get the top 'temporal' eigen-
    r <- eigen(WtW, T)          # vectors of W
    U <- r$vectors[,1:(2*t-1)]  # (note: we could use eigs function instead)
    W <- crossprod(U, W)        # and multiply these into W to get weighted 
                                # spatial eigenvectors
  }
  
  W[1:n,]                       # output just the required number of strongest spatial eigenvectors
}

tcpca <- function(Y, scale=T) {
  require(rARPACK)
  require(plyr)
  
  s <- length(Y)                # number of subjects
  t <- nrow(Y[[1]])             # number of time-points 
  v <- ncol(Y[[1]])             # number of voxels
  n <- t                        # number of components to output
  
  if (scale) {
    for (i in 1:s) {
      Y[[i]] <- scale(Y[[i]])
    }
  }
  
  # collapse all the subjects together
  Z <- laply(Y, function(y) y)
  dim(Z) <- c(dim(Z)[1]*dim(Z)[2], dim(Z)[3])
  
  ZtZ <- tcrossprod(Z)  # covariance of (subjs x time) x (subjs x time)
  W <- eigs(ZtZ, n)     # keep only n components
  W <- crossprod(W$vectors, Z) # now get the spatial eigenvectors
  
  W
}

tconcat <- function(Y, scale=T) {
  require(plyr)
  
  s <- length(Y)                # number of subjects
  t <- nrow(Y[[1]])             # number of time-points 
  v <- ncol(Y[[1]])             # number of voxels
  n <- t                        # number of components to output
  
  if (scale) {
    for (i in 1:s) {
      Y[[i]] <- scale(Y[[i]])
    }
  }
  
  # collapse all the subjects together
  Z <- laply(Y, function(y) y)
  dim(Z) <- c(dim(Z)[1]*dim(Z)[2], dim(Z)[3])
  
  return(Z)
}
  

cscov <- function(Y, scale=T) {
  require(plyr)
  
  s <- length(Y)                # number of subjects
  t <- nrow(Y[[1]])             # number of time-points 
  v <- ncol(Y[[1]])             # number of voxels
  n <- t                        # number of components to output
  
  if (scale) {
    for (i in 1:s) {
      Y[[i]] <- scale(Y[[i]])
    }
  }
  
  Z <- laply(Y, cov)
  Z <- apply(Z, c(2,3), mean)
  
  Z
}


Y <- lapply(1:10, function(i) matrix(rnorm(100*500), 100, 500))
compW <- migp(Y, T)
refW <- tcpca(Y, T) # doing something wrong here

cY <- tconcat(Y, T)
cWY <- cov(cY)

cW <- cscov(Y, T)

sqrt(mean((cWY - cW)^2)) # note that average of the covariance and covariance of temporally concatenated data are one and the same

tmp1 <- cor(compW)[1,]
tmp2 <- cor(refW)[1,]
tmp3 <- cov2cor(cWY)[1,]
tmp4 <- cov2cor(cW)[1,]
which(tmp1>0.2)
which(tmp2>0.2)
which(tmp3>0.06)
as.numeric(which(tmp4>0.06))

# want temporal stuff so flip
#p1 <- princomp(t(cY))
#dim(p1$scores)
p2 <- prcomp(t(cY))
dim(p2$x) # note: only gives up to 500 components because can't be larger than # of voxs
cor(t(p2$x))[1:5,1:5]
cW[1:5,1:5]
cor(as.numeric(cor(t(p2$x[,1:100]))), as.numeric(cW))
cor(as.numeric(cor(compW)), as.numeric(cW)) # exact same result as above
cor(as.numeric(cor(refW)), as.numeric(cW)) # exact same result as above

barplot((p2$sdev/sum(p2$sdev))*100)
sum((p2$sdev/sum(p2$sdev))[1:450]) # so maybe use number of components that estimate 95% of the variance
# or i could do the simpler thing and not do this temporal reduction
