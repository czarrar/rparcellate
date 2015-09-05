# TODO: add some amount of randomization in ordering?
# CITE: Group-PCA for very large fMRI datasets

#' Implements MELODIC's Incremental Group-PCA (MIGP)
#'
#' Returns the top spatial eigenvectors for the group-level PCA
#'
#' @param Y list of subject time-series matrices each tpts x regions
#' @param scale whether or not to scale each subject's / regions time-series 
#'  first (default = TRUE)
#'
#' @export
#' 
#' @return matrix of strongest spatial eigenvectors
#'
#' @examples
#' Y <- lapply(1:10, function(i) matrix(rnorm(100*500), 100, 500))
#' compW <- migp(Y, T)
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

migp_files <- function(func_files, mask, scale=TRUE) {
  # require(rARPACK)
  # note: eigs function in R requires k to not be as large as data dim - 1
  
  read_data <- function(func_file) {
    func  <- read.big.nifti4d(func_file)
    func  <- do.mask(func, mask)
    if (scale) {
      func <- scale_fast(func, to.copy=F)
    }
    as.matrix(func)
  }
  
  s   <- length(func_files)       # number of subjects
  hdr <- read.nifti.header(func_files[1])
  t   <- hdr$dim[4]               # number of time-points 
  v   <- sum(mask)                # number of voxels
  n   <- t                        # number of components to output
  
  cat("1.")
  o   <- sample(s)                # randomise the order subjects are processed
  ff  <- func_files[o[1]]         # copy 1st (randomly chosen) subject Y into W
  W   <- read_data(ff)
  for (i in 2:s) {                # main loop over all other subjects
    cat(i, ".", sep="")
    ff<- func_files[o[i]]         # concatenate W with the next subject
    W <- rbind(W, read_data(ff))      
    WtW <- tcrossprod(W)          # efficiently get the top 'temporal' eigen-
    r <- eigen(WtW, T)            # vectors of W
    U <- r$vectors[,1:(2*t-1)]    # (note: we could use eigs function instead)
    W <- crossprod(U, W)          # and multiply these into W to get weighted 
                                  # spatial eigenvectors
  }
  cat("\n")
  
  W[1:n,]                         # output just the required number of strongest spatial eigenvectors
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

