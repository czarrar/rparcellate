#!/usr/bin/env Rscript

library(XML)
library(plyr)
suppressMessages(library(niftir))
library(biganalytics)
library(stringr)

fsldir	<- Sys.getenv("FSLDIR")

hemisphere_rois <- function(hdr) {
  # Get the center coordinate
  center.coord <- as.numeric(hdr$qto.ijk %*% c(0,0,0,1))[-4]
  cat("xyz: 0 0 0 => ijk:", center.coord, "\n")
  
  # Note the above is based on the first element being 0
  # need to convert 0 => 1
  cx <- center.coord[1] + 1 # center of x-axis
  ex <- hdr$dim[1] # end of x-axis
  
  # Create images
  lh_roi <- array(1, hdr$dim)
  lh_roi[1:cx,,] <- 0	
  rh_roi <- array(1, hdr$dim)
  rh_roi[cx:ex,,] <- 0
  
  return(list(lh=as.vector(lh_roi)!=0, rh=as.vector(rh_roi)!=0))
}

# 1. Information on cortical/subcortical and cerebellum units
## read in harvard-oxform parcellation info
ho.xml.file <- file.path(fsldir, "data/atlases/HarvardOxford-Subcortical.xml")
ho.xml <- xmlInternalTreeParse(ho.xml.file)
ho.xml.data <- xmlRoot(ho.xml)[2]$data
ho.lst <- xmlToList(ho.xml.data)
ho.df0 <- ldply(ho.lst, function(lst) {
  return(c(
    index = as.numeric(lst$.attrs[1]), 
    name = as.character(lst$text)
  ))
})
## split left/right
ho.df <- data.frame(
  index = as.numeric(ho.df0$index) + 1, 
  hemi = str_extract(ho.df0$name, "Left|Right"), 
  name = str_replace(ho.df0$name, "Left |Right ", ""), 
  full.name = ho.df0$name
)
# ## read in cerebellum information
# ce.xml.file <- file.path(fsldir, "data/atlases/Cerebellum_MNIfnirt.xml")
# ce.xml <- xmlInternalTreeParse(ce.xml.file)
# ce.xml.data <- xmlRoot(ce.xml)[2]$data
# ce.lst <- xmlToList(ce.xml.data)
# ce.df0 <- ldply(ce.lst, function(lst) {
#   return(c(
#     index = as.numeric(lst$.attrs[1]), 
#     name = as.character(lst$text)
#   ))
# })
# ## split left/right
# ce.df <- data.frame(
#   index = as.numeric(ce.df0$index) + 1, 
#   hemi = str_extract(ce.df0$name, "Left|Right"), 
#   name = str_replace(ce.df0$name, "Left |Right ", ""), 
#   full.name = ce.df0$name
# )



# 2. Read in probability maps
## read in header
hdr       <- read.nifti.header(file.path(fsldir, "data/standard/MNI152_T1_2mm.nii.gz"))
## read in probability maps
ho.file   <- file.path(fsldir, "data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz")
ce.file   <- file.path(fsldir, "data/atlases/Cerebellum/Cerebellum-MNIfnirt-prob-2mm.nii.gz")
ho.probs  <- read.big.nifti(ho.file)
ce.probs  <- read.big.nifti(ce.file)
## since the cerebellum will be considered as one unit collapse
ce.mprobs <- colmax(ce.probs)
## add the cerebellum probs onto the others
probs     <- rbind(ho.probs[,], ce.mprobs)
## mask data
mask      <- apply(probs, 2, max) != 0
probs     <- probs[,mask] # total of 286,933 voxels

# 3. Max probabilities / assign labels
## assign each voxel to a label based on some max probability threshold
thrs      <- c(0, 10, 20, 25, 50)
#thrs <- c(20,25)
all.max.probs <- t(laply(thrs, function(thr) {
  cat("thr:", thr, "\n")
  aaply(probs, 2, function(x) {
    ## exclude white-matter to allow more grey-matter
    x[c(1,12)] <- 0
    ## set white-matter threshold to two times thr
    ## so allow more grey matter
    #x[c(1,12)] <- x[c(1,12)] * (x[c(1,12)] > (thr*2))
    ind <- which.max(x)
    ifelse(x[ind]>thr, ind, 0)
  }, .progress="text")
}))

# 4. Determine maps for three levels
## get left/right hemisphere ROIs
hemis <- hemisphere_rois(hdr)
hemis$lh <- hemis$lh[mask]
hemis$rh <- hemis$rh[mask]
## split cortex, subcortical, cerebellum, and subcortical and left/right
mappings <- list(
  cortex = c(2,13), 
  subcortical = c(4:11,15:21), 
  cerebellum = 22
)
## recombine max.probs based on new mapping above
all.max.probs2 <- t(aaply(all.max.probs, 2, function(x) {
  relabel.x <- vector("integer", length(x))
  for (i in 1:length(mappings)) {
    relabel.x[x %in% mappings[[i]]] <- i
  }
  relabel.x
}))
## make the right hemisphere have unique values +10
## and remove the midline
all.max.probs3 <- t(aaply(all.max.probs2, 2, function(x) {
  relabel.x <- vector("integer", length(x))
  relabel.x[hemis$lh] <- x[hemis$lh]
  relabel.x[hemis$rh] <- (x[hemis$rh] + 10) * (x[hemis$rh] != 0)
  relabel.x
}))

# 5. Save
# all done! save this out
dir.create("anatomical_priors", FALSE)
for (i in 1:ncol(all.max.probs3)) {
  outfile <- sprintf("anatomical_priors/fsl_maxprob%02i.nii.gz", thrs[i])
  cat("saving", outfile, "\n")
  write.nifti(all.max.probs3[,i], hdr, mask, outfile=outfile, overwrite=T)
}
