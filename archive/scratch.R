basedir="/Users/zarrar/Dropbox/Code/rparcellate"

# Settings
cor.method="pearson"
solutions=2:20

# Read image (4D and cingulate)
library(Rniftilib)

# Have input, roi, a, b, c, replace_a=..., replace_b=..., replace_c=..., outputdir
# cor_file, clust_file, 

ts.file = file.path(basedir, "example", "lfo_residuals2standard_4mm.nii.gz")
roi.file = file.path(basedir, "example", "cingulate_sag46_4mm.nii.gz")
mask.file = file.path(basedir, "example", "lfo_mask_4mm.nii.gz")

cor.file = file.path(basedir, "example", "cingulate_cor.Rda")
cluster.file = file.path(basedir, "example", "cingulate_cluster.Rda")
validation.file = file.path(basedir, "example", "cingulate_valid.Rda")

## Load data
ts.img = nifti.image.read(ts.file)
roi.img = nifti.image.read(roi.file)
mask.img = nifti.image.read(mask.file)

## Confirm that first 3 dimensions are the same
if (sum(dim(roi.img)==dim(ts.img)[1:3])!=3) {
    stop("Dimensions of functional and ROI image to not match\n")
}

## Get coords for 3D ROI
roi.coords = nifti.get.coords(roi.img)

# Convert the 4D data into 2D matrix based on regions in ROI>0
ts.mat = nifti.4Dto2D(ts.img, roi.coords)

# Correlate
ts.cor = cor(ts.mat, method=cor.method)

# Reduce based on kNN
# DECIDE IF YOU WANT THIS
#ts.cor.knn = cluster.kNN(1-ts.cor, knn)

# Compute spectral clustering
tmp.specc = cluster.spectral(ts.cor, solutions)
ts.clusters = tmp.specc$clusters

# Validation measures
dist.mat = ts.cor*-1
sim.mat = ts.cor
cluster.validation = cluster.valid.all(dist.mat, sim.mat, ts.clusters)

# Save
save(ts.cor, file=cor.file)
save(ts.clusters, file=cluster.file)
save(cluster.validation, file=validation.file)

