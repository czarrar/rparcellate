# Create matrix with list of regions and corresponding location in matrix/array
## If you convert a 3D array into a 1D vector, this function allows you to go back to 3D space
## Arguments:
## - nii.img: 3D nifti array
## - to.mask: TRUE or FALSE, mask the coordinates using only non-zero image values
## - use.img.val: TRUE or FALSE, have a column indicating the value of that coordinate
## Return:
## - matrix of coordinates
##   - each row is a voxel
##   - columns are x, y, and z coordinate and if use.img.val=T then value is a 4th column
nifti.get.coords = function(nii.img, to.mask=TRUE, use.img.val=FALSE) {
    # Get elements in each dimension
    tmp.elem = sapply(dim(nii.img), seq)
    ## fix in cases where symmetrical array (e.g. 3x3x3)
    if (is.matrix(tmp.elem)) {
        tmp.new.elem = list()
        for (i in 1:ncol(tmp.elem)) {
            tmp.new.elem[[i]] = tmp.elem[,i]
        }
        tmp.elem = tmp.new.elem
    }
    
    # Get coordinates
    tmp.coord = expand.grid(tmp.elem)
    
    # Get values
    img.vals = as.vector(nii.img[,,])
    
    # Add values in nii.img to coordinates matrix
    if (use.img.val) {
        tmp.coord = cbind(tmp.coord, img.vals)
        colnames(tmp.coord) = c("x", "y", "z", "value")
    } else {
        colnames(tmp.coord) = c("x", "y", "z")
    }
    
    # Get only non-zero coordinates in ROI
    if (to.mask) {
       tmp.coord = tmp.coord[img.vals!=0,] 
    }
        
    return(tmp.coord)
}

# Convert 4D fMRI dataset to a 2D matrix
## Arguments
## - nii.img: Your 4D nifti image
## - coord.mat: Matrix of coordinates to transplant into 2D (use nifti.get.coordinates)
## Results
## - matrix of time-points (rows) x voxels (columns)
nifti.4Dto2D = function(nii.img, coord.mat) {
	# Number of voxels
	num.vox = nrow(coord.mat)

	# Number of Time Points
	num.tp = dim(nii.img)[4]

	# Create new matrix
	new.ts = matrix(NA, num.tp, num.vox)
	
	# Transfer over data
	for (i in 1:num.vox) {
		## Get coordinate
		x = coord.mat[i,1]
		y = coord.mat[i,2]
		z = coord.mat[i,3]
                
		## Copy ts at voxel
		new.ts[,i] = nii.img[x,y,z,]
	}
	
	names(dim(new.ts)) = c("time-points", "voxels")
	rownames(new.ts) = 1:num.tp
	colnames(new.ts) = rownames(coord.mat)

	return(new.ts)
}

# Convert a vector with values for each voxel into a 3D array to be saved as a nifti file
## Arguments
## - vec.data: vector with data
## - coord.mat: Matrix of coordinates to transplant into 3D
## - ref.img: reference image (uses it to get array dimensions)
## - use.rownames: TRUE or FALSE (default=TRUE), rownames have numbers to indicate vector number of voxel in 3D => use them
## Return:
## - 3D array
nifti.1Dto3D = function(vec.data, coord.mat, ref.img, use.rownames=TRUE) {
	# Number of voxels with values
	num.vox = nrow(coord.mat)
	
	# Dimensions
	ref.dim = dim(ref.img)
	
	# Create the new matrix
	new.mat = array(0, ref.dim)
        
	# Use rownames to identify voxel coordinate?
	if (use.rownames)
		coord.vec = as.numeric(rownames(coord.mat))
        
	# Loop through and save data
	for (i in 1:num.vox) {
		if (use.rownames) {
			## Quick copy
			new.mat[coord.vec[i]] = vec.data[i]
		} else {
			## Get coordinates
			x = coord.mat[i,1]
			y = coord.mat[i,2]
			z = coord.mat[i,3]
			## Save
			new.mat[x,y,z] = vec.data[i]
		}
	}
	
	return(new.mat)
}
