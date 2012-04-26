shapeKernel <- function (width, dim = length(width), type = c("box","disc","diamond"), brush = TRUE, binary = TRUE, normalised = FALSE)
{
    type <- match.arg(type)
    bg <- ifelse(brush, NA, 0)
    
    if (dim > length(width))
        width <- rep(width, length.out=dim)
    else if (dim < length(width))
        width
    
    widthCeiling <- ceiling(width)
    scaleFactors <- max(width) / width
    
    size <- ifelse(widthCeiling %% 2 == 1, widthCeiling, widthCeiling+1)
    kernel <- array(bg, dim=size)
    
    x <- lapply(1:dim, function(i) 1:size[i] - (size[i]+1)/2)
    nearEdges <- lapply(1:dim, function(i) scaleFactors[i] * (x[[i]] - 0.5*sign(x[[i]])))
    farEdges <- lapply(1:dim, function(i) scaleFactors[i] * (x[[i]] + 0.5*sign(x[[i]])))
    
    normFun <- switch(type, box=function(a,b) pmax(abs(a),abs(b)),
                            disc=function(a,b) sqrt(a^2 + b^2),
                            diamond=function(a,b) (abs(a) + abs(b)))
    
    if (dim == 1)
    {
        minNorms <- abs(nearEdges[[1]])
        maxNorms <- abs(farEdges[[1]])
    }
    else
    {
        minNorms <- Reduce(function(a,b) outer(a,b,FUN=normFun), nearEdges)
        maxNorms <- Reduce(function(a,b) outer(a,b,FUN=normFun), farEdges)
    }
    
    maxDistance <- max(width) / 2
    indices <- minNorms < maxDistance
    kernel[indices] <- pmin(1, (maxDistance - minNorms[indices]) / (maxNorms[indices] - minNorms[indices]))
    
    if (binary)
        kernel <- ifelse(kernel < 0.5, as.integer(bg), 1L)
    else if (normalised)
        kernel <- kernel / sum(kernel, na.rm=TRUE)
    
    return (kernel)
}

gaussianKernel <- function (sigma, dim = length(sigma), size = NULL, normalised = TRUE)
{
    if (dim > length(sigma))
        sigma <- rep(sigma, length.out=dim)
    
    if (is.null(size))
        size <- ceiling(4 * sigma)
    else if (dim > length(size))
        size <- rep(ceiling(size), length.out=dim)
    else
        size <- ceiling(size)
    
    size <- ifelse(size %% 2 == 1, size, size+1)
    
    scaleFactors <- max(sigma) / sigma
    x <- lapply(1:dim, function(i) 1:size[i] - (size[i]+1)/2)
    centres <- lapply(1:dim, function(i) scaleFactors[i] * x[[i]])
    
    normFun <- function(a,b) sqrt(a^2 + b^2)
    norms <- Reduce(function(a,b) outer(a,b,FUN=normFun), centres)
    
    kernel <- array(dnorm(norms,sd=max(sigma)), dim=size)
    
    if (normalised)
        kernel <- kernel / sum(kernel, na.rm=TRUE)
    
    return (kernel)
}
