shapeKernel <- function (width, dim = length(width), type = c("box","disc","diamond"), brush = TRUE, binary = TRUE, normalised = FALSE)
{
    type <- match.arg(type)
    bg <- ifelse(brush, NA, 0)
    
    if (dim > length(width))
        width <- rep(width, length.out=dim)
    
    maxWidth <- ceiling(max(width))
    scaleFactors <- max(width) / width
    
    # This corresponds to the smallest isotropic kernel of odd width to fit the shape
    size <- ifelse(maxWidth %% 2 == 1, maxWidth, maxWidth+1)
    kernel <- array(bg, dim=rep(size,dim))
    
    x <- 1:size - (size+1)/2
    nearEdges <- lapply(scaleFactors, "*", x - 0.5*sign(x))
    farEdges <- lapply(scaleFactors, "*", x + 0.5*sign(x))
    
    normFun <- switch(type, box=function(a,b) pmax(abs(a),abs(b)),
                            disc=function(a,b) sqrt(a^2 + b^2),
                            diamond=function(a,b) (abs(a) + abs(b)))
    
    minNorms <- Reduce(function(a,b) outer(a,b,FUN=normFun), nearEdges)
    maxNorms <- Reduce(function(a,b) outer(a,b,FUN=normFun), farEdges)
    
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
        size <- ceiling(4 * max(sigma))
    else if (length(size) > 1)
        report(OL$Error, "The kernel size must be given as a single integer")
    else
        size <- ceiling(size)
    
    size <- ifelse(size %% 2 == 1, size, size+1)
    
    scaleFactors <- max(sigma) / sigma
    x <- 1:size - (size+1)/2
    centres <- lapply(scaleFactors, "*", x)
    
    normFun <- function(a,b) sqrt(a^2 + b^2)
    norms <- Reduce(function(a,b) outer(a,b,FUN=normFun), centres)
    
    kernel <- array(dnorm(norms,sd=max(sigma)), dim=rep(size,dim))
    
    if (normalised)
        kernel <- kernel / sum(kernel, na.rm=TRUE)
    
    return (kernel)
}
