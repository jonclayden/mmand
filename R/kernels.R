isKernel <- function (object)
{
    return (is.list(object) && "kernel" %in% class(object))
}

sampleKernel <- function (kernel, values)
{
    if (!isKernel(kernel))
        output(OL$Error, "Specified kernel is invalid")
    if (kernel$name == "discrete")
        output(OL$Warning, "The \"sampleKernel\" function does not produce useful output for discrete kernels")
    
    return (.Call("sample_kernel", kernel, as.numeric(values), PACKAGE="mmand"))
}

plot.kernel <- function (x, y, xlim = c(-2,2), ...)
{
    values <- seq(xlim[1], xlim[2], length.out=101)
    plot(values, sampleKernel(x,values), xlab="x", ylab="k(x)", xlim=xlim, ...)
}

discreteKernel <- function (values)
{
    if (isKernel(values))
        return (values)
    else
    {
        values <- as.array(values)
        if (!is.numeric(values))
            report(OL$Error, "Kernel must be numeric")
        storage.mode(values) <- "double"
        return (structure(list(name="discrete", values=values), class="kernel"))
    }
}

shapeKernel <- function (width, dim = length(width), type = c("box","disc","diamond"), brush = TRUE, normalised = FALSE)
{
    type <- match.arg(type)
    
    if (dim > length(width))
        width <- rep(width, length.out=dim)
    else if (dim < length(width))
        width
    
    widthCeiling <- ceiling(width)
    scaleFactors <- max(width) / width
    
    size <- ifelse(widthCeiling %% 2 == 1, widthCeiling, widthCeiling+1)
    kernel <- array(0, dim=size)
    
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
        kernel <- ifelse(kernel < 0.5, 0L, 1L)
    else if (normalised)
        kernel <- kernel / sum(kernel, na.rm=TRUE)
    
    return (discreteKernel(kernel))
}

gaussianKernel <- function (sigma, dim = length(sigma), size = 4*sigma, normalised = TRUE)
{
    if (dim > length(sigma))
        sigma <- rep(sigma, length.out=dim)
    
    if (dim > length(size))
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
    
    return (discreteKernel(kernel))
}

boxKernel <- function ()
{
    return (structure(list(name="box"), class="kernel"))
}

triangleKernel <- function ()
{
    return (structure(list(name="triangle"), class="kernel"))
}

mitchellNetravaliKernel <- mnKernel <- function (B = 1/3, C = 1/3)
{
    return (structure(list(name="mitchell-netravali", B=B, C=C), class="kernel"))
}
