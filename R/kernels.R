isKernel <- function (object)
{
    return ("kernel" %in% class(object))
}

isKernelArray <- function (object)
{
    return ("kernelArray" %in% class(object))
}

isKernelFunction <- function (object)
{
    return ("kernelFunction" %in% class(object))
}

plot.kernelArray <- function (x, y, axis = 1, lwd = 2, col = "red", ...)
{
    indices <- as.list(ceiling(dim(x) / 2))
    indices[[axis]] <- 1:(dim(x)[axis])
    line <- do.call("[", c(list(x),indices))
    
    limit <- dim(x)[axis] / 2
    xlim <- limit * c(-1,1)
    xloc <- rep(seq(xlim[1],xlim[2],1), each=2)
    yloc <- c(0, rep(line, each=2), 0)
    
    plot(xloc, yloc, xlab="x", ylab="k(x)", type="l", lwd=lwd, col=col, xlim=xlim, ...)
    abline(v=0, lty=2, col="grey60")
    abline(h=0, lty=2, col="grey60")
}

sampleKernelFunction <- function (kernel, values)
{
    if (!isKernelFunction(kernel))
        report(OL$Error, "Specified kernel is not a valid kernel function")
    
    return (.Call("sample_kernel", kernel, as.numeric(values), PACKAGE="mmand"))
}

plot.kernelFunction <- function (x, y, xlim = c(-2,2), lwd = 2, col = "red", ...)
{
    values <- seq(xlim[1], xlim[2], length.out=101)
    plot(values, sampleKernelFunction(x,values), xlab="x", ylab="k(x)", type="l", lwd=lwd, col=col, xlim=xlim, ...)
    abline(v=0, lty=2, col="grey60")
    abline(h=0, lty=2, col="grey60")
}

kernelArray <- function (values)
{
    if (isKernelArray(values))
        return (values)
    else if (isKernelFunction(values))
        report(OL$Error, "Kernel function cannot be converted to a kernel array")
    else
    {
        values <- as.array(values)
        if (!is.numeric(values))
            report(OL$Error, "Kernel must be numeric")
        storage.mode(values) <- "double"
        return (structure(values, class=c("kernelArray","kernel")))
    }
}

shapeKernel <- function (width, dim = length(width), type = c("box","disc","diamond"), binary = TRUE, normalised = FALSE)
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
    
    return (kernelArray(kernel))
}

gaussianKernel <- function (sigma, dim = length(sigma), size = 6*sigma, normalised = TRUE)
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
    
    return (kernelArray(kernel))
}

kernelFunction <- function (name = c("box","triangle","mitchell-netravali"), ...)
{
    if (is.character(name))
        name <- match.arg(name)
    else if (isKernelFunction(name))
        return (name)
    else if (isKernelArray(name))
        report(OL$Error, "Kernel array cannot be converted to a kernel function")
    else
        report(OL$Error, "Kernel function specification is not valid")
    
    return (structure(list(name=name, ...), class=c("kernelFunction","kernel")))
}

boxKernel <- function ()
{
    return (kernelFunction("box"))
}

triangleKernel <- function ()
{
    return (kernelFunction("triangle"))
}

mitchellNetravaliKernel <- mnKernel <- function (B = 1/3, C = 1/3)
{
    return (kernelFunction("mitchell-netravali", B=B, C=C))
}
