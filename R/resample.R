resample <- function (x, points, kernel, ...)
{
    UseMethod("resample")
}

resample.default <- function (x, points, kernel, pointType = c("auto","general","grid"), ...)
{
    x <- as.array(x)
    if (!is.numeric(x))
        report(OL$Error, "Target array must be numeric")
    
    nDims <- length(dim(x))
    
    if (nDims == 1 && !is.matrix(points) && !is.list(points))
        points <- list(points)
    
    pointType <- match.arg(pointType)
    if (pointType == "general" && (!is.matrix(points) || ncol(points) != nDims))
        report(OL$Error, "Points must be specified as a matrix with #{nDims} columns")
    else if (pointType == "grid" && (!is.list(points) || length(points) != nDims))
        report(OL$Error, "Points must be specified as a list of length #{nDims}")
    else if (pointType == "auto")
    {
        if (is.matrix(points) && ncol(points) == nDims)
            pointType <- "general"
        else if (is.list(points) && length(points) == nDims)
            pointType <- "grid"
        else
            report(OL$Error, "Point specification is not valid")
    }
    
    if (is.matrix(points))
        points <- points - 1
    else if (is.list(points))
        points <- lapply(points, "-", 1)
    
    if (!isKernel(kernel))
        report(OL$Error, "Specified kernel is invalid")
    
    result <- .Call("resample", x, kernel, list(type=pointType,points=points), PACKAGE="mmand")
    
    if (is.list(points) && nDims > 1)
        dim(result) <- sapply(points, length)
    
    return (result)
}

rescale <- function (x, factor, kernel, ...)
{
    x <- as.array(x)
    dims <- dim(x)
    nDims <- length(dims)
    
    if (length(factor) < nDims)
        factor <- rep(factor, length.out=nDims)
    
    points <- lapply(seq_len(nDims), function(i) {
        newLength <- ceiling(dims[i] * factor[i])
        locs <- seq(0.5, dims[i]+0.5, length.out=newLength+1)
        locs <- locs + diff(locs[1:2]) / 2
        locs <- locs[1:newLength]
    })
    
    resample(x, points, kernel, ...)
}

neighbourhood <- function (x, width)
{
    return (.Call("get_neighbourhood", as.array(x), as.integer(width[1]), PACKAGE="mmand"))
}
