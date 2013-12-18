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
        output(OL$Error, "Points must be specified as a matrix with #{nDims} columns")
    else if (pointType == "grid" && (!is.list(points) || length(points) != nDims))
        output(OL$Error, "Points must be specified as a list of length #{nDims}")
    else if (pointType == "auto")
    {
        if (is.matrix(points) && ncol(points) == nDims)
            pointType <- "general"
        else if (is.list(points) && length(points) == nDims)
            pointType <- "grid"
        else
            output(OL$Error, "Point specification is not valid")
    }
    
    if (is.matrix(points))
        points <- points - 1
    else if (is.list(points))
        points <- lapply(points, "-", 1)
    
    if (!is.list(kernel) || !("kernel" %in% class(kernel)))
        output(OL$Error, "Specified kernel is invalid")
    
    result <- .Call("resample", x, kernel, list(type=pointType,points=points), PACKAGE="irk")
    
    if (is.list(points) && nDims > 1)
        dim(result) <- sapply(points, length)
    
    return (result)
}

neighbourhood <- function (x, width)
{
    return (.Call("get_neighbourhood", as.array(x), as.integer(width[1]), PACKAGE="irk"))
}
