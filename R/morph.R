morph <- function (x, kernel, ...)
{
    UseMethod("morph")
}

morph.default <- function (x, kernel, value = NULL, valueNot = NULL, nNeighbours = NULL, nNeighboursNot = NULL, ...)
{
    x <- as.array(x)
    if (!is.numeric(x))
        report(OL$Error, "Target array must be numeric")
    
    if (is.array(kernel))
        kernel <- discreteKernel(kernel)
    else if (!isKernel(kernel))
        output(OL$Error, "Specified kernel is invalid")
    
    if (any(dim(kernel$values) %% 2 != 1))
        report(OL$Error, "Kernel must have odd width in all dimensions")
    
    if (length(dim(kernel$values)) < length(dim(x)))
        dim(kernel$values) <- c(dim(kernel$values), rep(1,length(dim(x))-length(dim(kernel$values))))
    else if (length(dim(kernel$values)) > length(dim(x)))
        report(OL$Error, "Kernel has greater dimensionality than the target array")
    
    storage.mode(x) <- "double"
    
    restrictions = list(value=as.double(value), valueNot=as.double(valueNot), nNeighbours=as.integer(nNeighbours), nNeighboursNot=as.integer(nNeighboursNot))
    
    returnValue <- .Call("morph", x, kernel, restrictions, PACKAGE="mmand")
    dim(returnValue) <- dim(x)
    return (returnValue)
}

binarise <- function (x)
{
    kernel <- discreteKernel(1, brush=TRUE)
    return (morph(x, kernel=kernel, valueNot=0))
}

gaussianSmooth <- function (x, sigma)
{
    kernel <- gaussianKernel(sigma, normalised=TRUE)
    return (morph(x, kernel))
}

erode <- function (x, kernel)
{
    kernel <- discreteKernel(kernel, brush=TRUE, eraser=TRUE)
    
    if (all(dim(kernel$values) <= 3))
        return (morph(x, kernel, value=0, nNeighboursNot=0))
    else
        return (morph(x, kernel, value=0))
}

dilate <- function (x, kernel)
{
    kernel <- discreteKernel(kernel, brush=TRUE, eraser=FALSE)
    
    if (all(dim(kernel$values) <= 3))
    {
        neighbourCount <- 3^length(dim(x)) - 1
        return (morph(x, kernel, valueNot=0, nNeighboursNot=neighbourCount))
    }
    else
        return (morph(x, kernel, valueNot=0))
}

opening <- function (x, kernel)
{
    return (dilate(erode(x, kernel), kernel))
}

closing <- function (x, kernel)
{
    return (erode(dilate(x, kernel), kernel))
}
