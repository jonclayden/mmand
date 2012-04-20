morph <- function (x, ...)
{
    UseMethod("morph")
}

morph.default <- function (x, kernel, brush = TRUE, value = NULL, valueNot = NULL, nNeighbours = NULL, nNeighboursNot = NULL)
{
    x <- as.array(x)
    if (!is.numeric(x))
        report(OL$Error, "Initial array must be numeric")
    
    kernel <- as.array(kernel)
    if (!is.numeric(kernel))
        report(OL$Error, "Kernel must be numeric")
    
    if (length(dim(kernel)) < length(dim(x)))
        dim(kernel) <- c(dim(kernel), rep(1,length(dim(x))-length(dim(kernel))))
    
    if (storage.mode(kernel) != storage.mode(x))
    {
        report(OL$Info, "Converting kernel to \"", storage.mode(x), "\" mode to match target array")
        storage.mode(kernel) <- storage.mode(x)
    }
    
    returnValue <- .Call("morph_R", x, kernel, as.double(value), as.double(valueNot), as.integer(nNeighbours), as.integer(nNeighboursNot), as.logical(brush), PACKAGE="mmand")
    
    return (returnValue)
}

binarise <- function (x)
{
    return (morph(x, kernel=1, brush=TRUE, valueNot=0))
}

erode <- function (x, kernel)
{
    return (morph(x, kernel, brush=TRUE, value=0))
}

dilate <- function (x, kernel)
{
    return (morph(x, kernel, brush=TRUE, valueNot=0))
}

opening <- function (x, kernel)
{
    return (dilate(erode(x, kernel), kernel))
}

closing <- function (x, kernel)
{
    return (erode(dilate(x, kernel), kernel))
}
