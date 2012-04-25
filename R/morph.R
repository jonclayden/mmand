morph <- function (x, ...)
{
    UseMethod("morph")
}

morph.default <- function (x, kernel, brush = TRUE, eraser = FALSE, value = NULL, valueNot = NULL, nNeighbours = NULL, nNeighboursNot = NULL)
{
    x <- as.array(x)
    if (!is.numeric(x))
        report(OL$Error, "Target array must be numeric")
    
    kernel <- as.array(kernel)
    if (!is.numeric(kernel))
        report(OL$Error, "Kernel must be numeric")
    
    if (!all(dim(kernel) == dim(kernel)[1]))
        report(OL$Error, "Kernel must have the same size in all dimensions")
    
    if (length(dim(kernel)) < length(dim(x)))
        dim(kernel) <- c(dim(kernel), rep(1,length(dim(x))-length(dim(kernel))))
    else if (length(dim(kernel)) > length(dim(x)))
        report(OL$Error, "Kernel has greater dimensionality than the target array")
    
    if (storage.mode(kernel) == "integer" && storage.mode(x) == "double")
    {
        report(OL$Verbose, "Converting kernel to \"double\" mode to match target array")
        storage.mode(kernel) <- "double"
    }
    else if (storage.mode(kernel) == "double" && storage.mode(x) == "integer")
    {
        kernelCopy <- kernel
        storage.mode(kernelCopy) <- "integer"
        if (isTRUE(all.equal(kernel, kernelCopy)))
        {
            report(OL$Verbose, "Converting kernel to \"integer\" mode to match target array")
            storage.mode(kernel) <- "integer"
        }
        else
        {
            # Kernel cannot be accurately represented in integer mode, so result won't be able to either
            # Hence, we need to modify the storage mode of the target array
            report(OL$Verbose, "Converting target array to \"double\" mode to match kernel")
            storage.mode(x) <- "double"
        }
    }
    
    returnValue <- .Call("morph_R", x, kernel, as.double(value), as.double(valueNot), as.integer(nNeighbours), as.integer(nNeighboursNot), as.logical(brush), as.logical(eraser), PACKAGE="mmand")
    dim(returnValue) <- dim(x)
    
    return (returnValue)
}

binarise <- function (x)
{
    kernel <- 1
    storage.mode(kernel) <- storage.mode(x)
    return (morph(x, kernel=kernel, brush=TRUE, valueNot=0))
}

gaussianSmooth <- function (x, sigma)
{
    kernel <- gaussianKernel(sigma, normalised=TRUE)
    return (morph(x, kernel, brush=FALSE))
}

erode <- function (x, kernel)
{
    if (is.array(kernel) && all(dim(kernel) == 3))
        return (morph(x, kernel, brush=TRUE, eraser=TRUE, value=0, nNeighboursNot=0))
    else
        return (morph(x, kernel, brush=TRUE, eraser=TRUE, value=0))
}

dilate <- function (x, kernel)
{
    if (is.array(kernel) && all(dim(kernel) == 3))
    {
        neighbourCount <- 3^length(dim(x)) - 1
        return (morph(x, kernel, brush=TRUE, valueNot=0, nNeighboursNot=neighbourCount))
    }
    else
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
