morph <- function (x, kernel, ...)
{
    UseMethod("morph")
}

morph.default <- function (x, kernel, operator = c("+","-","*","i","1","0"), merge = c("sum","min","max","mean","median"), value = NULL, valueNot = NULL, nNeighbours = NULL, nNeighboursNot = NULL, ...)
{
    x <- as.array(x)
    if (!is.numeric(x))
        report(OL$Error, "Target array must be numeric")
    
    kernel <- as.array(kernel)
    if (!is.numeric(kernel))
        report(OL$Error, "Kernel must be numeric")
    
    operator <- match.arg(operator)
    merge <- match.arg(merge)
    
    if (any(dim(kernel) %% 2 != 1))
        report(OL$Error, "Kernel must have odd width in all dimensions")
    
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
    
    returnValue <- .Call("morph_R", x, kernel, as.double(value), as.double(valueNot), as.integer(nNeighbours), as.integer(nNeighboursNot), operator, merge, PACKAGE="mmand")
    dim(returnValue) <- dim(x)
    
    return (returnValue)
}

binarise <- function (x)
{
    kernel <- as(1, storage.mode(x))
    return (morph(x, kernel=kernel, operator="1", valueNot=0))
}

gaussianSmooth <- function (x, sigma)
{
    kernel <- gaussianKernel(sigma, normalised=TRUE)
    return (morph(x, kernel, operator="*", merge="sum"))
}

meanFilter <- function (x, kernel)
{
    return (morph(x, kernel, operator="i", merge="mean"))
}

medianFilter <- function (x, kernel)
{
    return (morph(x, kernel, operator="i", merge="median"))
}

erode <- function (x, kernel, greyscale = FALSE)
{
    operator <- ifelse(greyscale, "-", "i")
    valueNot <- (if (greyscale) NULL else 0)
    
    if (is.array(x) && is.array(kernel) && all(dim(kernel) <= 3))
    {
        nNeighboursNot <- (if (greyscale) NULL else 3^length(dim(x))-1)
        return (morph(x, kernel, operator=operator, merge="min", valueNot=valueNot, nNeighboursNot=nNeighboursNot))
    }
    else
        return (morph(x, kernel, operator=operator, merge="min", valueNot=valueNot))
}

dilate <- function (x, kernel, greyscale = FALSE)
{
    operator <- ifelse(greyscale, "+", "i")
    value <- nNeighboursNot <- (if (greyscale) NULL else 0)
    
    if (is.array(kernel) && all(dim(kernel) <= 3))
        return (morph(x, kernel, operator=operator, merge="max", value=value, nNeighboursNot=nNeighboursNot))
    else
        return (morph(x, kernel, operator=operator, merge="max", value=value))
}

opening <- function (x, kernel, greyscale = FALSE)
{
    return (dilate(erode(x, kernel, greyscale), kernel, greyscale))
}

closing <- function (x, kernel, greyscale = FALSE)
{
    return (erode(dilate(x, kernel, greyscale), kernel, greyscale))
}
