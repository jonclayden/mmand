morph <- function (x, kernel, ...)
{
    UseMethod("morph")
}

morph.default <- function (x, kernel, operator = c("+","-","*","i","1","0"), merge = c("sum","min","max","mean","median"), value = NULL, valueNot = NULL, nNeighbours = NULL, nNeighboursNot = NULL, ...)
{
    x <- as.array(x)
    if (!is.numeric(x))
        report(OL$Error, "Target array must be numeric")
    
    if (!isKernel(kernel))
        kernel <- discreteKernel(kernel)
    
    if (any(dim(kernel$values) %% 2 != 1))
        report(OL$Error, "Kernel must have odd width in all dimensions")
    
    if (length(dim(kernel$values)) < length(dim(x)))
        dim(kernel$values) <- c(dim(kernel$values), rep(1,length(dim(x))-length(dim(kernel$values))))
    else if (length(dim(kernel$values)) > length(dim(x)))
        report(OL$Error, "Kernel has greater dimensionality than the target array")
    
    operator <- match.arg(operator)
    merge <- match.arg(merge)
    
    storage.mode(x) <- "double"
    
    restrictions <- list(value=as.double(value), valueNot=as.double(valueNot), nNeighbours=as.integer(nNeighbours), nNeighboursNot=as.integer(nNeighboursNot))
    
    returnValue <- .Call("morph", x, kernel, operator, merge, restrictions, PACKAGE="mmand")
    dim(returnValue) <- dim(x)
    return (returnValue)
}

binary <- function (x)
{
    x <- as.array(x)
    if (!is.numeric(x))
        report(OL$Error, "Array must be numeric")
    
    return (.Call("is_binary", x, PACKAGE="mmand"))
}

binarise <- function (x)
{
    return (morph(x, kernel=1, operator="1", valueNot=0))
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

erode <- function (x, kernel)
{
    greyscaleImage <- !binary(x)
    if (greyscaleImage)
    {
        flatKernel <- binary(kernel)
        operator <- ifelse(!flatKernel, "-", "i")
        valueNot <- NULL
    }
    else
    {
        flatKernel <- TRUE
        operator <- "i"
        valueNot <- 0
    }
    
    if (is.array(x) && is.array(kernel) && all(dim(kernel) <= 3))
    {
        nNeighboursNot <- (if (greyscaleImage) NULL else 3^length(dim(x))-1)
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

opening <- function (x, kernel)
{
    return (dilate(erode(x, kernel), kernel))
}

closing <- function (x, kernel)
{
    return (erode(dilate(x, kernel), kernel))
}
