#' Morph an array with a kernel
#' 
#' The \code{morph} function applies a kernel to a target array. Optionally,
#' applying the kernel to a particular array element can be made conditional on
#' its value, or the number of nonzero immediate neighbours that it has. The
#' \code{morph} function is (S3) generic.
#' 
#' 
#' @aliases morph morph.default
#' @param x Any object. For the default method, this must be coercible to an
#' array.
#' @param kernel An object representing the kernel to be applied, which must be
#' coercible to an array. It must have odd width in all dimensions, but does
#' not have to be isotropic in size. The kernel's dimensionality may be less
#' than that of the target array, \code{x}. See \code{\link{kernels}} for
#' kernel-generating functions.
#' @param operator The operator applied elementwise within the kernel, as a
#' function of the original image value and the kernel value. Arithmetic
#' operators are as usual; \code{"i"} is the identity operator, where every
#' value within the kernel will be included as-is; and \code{1} and \code{0}
#' include a 1 or 0 for each element within the kernel's nonzero region.
#' @param merge The operator applied to combine the elements into a final value
#' for the centre pixel. All have their usual meanings.
#' @param value An optional vector of values in the target array for which to
#' apply the kernel. Takes priority over \code{valueNot} if both are specified.
#' @param valueNot An optional vector of values in the target array for which
#' not to apply the kernel.
#' @param nNeighbours An optional numeric vector giving allowable numbers of
#' nonzero neighbours (including diagonal neighbours) for array elements where
#' the kernel will be applied. Takes priority over \code{nNeighboursNot} if
#' both are specified.
#' @param nNeighboursNot An optional numeric vector giving nonallowable numbers
#' of nonzero neighbours (including diagonal neighbours) for array elements
#' where the kernel will be applied.
#' @param \dots Additional arguments to methods.
#' @return A morphed array with the same dimensions as the original array.
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{kernels}} for kernel-generating functions, and
#' \code{\link{morphology}} for more specific mathematical morphology
#' functions. \code{\link{gameOfLife}} shows how this function can be used for
#' non-morphological purposes, in that case to power a cellular automaton. See
#' also the \code{kernel} and \code{kernapply} functions in the \code{stats}
#' package, particularly if you want to smooth time series.
morph <- function (x, kernel, ...)
{
    UseMethod("morph")
}

morph.default <- function (x, kernel, operator = c("+","-","*","i","1","0"), merge = c("sum","min","max","mean","median"), value = NULL, valueNot = NULL, nNeighbours = NULL, nNeighboursNot = NULL, ...)
{
    x <- as.array(x)
    if (!is.numeric(x))
        report(OL$Error, "Target array must be numeric")
    
    if (!isKernelArray(kernel))
        kernel <- kernelArray(kernel)
    
    if (any(dim(kernel) %% 2 != 1))
        report(OL$Error, "Kernel must have odd width in all dimensions")
    
    if (length(dim(kernel)) < length(dim(x)))
        dim(kernel) <- c(dim(kernel), rep(1,length(dim(x))-length(dim(kernel))))
    else if (length(dim(kernel)) > length(dim(x)))
        report(OL$Error, "Kernel has greater dimensionality than the target array")
    
    operator <- match.arg(operator)
    merge <- match.arg(merge)
    
    storage.mode(x) <- "double"
    
    restrictions <- list(value=as.double(value), valueNot=as.double(valueNot), nNeighbours=as.integer(nNeighbours), nNeighboursNot=as.integer(nNeighboursNot))
    
    returnValue <- .Call("morph", x, kernel, operator, merge, restrictions, PACKAGE="mmand")
    
    if (length(dim(x)) > 1)
        dim(returnValue) <- dim(x)
    
    return (returnValue)
}


#' Check for a binary array
#' 
#' This function checks whether a numeric array is binary, with only one unique
#' nonzero value, or not.
#' 
#' 
#' @param x An object that can be coerced to a numeric array
#' @return A logical value indicating whether the array is binary or not.
#' Binary in this case means that the array contains only one unique nonzero
#' value.
#' @author Jon Clayden <code@@clayden.org>
binary <- function (x)
{
    x <- as.array(x)
    if (!is.numeric(x))
        report(OL$Error, "Array must be numeric")
    
    return (.Call("is_binary", x, PACKAGE="mmand"))
}


#' Binarise a numeric array
#' 
#' This function binarises an array, setting all nonzero elements to unity.
#' 
#' 
#' @param x An object that can be coerced to an array, or for which a
#' \code{\link{morph}} method exists.
#' @return A morphed array with the same dimensions as the original array.
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{morph}} for the function underlying this operation, and
#' \code{\link{erode}} for mathematical morphology functions.
binarise <- function (x)
{
    return (morph(x, kernel=1, operator="1", valueNot=0))
}


#' Smooth a numeric array with a Gaussian kernel
#' 
#' This function smoothes an array using a Gaussian kernel with a specified
#' standard deviation.
#' 
#' 
#' @param x An object that can be coerced to an array, or for which a
#' \code{\link{morph}} method exists.
#' @param sigma A numeric vector giving the standard deviation of the kernel in
#' each dimension. Can have lower dimensionality than the target array.
#' @return A morphed array with the same dimensions as the original array.
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{morph}} for the function underlying this operation,
#' \code{\link{gaussianKernel}} for generating Gaussian kernels (which is also
#' used by this function), and \code{\link{erode}} for mathematical morphology
#' functions.
gaussianSmooth <- function (x, sigma)
{
    kernel <- gaussianKernel(sigma, normalised=TRUE)
    return (morph(x, kernel, operator="*", merge="sum"))
}


#' Apply a filter to an array
#' 
#' These functions apply mean or median filters to an array.
#' 
#' 
#' @aliases filters medianFilter meanFilter
#' @param x An object that can be coerced to an array, or for which a
#' \code{\link{morph}} method exists.
#' @param kernel A kernel array, indicating the scope of the filter.
#' @return A morphed array with the same dimensions as the original array.
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{morph}} for the function underlying this operation, and
#' \code{\link{kernels}} for kernel-generating functions.
meanFilter <- function (x, kernel)
{
    return (morph(x, kernel, operator="i", merge="mean"))
}

medianFilter <- function (x, kernel)
{
    return (morph(x, kernel, operator="i", merge="median"))
}


#' Standard mathematical morphology operations
#' 
#' These functions provide standard mathematical morphology operations, which
#' can be applied to array data with any number of dimensions. Binary and
#' greyscale morphology is supported.
#' 
#' The \code{erode} function uses the kernel as an eraser, centring it on each
#' zero-valued pixel, which has the effect of eroding the extent of nonzero
#' areas. Dilation has the opposite effect, extending the nonzero regions in
#' the array. Opening is an erosion followed by a dilation, and closing is a
#' dilation followed by an erosion, using the same kernel in both cases.
#' 
#' If the kernel has only one unique nonzero value, it is described as
#' ``flat''. For a flat kernel, the erosion is the minimum value of \code{x}
#' within the nonzero region of \code{kernel}. For a nonflat kernel, this
#' becomes minimum value of \code{x - kernel}. Dilation is the opposite
#' operation, taking the maximum within the kernel.
#' 
#' @aliases morphology erode dilate opening closing
#' @param x An object that can be coerced to an array, or for which a
#' \code{\link{morph}} method exists.
#' @param kernel An array representing the kernel to be used. See
#' \code{\link{shapeKernel}} for functions to generate a suitable kernel.
#' @return A morphed array with the same dimensions as the original array.
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{morph}} for the function underlying all of these
#' operations, \code{\link{kernels}} for kernel-generating functions,
#' \code{\link{binarise}} for binarising an array, and
#' \code{\link{gaussianSmooth}} for smoothing. The \code{EBImage} Bioconductor
#' package also supplies functions to perform these operations, and may be
#' slightly faster, but only works in two dimensions.
#' @examples
#' x <- c(0,0,1,0,0,0,1,1,1,0,0)
#' k <- c(1,1,1)
#' erode(x,k)
#' dilate(x,k)
erode <- function (x, kernel)
{
    x <- as.array(x)
    if (!isKernelArray(kernel))
        kernel <- kernelArray(kernel)
    
    greyscaleImage <- !binary(x)
    nNeighboursNot <- NULL
    
    if (greyscaleImage)
    {
        operator <- ifelse(binary(kernel), "i", "-")
        valueNot <- NULL
    }
    else
    {
        if (all(dim(kernel) <= 3))
            nNeighboursNot <- 3^length(dim(x)) - 1
        operator <- "i"
        valueNot <- 0
    }
    
    return (morph(x, kernel, operator=operator, merge="min", valueNot=valueNot, nNeighboursNot=nNeighboursNot))
}

dilate <- function (x, kernel)
{
    x <- as.array(x)
    if (!isKernelArray(kernel))
        kernel <- kernelArray(kernel)
    
    greyscaleImage <- !binary(x)
    nNeighboursNot <- NULL
    
    if (greyscaleImage)
    {
        operator <- ifelse(binary(kernel), "i", "+")
        value <- NULL
    }
    else
    {
        if (all(dim(kernel) <= 3))
            nNeighboursNot <- 0
        operator <- "i"
        value <- 0
    }
    
    return (morph(x, kernel, operator=operator, merge="max", value=value, nNeighboursNot=nNeighboursNot))
}

opening <- function (x, kernel)
{
    x <- as.array(x)
    if (!isKernelArray(kernel))
        kernel <- kernelArray(kernel)
    
    return (dilate(erode(x, kernel), kernel))
}

closing <- function (x, kernel)
{
    x <- as.array(x)
    if (!isKernelArray(kernel))
        kernel <- kernelArray(kernel)
    
    return (erode(dilate(x, kernel), kernel))
}
