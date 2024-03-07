#' Resample an array
#' 
#' The \code{resample} function uses a kernel function to resample a target
#' array. This can be thought of as a generalisation of array indexing which
#' allows fractional indices. It is (S3) generic. The \code{rescale} function
#' is an alternative interface for the common case where the image is being
#' scaled to a new size.
#' 
#' @param x Any object. For the default method, this must be coercible to an
#'   array.
#' @param points Either a matrix giving the points to sample at, one per row,
#'   or a list giving the locations on each axis, which will be made into a grid.
#' @param kernel A kernel function object, used to provide coefficients for
#'   each resampled value, or the name of one.
#' @param pointType A string giving the type of the point specification being
#'   used. Usually can be left as \code{"auto"}.
#' @param threads If a positive integer, and the package is compiled with
#'   OpenMP support, the number of threads to use during the calculation.
#' @param factor A vector of scale factors, which will be recycled to the
#'   dimensionality of \code{x}.
#' @param \dots Additional options, such as kernel parameters.
#' @return If a generalised sampling scheme is used (i.e. with \code{points} a
#'   matrix), the result is a vector of sampled values. For a grid scheme (i.e.
#'   with \code{points} a list, including for \code{rescale}), it is a
#'   resampled array.
#' 
#' @examples
#' resample(c(0,0,1,0,0), seq(0.75,5.25,0.5), triangleKernel())
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{kernels}} for kernel-generating functions.
#' @export
resample <- function (x, points, kernel, ...)
{
    UseMethod("resample")
}

#' @rdname resample
#' @export
resample.default <- function (x, points, kernel, pointType = c("auto","general","grid"), threads = getOption("mmand.threads"), ...)
{
    x <- as.array(x)
    if (!is.numeric(x) && !is.logical(x))
        stop("Target array must be numeric")
    
    if (!isKernelFunction(kernel))
        kernel <- kernelFunction(kernel, ...)
    
    nDims <- length(dim(x))
    
    if (nDims == 1 && !is.matrix(points) && !is.list(points))
        points <- list(points)
    
    pointType <- match.arg(pointType)
    if (pointType == "general" && (!is.matrix(points) || ncol(points) != nDims))
        stop("Points must be specified as a matrix with #{nDims} columns")
    else if (pointType == "grid" && (!is.list(points) || length(points) != nDims))
        stop("Points must be specified as a list of length #{nDims}")
    else if (pointType == "auto")
    {
        if (is.matrix(points) && ncol(points) == nDims)
            pointType <- "general"
        else if (is.list(points) && length(points) == nDims)
            pointType <- "grid"
        else
            stop("Point specification is not valid")
    }
    
    if (is.matrix(points))
        points <- points - 1
    else if (is.list(points))
        points <- lapply(points, "-", 1)
    
    result <- .Call(C_resample, x, kernel, list(type=pointType,points=points), threads)
    
    if (is.list(points) && nDims > 1)
        dim(result) <- sapply(points, length)
    
    return (result)
}

#' @rdname resample
#' @export
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

#' Get neighbourhood information for an array
#' 
#' This function provides information about the structure of a neighbourhood of
#' a given width within a specified array.
#' 
#' @param x An object that can be coerced to an array.
#' @param width A vector giving the width of the neighbourhood in each
#'   dimension, which will be recycled if necessary. Must not be greater than
#'   the size of the array. Even values are rounded up to the next odd integer.
#' @return A list with the following elements.
#'     \item{widths}{The width of the neighbourhood along each dimension.
#'       Currently all elements of this vector will be the same.}
#'     \item{size}{The number of pixels within the neighbourhood.}
#'     \item{locs}{A matrix giving the coordinates of each neighbourhood pixel
#'       relative to the centre pixel, one per row.}
#'     \item{offsets}{Vector offsets of the neighbourhood values within
#'       \code{x}.}
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @export
neighbourhood <- function (x, width)
{
    x <- as.array(x)
    
    nDims <- length(dim(x))
    if (length(width) < nDims)
        width <- rep(width, length.out=nDims)
    if (any(width > dim(x)))
        stop("Requested neighbourhood is larger than the data")
    
    return (.Call(C_get_neighbourhood, x, as.integer(width)))
}
