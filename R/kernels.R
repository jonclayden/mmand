#' Kernel-generating functions
#' 
#' These functions can be used to generate kernels for morphological, smoothing
#' or resampling operations. There are two types of kernels: kernel arrays,
#' which are used with \code{\link{morph}}, and kernel functions, which are
#' used with \code{\link{resample}}.
#' 
#' There are two forms of kernel used by this package. Kernel arrays, otherwise
#' known in mathematical morphology as structuring elements, are numeric arrays
#' with class \code{kernelArray}. They are defined on a grid of odd width, and
#' are used by \code{\link{morph}} and related functions. Kernel functions, by
#' contrast, are represented in R as a list containing a name and, optionally,
#' some parameters. The real implementation is in C++. They are defined
#' everywhere within the support of the kernel, and are used by
#' \code{\link{resample}} and friends. The key distinction is in whether the
#' kernel will always be centred exactly on the location of an existing value
#' in the data (for kernel arrays) or not (for kernel functions).
#' 
#' The \code{kernelArray} and \code{kernelFunction} functions create objects of
#' the corresponding classes, while \code{isKernelArray} and
#' \code{isKernelFunction} test for them. In addition, \code{isKernel} returns
#' \code{TRUE} if its argument is of either kernel class.
#' 
#' The remaining functions generate special-case kernels: \code{shapeKernel}
#' generates arrays with nonzero elements in a box, disc or diamond shape for
#' use with \code{\link{morphology}} functions; \code{gaussianKernel} generates
#' Gaussian coefficients and is used by \code{\link{gaussianSmooth}};
#' \code{sobelKernel} generates the Sobel-Feldman gradient operator, for use by
#' \code{\link{sobelFilter}}; \code{boxKernel} is used for ``nearest
#' neighbour'' resampling, and \code{triangleKernel} for linear, bilinear, etc.
#' The Mitchell-Netravali kernel, a.k.a. BC-spline, is based on a family of
#' piecewise-cubic polynomial functions, with support of four times the pixel
#' separation in each dimension. The default parameters are the ones
#' recommended by Mitchell and Netravali as a good trade-off between various
#' artefacts, but other well-known special cases include B=1, C=0 (the cubic
#' B-spline) and B=0, C=0.5 (the Catmull-Rom spline). \code{mnKernel} is a
#' shorter alias for \code{mitchellNetravaliKernel}.
#' 
#' @param object Any object.
#' @param values A numeric vector or array, containing the values of the kernel
#'   array.
#' @param width A numeric vector giving the width of the shape in each
#'   dimension, in array elements. Does not need to be integer-valued, or equal
#'   for all dimensions. Will be recycled to length \code{dim} if that
#'   parameter is also specified.
#' @param sigma A numeric vector giving the standard deviation of the
#'   underlying Gaussian distribution in each dimension, in array elements.
#'   Does not need to be equal for all dimensions. Will be recycled to length
#'   \code{dim} if that parameter is also specified.
#' @param dim An integer value giving the dimensionality of the kernel.
#'   Defaults to the length of \code{width} or \code{sigma}, where available.
#' @param size A numeric vector giving the width of the kernel in each
#'   dimension, which will be rounded up to the nearest odd integer. Defaults
#'   to four times the corresponding \code{sigma} value.
#' @param type A string giving the type of shape to produce. In one dimension,
#'   these shapes are all equivalent.
#' @param binary If \code{FALSE}, the value of the kernel at each point
#'   represents the proportion of the array element within the shape. If
#'   \code{TRUE}, these values are binarised to be 1 if at least half of the
#'   element is within the shape, and 0 otherwise.
#' @param normalised If \code{TRUE}, the sum of non-missing elements of the
#'   kernel will be unity. Note that this is the default for
#'   \code{gaussianKernel}, but not for \code{shapeKernel}.
#' @param axis The axis along which the gradient operator will be applied.
#' @param name A string giving the name of the kernel function required.
#' @param \dots Parameters for the kernel function.
#' @param B,C Mitchell-Netravali coefficients, each of which must be between 0
#'   and 1.
#' @return For \code{isKernel}, \code{isKernelArray} and
#'   \code{isKernelFunction}, a logical value. For \code{kernelArray},
#'   \code{shapeKernel}, \code{gaussianKernel} and \code{sobelKernel}, a kernel
#'   array. For \code{kernelFunction}, \code{boxKernel}, \code{triangleKernel},
#'   \code{mitchellNetravaliKernel} and \code{mnKernel}, a kernel function.
#' 
#' @examples
#' shapeKernel(c(3,5), type="diamond")
#' gaussianKernel(c(0.3,0.3))
#' mnKernel()
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{morph}} for general application of kernel arrays to
#'   data, \code{\link{morphology}} for mathematical morphology functions,
#'   \code{\link{resample}} for resampling, and \code{\link{gaussianSmooth}}
#'   for smoothing. Also see \code{\link{sampleKernelFunction}} for kernel
#'   sampling and plotting.
#' @references The Mitchell-Netravali kernel is described in the following
#'   paper.
#'   
#'   D.P. Mitchell & A.N. Netravali (1988). Reconstruction filters in computer
#'   graphics. Computer Graphics 22(4):221-228.
#' @rdname kernels
#' @aliases kernels
#' @export
isKernel <- function (object)
{
    return ("kernel" %in% class(object))
}

#' @rdname kernels
#' @export
isKernelArray <- function (object)
{
    return ("kernelArray" %in% class(object))
}

#' @rdname kernels
#' @export
isKernelFunction <- function (object)
{
    return ("kernelFunction" %in% class(object))
}

#' Sampling and plotting kernels
#' 
#' These functions can be used to sample and plot kernel profiles.
#' 
#' @param kernel A kernel function object.
#' @param values A vector of values to sample the function at. These are in
#'   units of pixels, with zero representing the centre of the kernel.
#' @param x A kernel object of the appropriate class.
#' @param y Ignored.
#' @param xlim The limits of the range used to profile the kernel.
#' @param lwd The line width to use for the kernel profile.
#' @param col The line colour to use for the kernel profile.
#' @param axis The axis to profile along.
#' @param \dots Additional plot parameters.
#' @return For \code{sampleKernelFunction} a vector of kernel values at the
#'   locations requested. The \code{plot} methods are called for their
#'   side-effects.
#' 
#' @examples
#' 
#' sampleKernelFunction(mnKernel(), -2:2)
#' plot(mnKernel())
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{kernels}} for kernel-generating functions.
#' @export 
sampleKernelFunction <- function (kernel, values)
{
    if (!isKernelFunction(kernel))
        report(OL$Error, "Specified kernel is not a valid kernel function")
    
    return (.Call("sample_kernel", kernel, as.numeric(values), PACKAGE="mmand"))
}

#' @rdname sampleKernelFunction
#' @export
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

#' @rdname sampleKernelFunction
#' @export
plot.kernelFunction <- function (x, y, xlim = c(-2,2), lwd = 2, col = "red", ...)
{
    values <- seq(xlim[1], xlim[2], length.out=101)
    plot(values, sampleKernelFunction(x,values), xlab="x", ylab="k(x)", type="l", lwd=lwd, col=col, xlim=xlim, ...)
    abline(v=0, lty=2, col="grey60")
    abline(h=0, lty=2, col="grey60")
}

#' @rdname kernels
#' @export
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

#' @rdname kernels
#' @export
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

#' @rdname kernels
#' @export
gaussianKernel <- function (sigma, dim = length(sigma), size = 6*sigma, normalised = TRUE)
{
    if (dim > length(sigma))
        sigma <- rep(sigma, length.out=dim)
    
    if (dim > length(size))
        size <- rep(ceiling(size), length.out=dim)
    else
        size <- ceiling(size)
    
    size <- ifelse(size %% 2 == 1, size, size+1)
    
    scaleFactors <- max(sigma) / ifelse(sigma==0, 1, sigma)
    x <- lapply(1:dim, function(i) 1:size[i] - (size[i]+1)/2)
    centres <- lapply(1:dim, function(i) scaleFactors[i] * x[[i]])
    
    normFun <- function(a,b) sqrt(a^2 + b^2)
    norms <- Reduce(function(a,b) outer(a,b,FUN=normFun), centres)
    
    kernel <- array(dnorm(norms,sd=max(sigma)), dim=size)
    
    if (normalised)
        kernel <- kernel / sum(kernel, na.rm=TRUE)
    
    return (kernelArray(kernel))
}

#' @rdname kernels
#' @export
sobelKernel <- function (dim, axis = 1)
{
    if (!(axis %in% 0:dim))
        report(OL$Error, "The axis should be between 0 and #{dim}")
    
    parts <- rep(list(c(1,2,1) / 4), dim)
    if (axis > 0)
        parts[[axis]] <- c(1, 0, -1)
    
    kernel <- Reduce(outer, parts)
    return (kernelArray(kernel))
}

#' @rdname kernels
#' @export
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

#' @rdname kernels
#' @export
boxKernel <- function ()
{
    return (kernelFunction("box"))
}

#' @rdname kernels
#' @export
triangleKernel <- function ()
{
    return (kernelFunction("triangle"))
}

#' @rdname kernels
#' @export
mitchellNetravaliKernel <- function (B = 1/3, C = 1/3)
{
    return (kernelFunction("mitchell-netravali", B=B, C=C))
}

#' @rdname kernels
#' @export
mnKernel <- mitchellNetravaliKernel
