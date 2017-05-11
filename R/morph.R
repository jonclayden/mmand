#' Morph an array with a kernel
#' 
#' The \code{morph} function applies a kernel to a target array. Optionally,
#' applying the kernel to a particular array element can be made conditional on
#' its value, or the number of nonzero immediate neighbours that it has. The
#' \code{morph} function is (S3) generic.
#' 
#' @param x Any object. For the default method, this must be coercible to an
#'   array.
#' @param kernel An object representing the kernel to be applied, which must be
#'   coercible to an array. It must have odd width in all dimensions, but does
#'   not have to be isotropic in size. The kernel's dimensionality may be less
#'   than that of the target array, \code{x}. See \code{\link{kernels}} for
#'   kernel-generating functions.
#' @param operator The operator applied elementwise within the kernel, as a
#'   function of the original image value and the kernel value. Arithmetic
#'   operators are as usual; \code{"i"} is the identity operator, where every
#'   value within the kernel will be included as-is; \code{"1"} and \code{"0"}
#'   include a 1 or 0 for each element within the kernel's nonzero region;
#'   \code{"=="} produces a 1 where the image matches the kernel, and 0
#'   elsewhere.
#' @param merge The operator applied to combine the elements into a final value
#'   for the centre pixel. All have their usual meanings.
#' @param value An optional vector of values in the target array for which to
#'   apply the kernel. Takes priority over \code{valueNot} if both are
#'   specified.
#' @param valueNot An optional vector of values in the target array for which
#'   not to apply the kernel.
#' @param nNeighbours An optional numeric vector giving allowable numbers of
#'   nonzero neighbours (including diagonal neighbours) for array elements
#'   where the kernel will be applied. Takes priority over
#'   \code{nNeighboursNot} if both are specified.
#' @param nNeighboursNot An optional numeric vector giving nonallowable numbers
#'   of nonzero neighbours (including diagonal neighbours) for array elements
#'   where the kernel will be applied.
#' @param renormalise If \code{TRUE}, the default, and \code{merge} is
#'   \code{"sum"}, the sum will be renormalised relative to the sum over the
#'   visited part of the kernel. This avoids low-intensity bands around the
#'   edges of a morphed image.
#' @param \dots Additional arguments to methods.
#' @return A morphed array with the same dimensions as the original array.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{kernels}} for kernel-generating functions, and
#'   \code{\link{morphology}} for more specific mathematical morphology
#'   functions. \code{\link{gameOfLife}} shows how this function can be used
#'   for non-morphological purposes, in that case to power a cellular
#'   automaton. See also the \code{kernel} and \code{kernapply} functions in
#'   the \code{stats} package, particularly if you want to smooth time series.
#' @export
morph <- function (x, kernel, ...)
{
    UseMethod("morph")
}

#' @rdname morph
#' @export
morph.default <- function (x, kernel, operator = c("+","-","*","i","1","0","=="), merge = c("sum","min","max","mean","median","all","any"), value = NULL, valueNot = NULL, nNeighbours = NULL, nNeighboursNot = NULL, renormalise = TRUE, ...)
{
    x <- as.array(x)
    if (!is.numeric(x) && !is.logical(x))
        stop("Target array must be numeric")
    
    if (!isKernelArray(kernel))
        kernel <- kernelArray(kernel)
    
    if (any(dim(kernel) %% 2 != 1))
        stop("Kernel must have odd width in all dimensions")
    
    if (length(dim(kernel)) < length(dim(x)))
        dim(kernel) <- c(dim(kernel), rep(1,length(dim(x))-length(dim(kernel))))
    else if (length(dim(kernel)) > length(dim(x)))
        stop("Kernel has greater dimensionality than the target array")
    
    operator <- match.arg(operator)
    merge <- match.arg(merge)
    
    storage.mode(x) <- "double"
    
    restrictions <- list(value=as.double(value), valueNot=as.double(valueNot), nNeighbours=as.integer(nNeighbours), nNeighboursNot=as.integer(nNeighboursNot))
    
    returnValue <- .Call(C_morph, x, kernel, operator, merge, restrictions, renormalise)
    
    if (length(dim(x)) > 1)
        dim(returnValue) <- dim(x)
    
    return (returnValue)
}

#' Check for a binary array
#' 
#' This function checks whether a numeric array is binary, with only one unique
#' nonzero value, or not.
#' 
#' @param x An object that can be coerced to a numeric array.
#' @return A logical value indicating whether the array is binary or not.
#'   Binary in this case means that the array contains only one unique nonzero
#'   value, which is stored with the return value in an attribute.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @export
binary <- function (x)
{
    x <- as.array(x)
    if (!is.numeric(x) && !is.logical(x))
        stop("Array must be numeric")
    
    return (.Call(C_is_binary, x))
}

#' Binarise a numeric array
#' 
#' This function binarises an array, setting all nonzero elements to unity.
#' 
#' @param x An object that can be coerced to an array, or for which a
#'   \code{\link{morph}} method exists.
#' @return A morphed array with the same dimensions as the original array.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{morph}} for the function underlying this operation, and
#'   \code{\link{erode}} for mathematical morphology functions.
#' @aliases binarize
#' @export binarise binarize
binarise <- binarize <- function (x)
{
    return (morph(x, kernel=1, operator="1", valueNot=0))
}

#' Threshold a numeric array or vector
#' 
#' This function thresholds an array or vector, setting elements below the
#' threshold value to zero. The threshold can be given literally or calculated
#' using k-means clustering.
#' 
#' @param x A numeric vector or array.
#' @param level The literal threshold level, if required.
#' @param method The method to use to calculate the threshold. If
#'   \code{"literal"} (the default) then the value of \code{level} will be
#'   used. If \code{"kmeans"} then the threshold value will be determined
#'   implicitly using k-means clustering.
#' @param binarise Whether to set suprathreshold elements to unity (if
#'   \code{TRUE}), or leave them at their original values (if \code{FALSE}).
#' 
#' @examples
#' x <- c(0.1, 0.05, 0.95, 0.85, 0.15, 0.9)
#' threshold(x, method="kmeans")
#' threshold(x, 0.5)
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{binarise}}
#' @export
threshold <- function (x, level, method = c("literal","kmeans"), binarise = TRUE)
{
    method <- match.arg(method)
    if (missing(level) && method == "literal")
        stop("A literal threshold level is required")
    
    if (method == "literal")
    {
        if (binarise)
            x <- ifelse(x < level, 0L, 1L)
        else
            x <- ifelse(x < level, 0L, x)
    }
    else if (method == "kmeans")
    {
        # Drop dims so that the clustering is done by intensity only
        dims <- dim(x)
        dim(x) <- NULL
        valid <- is.finite(x)
        
        kmeansResult <- kmeans(x[valid], 2)
        if (binarise)
            x[valid] <- ifelse(kmeansResult$cluster == which.min(kmeansResult$centers), 0L, 1L)
        else
            x[valid] <- ifelse(kmeansResult$cluster == which.min(kmeansResult$centers), 0L, x[valid])
        x[!valid] <- NA
        
        dim(x) <- dims
    }
    
    return (x)
}

#' Smooth a numeric array with a Gaussian kernel
#' 
#' This function smoothes an array using a Gaussian kernel with a specified
#' standard deviation.
#' 
#' This implementation takes advantage of the separability of the Gaussian
#' kernel for speed when working in multiple dimensions. It is therefore
#' equivalent to, but much faster than, directly applying a multidimensional
#' kernel.
#' 
#' @param x An object that can be coerced to an array, or for which a
#'   \code{\link{morph}} method exists.
#' @param sigma A numeric vector giving the standard deviation of the kernel in
#'   each dimension. Can have lower dimensionality than the target array.
#' @return A morphed array with the same dimensions as the original array.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{morph}} for the function underlying this operation,
#'   \code{\link{gaussianKernel}} for generating Gaussian kernels (which is
#'   also used by this function), and \code{\link{erode}} for mathematical
#'   morphology functions.
#' @export
gaussianSmooth <- function (x, sigma)
{
    morphFun <- function(y,k) morph(y, k, operator="*", merge="sum")
    
    kernels <- lapply(seq_along(sigma), function(i) {
        currentSigma <- replace(rep(0,length(sigma)), i, sigma[i])
        gaussianKernel(currentSigma, normalised=TRUE)
    })
    
    return (Reduce(morphFun, kernels, x))
}

#' Apply a filter to an array
#' 
#' These functions apply mean, median or Sobel filters to an array.
#' 
#' @param x An object that can be coerced to an array, or for which a
#'   \code{\link{morph}} method exists.
#' @param kernel A kernel array, indicating the scope of the filter.
#' @param dim For \code{sobelFilter}, the dimensionality of the kernel. If
#'   missing, this defaults to the dimensionality of \code{x}.
#' @param axis For \code{sobelFilter}, the axis along which to apply the
#'   operator, or 0 to apply it along all directions and generate a magnitude
#'   image. See also \code{\link{sobelKernel}}.
#' @return A morphed array with the same dimensions as the original array.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{morph}} for the function underlying these operations,
#'   and \code{\link{kernels}} for kernel-generating functions.
#' @rdname filters
#' @export
meanFilter <- function (x, kernel)
{
    return (morph(x, kernel, operator="i", merge="mean"))
}

#' @rdname filters
#' @export
medianFilter <- function (x, kernel)
{
    return (morph(x, kernel, operator="i", merge="median"))
}

#' @rdname filters
#' @export
sobelFilter <- function (x, dim, axis = 0)
{
    x <- as.array(x)
    if (missing(dim))
        dim <- length(base::dim(x))
    
    if (axis < 0 || axis > dim)
        stop("The axis should be between 0 and the dimensionality of the kernel")
    
    if (axis == 0)
    {
        squaredImages <- lapply(1:dim, function(i) sobelFilter(x,dim,i)^2)
        return (sqrt(Reduce("+", squaredImages)))
    }
    else
    {
        morphFun <- function(y,k) morph(y, k, operator="*", merge="sum")
        
        kernels <- lapply(1:dim, function(i) {
            if (i == axis)
                k <- sobelKernel(1)
            else
                k <- sobelKernel(1, 0)
            array(k, dim=replace(rep(1L,dim), i, 3L))
        })
        
        return (Reduce(morphFun, kernels, x))
    }
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
#' becomes the minimum value of \code{x - kernel}. Dilation is the opposite
#' operation, taking the maximum within the kernel.
#' 
#' @param x An object that can be coerced to an array, or for which a
#'   \code{\link{morph}} method exists.
#' @param kernel An array representing the kernel to be used. See
#'   \code{\link{shapeKernel}} for functions to generate a suitable kernel.
#' @return A morphed array with the same dimensions as the original array.
#' 
#' @examples
#' x <- c(0,0,1,0,0,0,1,1,1,0,0)
#' k <- c(1,1,1)
#' erode(x,k)
#' dilate(x,k)
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{morph}} for the function underlying all of these
#'   operations, \code{\link{kernels}} for kernel-generating functions,
#'   \code{\link{binarise}} for binarising an array, and
#'   \code{\link{gaussianSmooth}} for smoothing. The \code{EBImage}
#'   Bioconductor package also supplies functions to perform these operations,
#'   and may be slightly faster, but only works in two dimensions.
#' @rdname morphology
#' @aliases morphology
#' @export
erode <- function (x, kernel)
{
    x <- as.array(x)
    if (!isKernelArray(kernel))
        kernel <- kernelArray(kernel)
    
    greyscaleImage <- !binary(x)
    nNeighboursNot <- valueNot <- NULL
    
    if (greyscaleImage)
        operator <- ifelse(binary(kernel), "i", "-")
    else
    {
        # Kernels with zero at the origin mess up the heuristics
        operator <- "i"
        if (kernel[ceiling(length(kernel)/2)] != 0)
        {
            valueNot <- 0
            if (all(dim(kernel) <= 3))
                nNeighboursNot <- 3^length(dim(x)) - 1
        }
    }
    
    return (morph(x, kernel, operator=operator, merge="min", valueNot=valueNot, nNeighboursNot=nNeighboursNot))
}

#' @rdname morphology
#' @export
dilate <- function (x, kernel)
{
    x <- as.array(x)
    if (!symmetric(kernel))
        kernel <- array(rev(kernel), dim=dim(as.array(kernel)))
    if (!isKernelArray(kernel))
        kernel <- kernelArray(kernel)
    
    greyscaleImage <- !binary(x)
    nNeighboursNot <- value <- NULL
    
    if (greyscaleImage)
        operator <- ifelse(binary(kernel), "i", "+")
    else
    {
        operator <- "i"
        if (kernel[ceiling(length(kernel)/2)] != 0)
        {
            value <- 0
            if (all(dim(kernel) <= 3))
                nNeighboursNot <- 0
        }
    }
    
    return (morph(x, kernel, operator=operator, merge="max", value=value, nNeighboursNot=nNeighboursNot))
}

#' @rdname morphology
#' @export
opening <- function (x, kernel)
{
    x <- as.array(x)
    if (!isKernelArray(kernel))
        kernel <- kernelArray(kernel)
    
    return (dilate(erode(x, kernel), kernel))
}

#' @rdname morphology
#' @export
closing <- function (x, kernel)
{
    x <- as.array(x)
    if (!isKernelArray(kernel))
        kernel <- kernelArray(kernel)
    
    return (erode(dilate(x, kernel), kernel))
}

#' Skeletonise a numeric array
#' 
#' Skeletonisation is the process of thinning a shape to a medial line or
#' surface, and can be achieved using elementary mathematical morphology
#' operations in a number of ways. Three methods are available through this
#' function. They are all iterative and therefore relatively time-consuming.
#' 
#' The default method is Lantuéjoul's formula, a union across repeated
#' erosions, which works in any number of dimensions and may produce
#' reasonable results on greyscale images, but does not in general produce a
#' connected skeleton. Beucher introduced an alternative which may produce a
#' better result (although again the skeleton may not be connected), but this
#' implementation of the latter algorithm only applies to binary arrays. The
#' final method uses the so-called hit-or-miss transform, which searches for
#' exact patterns in the source array. This is guaranteed to produce a
#' connected skeleton, which is often desirable, but uses fixed kernels (so the
#' \code{kernel} argument is ignored) and is currently only implemented for 2D
#' binary arrays.
#' 
#' @param x An object that can be coerced to an array, or for which a
#'   \code{\link{morph}} method exists.
#' @param kernel An array representing the kernel to be used for the underlying
#'   morphology operations. The kernel is fixed for the \code{"hitormiss"}
#'   method, so this argument will be ignored.
#' @param method A string giving the method to use. See Details.
#' @return A skeletonised array with the same dimensions as the original array.
#' 
#' @examples
#' x <- c(0,0,1,0,0,0,1,1,1,0,0)
#' k <- c(1,1,1)
#' skeletonise(x,k)
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{morphology}}
#' @references
#' C. Lantuéjoul (1977). Sur le modèle de Johnson-Mehl généralisé. Technical
#' report, Centre de Morphologie Mathématique, Fontainebleau, France.
#' 
#' S. Beucher (1994). Digital skeletons in Euclidean and geodesic spaces.
#' Signal Processing 38(1):127-141.
#' \url{https://doi.org/10.1016/0165-1684(94)90061-2}.
#' @aliases skeletonize
#' @export skeletonise skeletonize
skeletonise <- skeletonize <- function (x, kernel = NULL, method = c("lantuejoul","beucher","hitormiss"))
{
    method <- match.arg(method)
    isBinary <- binary(x)
    x <- result <- as.array(x)
    nDims <- sum(dim(x) > 1)
    
    if (method == "lantuejoul")
    {
        result <- x - opening(x, kernel)
        eroded <- x
    
        repeat
        {
            eroded <- erode(eroded, kernel)
            if (all(eroded == 0))
                break
            result <- pmax(result, eroded - opening(eroded,kernel))
        }
    }
    else if (method == "beucher")
    {
        if (!isBinary)
            warning("The Beucher skeletonisation method will not preserve a greyscale")
        
        repeat
        {
            result <- x
            x <- erode(x,kernel) | (dilate(x - opening(x,kernel), kernel) & x)
            if (isTRUE(all.equal(x, result)))
                break
        }
    }
    else
    {
        if (!isBinary)
            stop("The hit-or-miss transform skeletonisation method only works with binary images")
        if (nDims != 2)
            stop("The hit-or-miss transform skeletonisation method is only implemented in 2D")
        
        x <- x / attr(isBinary, "value")
        k1 <- matrix(c(0,NA,1,0,1,1,0,NA,1), 3, 3)
        k2 <- matrix(c(NA,1,NA,0,1,1,0,0,NA), 3, 3)
        rotateKernel <- function(k) t(apply(k, 2, rev))
        hitOrMiss <- function(x,k) morph(x, k, operator="==", merge="all", value=1)
    
        repeat
        {
            result <- x
            for (i in 1:4)
            {
                x <- x & !hitOrMiss(x, k1)
                x <- x & !hitOrMiss(x, k2)
                k1 <- rotateKernel(k1)
                k2 <- rotateKernel(k2)
            }
        
            if (isTRUE(all.equal(x, result)))
                break
        }
    }
    
    result <- as.double(result)
    if (length(dim(x)) > 1)
        dim(result) <- dim(x)
    
    return (result)
}
