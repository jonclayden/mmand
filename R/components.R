#' Check for a symmetric array
#' 
#' This function checks whether a numeric array is symmetric, in the sense of
#' transposition. This is tested by comparing the reversed vectorised array to
#' the unreversed equivalent.
#' 
#' @param x An object that can be coerced to a numeric array.
#' @return A logical value indicating whether the array is symmetric or not.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @export
symmetric <- function (x)
{
    x <- as.array(x)
    if (!is.numeric(x))
        report(OL$Error, "Array must be numeric")
    
    return (.Call("is_symmetric", x, PACKAGE="mmand"))
}

#' Find connected components
#' 
#' The \code{components} function finds connected components in a numeric
#' array. The kernel determines which neighbours are considered connected (e.g.
#' including or excluding diagonal neighbours), and will usually have width 3
#' in each dimension.
#' 
#' @param x Any object. For the default method, this must be coercible to an
#'   array.
#' @param kernel An object representing the kernel to be used, which must be
#'   coercible to an array. It must have odd width in all dimensions, but does
#'   not have to be isotropic in size. The kernel's dimensionality may be less
#'   than that of the target array, \code{x}. See \code{\link{kernels}} for
#'   kernel-generating functions.
#' @param \dots Additional arguments to methods.
#' @return An array of the same dimension as the original, whose integer-valued
#'   elements identify the component to which each element in the array
#'   belongs. Zero values in the original array will result in NAs.
#' 
#' @examples
#' x <- c(0,0,1,0,0,0,1,1,1,0,0)
#' k <- c(1,1,1)
#' components(x,k)
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{kernels}} for kernel-generating functions.
#' @export
components <- function (x, kernel, ...)
{
    UseMethod("components")
}

#' @rdname components
#' @export
components.default <- function (x, kernel, ...)
{
    x <- as.array(x)
    if (!is.numeric(x))
        report(OL$Error, "Target array must be numeric")
    
    if (!isKernelArray(kernel))
        kernel <- kernelArray(kernel)
    
    if (any(dim(kernel) %% 2 != 1))
        report(OL$Error, "Kernel must have odd width in all dimensions")
    
    if (!symmetric(kernel))
        report(OL$Error, "Kernel must be symmetric")
    
    if (length(dim(kernel)) < length(dim(x)))
        dim(kernel) <- c(dim(kernel), rep(1,length(dim(x))-length(dim(kernel))))
    else if (length(dim(kernel)) > length(dim(x)))
        report(OL$Error, "Kernel has greater dimensionality than the target array")
    
    storage.mode(x) <- "double"
    
    returnValue <- .Call("connected_components", x, kernel, PACKAGE="mmand") + 1
    
    if (length(dim(x)) > 1)
        dim(returnValue) <- dim(x)
    
    return (returnValue)
}
