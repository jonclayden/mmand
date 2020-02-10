#' Distance transforms
#' 
#' The Euclidean distance transform produces an array like its argument, but
#' with element values representing the Euclidean distance to the nearest
#' nonzero element. The input is treated as logically binary, with all nonzero
#' values treated as "on", and all zeroes as "off".
#' 
#' @param x Any object. For the default method, this must be coercible to an
#'   array.
#' @param \dots Additional arguments to methods.
#' @return An array of the same dimension as the original, whose elements give
#'   the Euclidean distance from that element to the nearest "on" element in
#'   the original.
#' 
#' @examples
#' x <- c(0,0,1,0,0,0,1,1,1,0,0)
#' distanceTransform(x)
#' @author Jon Clayden <code@@clayden.org>
#' @export
distanceTransform <- function (x, ...)
{
    UseMethod("distanceTransform")
}

#' @rdname distanceTransform
#' @export
distanceTransform.default <- function (x, ...)
{
    x <- as.array(x)
    if (!is.numeric(x) && !is.logical(x))
        stop("Array must be numeric")
    
    returnValue <- .Call(C_distance_transform, x)
    
    if (length(dim(x)) > 1)
        dim(returnValue) <- dim(x)
    
    return (returnValue)
}
