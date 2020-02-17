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
#' @param pixdim An optional numeric vector or logical value. In the former
#'   case it will be taken as giving the physical size of the array elements of
#'   \code{x} along each dimension, and these will be incorporated into the
#'   distance calculation. If \code{TRUE}, the default, the \code{"pixdim"}
#'   attribute of \code{x} will be used for this purpose, if it is present. If
#'   \code{FALSE}, any such attribute will be ignored, and distances will
#'   always be counted in array elements, with all dimensions treated equally.
#' @param signed Logical value. If \code{TRUE}, the signed distance transform
#'   is returned, such that distances from the region boundary are negative
#'   within the region and positive outside. Otherwise, distances are zero
#'   within the region.
#' @return An array of the same dimension as the original, whose elements give
#'   the Euclidean distance from that element to the nearest "on" element in
#'   the original.
#' 
#' @examples
#' x <- c(0,0,1,0,0,0,1,1,1,0,0)
#' distanceTransform(x)
#' distanceTransform(x, pixdim=2)
#' @author Jon Clayden <code@@clayden.org>
#' @export
distanceTransform <- function (x, ...)
{
    UseMethod("distanceTransform")
}

#' @rdname distanceTransform
#' @export
distanceTransform.default <- function (x, pixdim = TRUE, signed = FALSE, ...)
{
    x <- as.array(x)
    if (!is.numeric(x) && !is.logical(x))
        stop("Array must be numeric")
    
    if (is.numeric(pixdim))
    {
        if (length(pixdim) == length(dim(x)))
        {
            attr(x, "pixdim") <- pixdim
            pixdim <- TRUE
        }
        else
        {
            warning("Specified pixdim vector is of the wrong length - ignoring it")
            pixdim <- FALSE
        }
    }
    
    isBinary <- binary(x)
    if (signed)
    {
        if (!isBinary)
        {
            x <- binarise(x)
            value <- 1
        }
        else
            value <- attr(isBinary, "value")
        returnValue <- .Call(C_distance_transform, x, pixdim) - .Call(C_distance_transform, value-x, pixdim)
    }
    else
        returnValue <- .Call(C_distance_transform, x, pixdim)
    
    if (length(dim(x)) > 1)
        dim(returnValue) <- dim(x)
    
    return (returnValue)
}
