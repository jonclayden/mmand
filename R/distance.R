#' @export
distanceTransform <- function (x)
{
    x <- as.array(x)
    if (!is.numeric(x) && !is.logical(x))
        stop("Array must be numeric")
    
    returnValue <- .Call(C_distance_transform, x)
    
    if (length(dim(x)) > 1)
        dim(returnValue) <- dim(x)
    
    return (returnValue)
}
