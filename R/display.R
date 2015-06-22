.Cache <- new.env()

.checkDpi <- function ()
{
    if (!exists("dpiDevice",.Cache) || !identical(.Cache$dpiDevice,getOption("device")))
    {
        .Cache$dpiDevice <- getOption("device")
        deviceFunction <- try(match.fun(.Cache$dpiDevice), silent=TRUE)
        if ((class(deviceFunction) == "function") && ("dpi" %in% names(formals(deviceFunction))))
            .Cache$dpi <- structure(c(72,72), explicit=TRUE)
        else
        {
            dev.new(width=1, height=1)
            .Cache$dpi <- structure(dev.size(units="px"), explicit=FALSE)
            dev.off()
        }
    }
    
    return (.Cache$dpi)
}

#' Display a 2D image
#' 
#' This function displays a 2D greyscale image. It is a wrapper around
#' \code{image}, with more sensible defaults for images. It is (S3) generic.
#' 
#' Relative to the defaults for \code{image} (from the \code{graphics}
#' package), this function transposes and then inverts the matrix along the
#' y-direction, uses a grey colour scale, fills the entire device with the
#' image, and tries to size the image correctly given the dot pitch of the
#' display. Unfortunately the latter is not always possible, due to downstream
#' limitations.
#' 
#' @param x An object that can be coerced to a numeric matrix.
#' @param transpose Whether to transpose the matrix before display. This is
#'   usually necessary due to the conventions of \code{image}.
#' @param useRaster Whether to use raster graphics if possible. This is
#'   generally preferred for speed. Passed to \code{image}.
#' @param add Whether to add the image to an existing plot.
#' @param col The colour scale to use. The default is 256 grey levels.
#' @param \dots Additional arguments to \code{image}.
#' @return This function is called for its side-effect of displaying an image
#'   on a new R device.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @export
display <- function (x, ...)
{
    UseMethod("display")
}

#' @rdname display
#' @export
display.default <- function (x, transpose = TRUE, useRaster = TRUE, add = FALSE, col = grey(0:255/255), ...)
{
    x <- as.matrix(x)
    if (transpose)
        x <- t(x)
    
    dpi <- .checkDpi()
    
    if (add)
        image(x[1:nrow(x),ncol(x):1], col=col, useRaster=useRaster, add=TRUE, ...)
    else
    {
        if (attr(dpi,"explicit") && dpi[1] == dpi[2])
            dev.new(width=nrow(x)/dpi[1], height=ncol(x)/dpi[2], dpi=dpi[1])
        else
            dev.new(width=nrow(x)/dpi[1], height=ncol(x)/dpi[2])
        
        oldPars <- par(mai=c(0,0,0,0), bg="grey70")
        image(x[1:nrow(x),ncol(x):1], col=col, asp=ncol(x)/nrow(x), useRaster=useRaster, add=FALSE, ...)
        par(oldPars)
    }
    
    invisible(NULL)
}
