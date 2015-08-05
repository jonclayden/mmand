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
#' This function displays a 2D greyscale or RGB colour image. It is a wrapper
#' around \code{image}, with more sensible defaults for images. It is (S3)
#' generic. A method for 3D arrays is provided, which assumes that the third
#' dimension corresponds to channel (grey/alpha for two channels, red/green/
#' blue for three, red/green/blue/alpha for four).
#' 
#' Relative to the defaults for \code{image} (from the \code{graphics}
#' package), this function transposes and then inverts the matrix along the
#' y-direction, uses a grey colour scale, fills the entire device with the
#' image, and tries to size the image correctly given the dot pitch of the
#' display. Unfortunately the latter is not always possible, due to downstream
#' limitations.
#' 
#' @param x An R object. For the default method, it must be coercible to a
#'   numeric matrix.
#' @param transpose Whether to transpose the matrix before display. This is
#'   usually necessary due to the conventions of \code{image}.
#' @param useRaster Whether to use raster graphics if possible. This is
#'   generally preferred for speed. Passed to \code{image}.
#' @param add Whether to add the image to an existing plot. If \code{TRUE},
#'   zero values in the image will be converted to \code{NA}s for plotting
#'   purposes, to make them transparent. This will not affect the original
#'   image data.
#' @param col The colour scale to use. The default is 256 grey levels. The
#'   array method overrides this appropriately.
#' @param max The maximum colour value for each channel. If the array has
#'   integer mode, this is fixed to 255. Passed to \code{\link{rgb}}.
#' @param \dots Additional arguments to \code{image}, or the default method.
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
    {
        x[x==0] <- NA
        image(x[1:nrow(x),ncol(x):1], col=col, useRaster=useRaster, add=TRUE, ...)
    }
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

#' @rdname display
#' @export
display.matrix <- function (x, ...)
{
    display.default(x, ...)
}

#' @rdname display
#' @export
display.array <- function (x, max = 1, ...)
{
    if (length(dim(x)) != 3)
        stop("Only three-dimensional arrays may be displayed")
    
    mode <- storage.mode(x)
    max <- ifelse(mode == "integer", 255L, as.double(max))
    
    dim3 <- dim(x)[3]
    if (dim3 == 1L)
        display.default(drop(x), ...)
    else
    {
        x[x < 0] <- as(0, mode)
        x[x > max] <- max
        
        if (dim3 == 2L)
            cols <- rgb(x[,,1], x[,,1], x[,,1], x[,,2], maxColorValue=max)
        else if (dim3 == 3L)
            cols <- rgb(x[,,1], x[,,2], x[,,3], maxColorValue=max)
        else if (dim3 == 4L)
            cols <- rgb(x[,,1], x[,,2], x[,,3], x[,,4], maxColorValue=max)
        else
            stop("Third dimension should not be greater than 4")
        
        uniqueCols <- unique(cols)
        indices <- match(cols, uniqueCols)
        dim(indices) <- dim(x)[1:2]
        display.default(indices, col=uniqueCols, ...)
    }
}
