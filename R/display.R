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

.checkAttribs <- function (x, ...)
{
    result <- list(...)
    attribs <- attributes(x)
    for (name in names(attribs))
        result[[name]] <- attribs[[name]]
    
    return (result)
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
#' If \code{x} has attributes \code{"range"}, \code{"background"}, \code{"asp"}
#' or \code{"dpi"}, these are respected.
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
#' @param max The maximum colour value for each channel. If \code{NULL}, the
#'   default, this is taken from the \code{"range"} attribute, if there is one,
#'   otherwise it is 255 for integer-mode arrays, and 1 otherwise. Passed to
#'   \code{\link{rgb}}.
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
    dpiDevice <- .checkDpi()
    attribs <- .checkAttribs(x, background="grey70", range=range(x[is.finite(x)]), asp=ncol(x)/nrow(x), dpi=dpiDevice)
    if (is.null(attr(x,"asp")) && !is.null(attr(x,"dpi")))
        attribs$asp <- (ncol(x) * attribs$dpi[1]) / (nrow(x) * attribs$dpi[2])
    
    x <- as.matrix(x)
    if (transpose)
        x <- t(x)
    
    if (add)
    {
        x[x==0] <- NA
        image(x[1:nrow(x),ncol(x):1], col=col, useRaster=useRaster, zlim=sort(attribs$range), add=TRUE, ...)
    }
    else
    {
        if (attr(dpiDevice,"explicit") && attribs$dpi[1] == attribs$dpi[2])
            dev.new(width=nrow(x)/attribs$dpi[1], height=ncol(x)/attribs$dpi[2], dpi=attribs$dpi[1])
        else
            dev.new(width=nrow(x)/attribs$dpi[1], height=ncol(x)/attribs$dpi[2])
        
        oldPars <- par(mai=c(0,0,0,0), bg=attribs$background)
        image(x[1:nrow(x),ncol(x):1], col=col, asp=attribs$asp, useRaster=useRaster, zlim=sort(attribs$range), add=FALSE, ...)
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
display.array <- function (x, max = NULL, ...)
{
    if (length(dim(x)) == 2)
        return (display.default(x, ...))
    else if (length(dim(x)) != 3)
        stop("Only three-dimensional arrays may be displayed")
    
    mode <- storage.mode(x)
    if (is.null(max))
    {
        if (is.null(attr(x, "range")))
            max <- ifelse(mode == "integer", 255L, 1)
        else
            max <- max(attr(x, "range"))
    }
    
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
        attr(indices, "range") <- c(1L, length(uniqueCols))
        for (attrib in c("background","asp","dpi"))
            attr(indices, attrib) <- attr(x, attrib)
        display.default(indices, col=uniqueCols, ...)
    }
}

#' Show an ASCII art representation of a 2D image or matrix
#' 
#' This function prints a rough, text-only representation of an image argument
#' to the R terminal, mapping image intensities to a 10-level pseudo-greyscale.
#' The image is first rescaled to fit into the terminal or other specified
#' width, and downsampled in the row direction to correct for nonsquare
#' character shapes.
#' 
#' The result is a compact representation of a matrix that can be used for
#' visualising kernel arrays, sparse matrices and other non-images.
#' 
#' @note If the terminal does not used a fixed-width font, the result is
#'   unlikely to be useful.
#' 
#' @param x An object that can be coerced to a numeric matrix or array. 3D
#'   arrays with third dimension no greater than 4 will be taken as
#'   multichannel 2D images, and their channels averaged before display.
#'   Plain vectors and 1D arrays will be treated as single-row matrices.
#' @param invert By default the mapping uses heavier type for brighter areas.
#'   If this option is \code{TRUE}, the sense of the scale will be reversed.
#' @param width The width of sketch to draw, in characters.
#' @param squash The factor by which to scale the row direction of the image.
#'   Generally this should be markedly less than one, to preserve the aspect
#'   ratio of the image, since most fixed-width font characters are taller than
#'   they are wide.
#' @return This function is called for the side-effect of printing an ASCII
#'   representation of its argument.
#' 
#' @examples
#' sketch(shapeKernel(c(9,15), type="diamond"))
#' sketch(shapeKernel(c(9,15), type="diamond"), squash=1)
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{display}}
#' @export
sketch <- function (x, invert = FALSE, width = getOption("width"), squash = 0.5)
{
    # Ref: http://paulbourke.net/dataformats/asciiart/
    scale <- unlist(strsplit(" .:-=+*#%@", "", fixed=TRUE))
    if (invert)
        scale <- rev(scale)
    n <- length(scale)
    
    attribs <- .checkAttribs(x, range=range(x[is.finite(x)]))
    
    x <- drop(as.array(x))
    if (length(dim(x)) == 1)
        dim(x) <- c(1, length(x))
    else if (length(dim(x)) == 3 && dim(x)[3] < 5)
        x <- apply(x, 1:2, mean)
    else if (length(dim(x)) > 3)
        stop("Only 2D single-channel or multichannel images may be sketched")
    
    scaleFactors <- c(squash, 1)
    if (ncol(x) > width)
        scaleFactors <- scaleFactors * width / ncol(x)
    x <- rescale(x, scaleFactors, triangleKernel())
    
    indices <- round((x - min(attribs$range)) * (n-1) / (max(attribs$range) - min(attribs$range))) + 1
    chars <- structure(scale[indices], dim=dim(x))
    
    apply(chars, 1, function(line) cat(paste0(paste(line,collapse=""),"\n")))
    
    invisible(NULL)
}
