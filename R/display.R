display <- function (x, ...)
{
    UseMethod("display")
}

display.default <- function (x, transpose = TRUE, useRaster = TRUE, add = FALSE, ...)
{
    x <- as.matrix(x)
    if (transpose)
        x <- t(x)
    
    col <- grey(0:255/255)
    dpi <- 72
    
    if (add)
        image(x[1:nrow(x),ncol(x):1], col=col, useRaster=useRaster, add=TRUE, ...)
    else
    {
        dev.new(width=nrow(x)/dpi, height=ncol(x)/dpi, dpi=dpi)
        oldPars <- par(mai=c(0,0,0,0), bg="grey70")
        image(x[1:nrow(x),ncol(x):1], col=col, asp=ncol(x)/nrow(x), useRaster=useRaster, add=FALSE, ...)
        par(oldPars)
    }
    
    invisible(NULL)
}
