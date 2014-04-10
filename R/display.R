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

display <- function (x, ...)
{
    UseMethod("display")
}

display.default <- function (x, transpose = TRUE, useRaster = TRUE, add = FALSE, ...)
{
    x <- as.matrix(x)
    if (transpose)
        x <- t(x)
    
    dpi <- .checkDpi()
    col <- grey(0:255/255)
    
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
