.Workspace <- new.env()

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
    
    if (!exists("dpi",.Workspace) || is.null(.Workspace$dpi))
    {
        dev.new(width=1, height=1)
        .Workspace$dpi <- dev.size(units="px")
        dev.off()
    }
    
    if (add)
        image(x[1:nrow(x),ncol(x):1], col=col, useRaster=useRaster, add=TRUE, ...)
    else
    {
        dev.new(width=nrow(x)/.Workspace$dpi[1], height=ncol(x)/.Workspace$dpi[2])
        oldPars <- par(mai=c(0,0,0,0), bg="grey70")
        image(x[1:nrow(x),ncol(x):1], col=col, asp=ncol(x)/nrow(x), useRaster=useRaster, add=FALSE, ...)
        par(oldPars)
    }
    
    invisible(NULL)
}
