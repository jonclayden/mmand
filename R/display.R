display <- function (x, ...)
{
    UseMethod("display")
}

display.default <- function (x, transpose = TRUE, add = FALSE, ...)
{
    x <- as.matrix(x)
    if (transpose)
        x <- t(x)
    
    col <- grey(0:255/255)
    
    if (add)
        image(x[1:nrow(x),ncol(x):1], col=col, add=TRUE, ...)
    else
    {
        dev.new(width=nrow(x)/72, height=ncol(x)/72, dpi=72)
        oldPars <- par(mai=c(0,0,0,0), bg="grey70")
        image(x[1:nrow(x),ncol(x):1], col=col, asp=ncol(x)/nrow(x), add=FALSE, ...)
        par(oldPars)
    }
    
    invisible(NULL)
}
