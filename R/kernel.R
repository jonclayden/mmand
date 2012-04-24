shapeKernel <- function (width, dim = length(width), type = c("box","disc","diamond"), brush = TRUE, binary = TRUE, normalised = FALSE)
{
    type <- match.arg(type)
    
    if (dim > length(width))
        width <- rep(width, length.out=dim)
    
    # This corresponds to the smallest isotropic kernel of odd width to fit the shape
    size <- max(ifelse(width %% 2 == 1, width, width+1))
    kernel <- array(ifelse(brush,NA,0), dim=rep(size,dim))
    
    x <- 1:size - (size+1)/2
    x <- lapply(max(width)/width, "*", x)
    normFun <- switch(type, box=function(a,b) pmax(abs(a),abs(b)),
                            disc=function(a,b) sqrt(a^2 + b^2),
                            diamond=function(a,b) (abs(a) + abs(b)))
    norms <- Reduce(function(a,b) outer(a,b,FUN=normFun), x)
    
    if (binary)
        kernel <- ifelse(kernel==0, 0L, 1L)
    else if (normalised)
        kernel <- kernel / sum(kernel)
    
    return (kernel)
}

gaussianKernel <- function (sigma, dim = length(sigma), size = NULL)
{
    
}
