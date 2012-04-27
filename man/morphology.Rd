\name{morphology}
\alias{morphology}
\alias{erode}
\alias{dilate}
\alias{opening}
\alias{closing}

\title{Standard mathematical morphology operations}
\description{
These functions provide standard mathematical morphology operations, which can be applied to array data with any number of dimensions.
}

\usage{
erode(x, kernel)
dilate(x, kernel)
opening(x, kernel)
closing(x, kernel)
}

\arguments{
  \item{x}{An object that can be coerced to an array, or for which a \code{\link{morph}} method exists.}
  \item{kernel}{An array representing the kernel to be used. See \code{\link{shapeKernel}} for functions to generate a suitable kernel.}
}

\details{
The \code{erode} function uses the kernel as an eraser, centring it on each zero-valued voxel, which has the effect of eroding the extent of nonzero areas. Dilation has the opposite effect, extending the nonzero regions in the array. Opening is an erosion followed by a dilation, and closing is a dilation followed by an erosion, using the same kernel in both cases.
}

\value{
A morphed array with the same dimensions as the original array.
}

\author{Jon Clayden <jon.clayden+mmand@gmail.com>}

\seealso{\code{\link{morph}} for the function underlying all of these operations, \code{\link{kernels}} for kernel-generating functions, \code{\link{binarise}} for binarising an array, and \code{\link{gaussianSmooth}} for smoothing. The \code{EBImage} Bioconductor package also supplies functions to perform these operations, and may be slightly faster, but only works in two dimensions.}

\examples{
x <- c(0,0,1,0,0,0,1,1,1,0,0)
k <- c(1,1,1)
erode(x,k)
dilate(x,k)
}
