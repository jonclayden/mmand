\name{gaussianSmooth}
\alias{gaussianSmooth}

\title{Smooth a numeric array with a Gaussian kernel}
\description{
This function smoothes an array using a Gaussian kernel with a specified standard deviation.
}

\usage{
gaussianSmooth(x, sigma)
}

\arguments{
  \item{x}{An object that can be coerced to an array, or for which a \code{\link{morph}} method exists.}
  \item{sigma}{A numeric vector giving the standard deviation of the kernel in each dimension. Can have lower dimensionality than the target array.}
}

\value{
A morphed array with the same dimensions as the original array.
}

\author{Jon Clayden <code@clayden.org>}

\seealso{\code{\link{morph}} for the function underlying this operation, \code{\link{gaussianKernel}} for generating Gaussian kernels (which is also used by this function), and \code{\link{erode}} for mathematical morphology functions.}
