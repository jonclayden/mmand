\name{resample}
\alias{resample}
\alias{resample.default}
\alias{rescale}

\title{Resample an array}
\description{
The \code{resample} function uses a kernel function to resample a target array. It is (S3) generic. The \code{rescale} is an alternative interface for the common case where the image is being scaled to a new size.
}

\usage{
resample(x, points, kernel, ...)

\method{resample}{default}(x, points, kernel, pointType = c("auto","general","grid"), ...)

rescale(x, factor, kernel, ...)
}

\arguments{
  \item{x}{Any object. For the default method, this must be coercible to an array.}
  \item{points}{Either a matrix giving the points to sample at, one per row, or a list giving the locations on each axis, which will be made into a grid.}
  \item{kernel}{A kernel function object, used to provide coefficients for each resampled value.}
  \item{pointType}{A string giving the type of the point specification being used. Usually can be left as \code{"auto"}.}
  \item{factor}{A vector of scale factors, which will be recycled to the dimensionality of \code{x}.}
}

\value{
If a generalised sampling scheme is used (i.e. with \code{points} a matrix), the result is a vector of sampled values. For a grid scheme (i.e. with \code{points} a list), it is a resampled array.
}

\author{Jon Clayden <code@clayden.org>}

\seealso{\code{\link{kernels}} for kernel-generating functions.}

\examples{
resample(c(0,0,1,0,0), seq(0.75,5.25,0.5), triangleKernel())
}
