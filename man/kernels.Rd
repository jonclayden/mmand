\name{kernels}
\alias{kernels}
\alias{shapeKernel}
\alias{gaussianKernel}

\title{Kernel-generating functions}
\description{
These functions can be used for generating kernels for morphological or smoothing operations.
}

\usage{
shapeKernel(width, dim = length(width), type = c("box","disc","diamond"),
      brush = TRUE, binary = TRUE, normalised = FALSE)

gaussianKernel(sigma, dim = length(sigma), size = 4*sigma, normalised = TRUE)
}

\arguments{
  \item{width}{A numeric vector giving the width of the shape in each dimension, in array elements. Does not need to be integer-valued, or equal for all dimensions. Will be recycled to length \code{dim} if that parameter is also specified.}
  \item{sigma}{A numeric vector giving the standard deviation of the underlying Gaussian distribution in each dimension, in array elements. Does not need to be equal for all dimensions. Will be recycled to length \code{dim} if that parameter is also specified.}
  \item{dim}{An integer value giving the dimensionality of the kernel. Defaults to the length of \code{width} or \code{sigma}.}
  \item{size}{A numeric vector giving the width of the kernel in each dimension, which will be rounded up to the nearest odd integer. Defaults to four times the corresponding \code{sigma} value.}
  \item{type}{A string giving the type of shape to produced. In one dimension, these shapes are all equivalent.}
  \item{brush}{If \code{TRUE}, elements of the kernel which are outside the shape will be set to \code{NA}, otherwise they will have zero value. This is to ensure appropriate behaviour in \code{\link{morph}}.}
  \item{binary}{If \code{FALSE}, the value of the kernel at each point represents the proportion of the array element within the shape. If \code{TRUE}, these values are binarised to be 1 if at least half of the element is within the shape, and 0 otherwise.}
  \item{normalised}{If \code{TRUE}, the sum of non-missing elements of the kernel will be unity. Note that this is the default for \code{gaussianKernel}, but not for \code{shapeKernel}.}
}

\value{
An array of the requested dimensionality, which can be used as a kernel with \code{\link{morph}} and related functions.
}

\author{Jon Clayden <jon.clayden+mmand@gmail.com>}

\seealso{\code{\link{morph}} for general application of kernels to data, \code{\link{erode}} for mathematical morphology functions, and \code{\link{gaussianSmooth}} for applying smoothing to data.}

\examples{
shapeKernel(c(3,5), type="diamond")
gaussianKernel(c(0.3,0.3))
}
