\name{kernels}
\alias{kernels}
\alias{isKernel}
\alias{kernelArray}
\alias{isKernelArray}
\alias{shapeKernel}
\alias{gaussianKernel}
\alias{kernelFunction}
\alias{isKernelFunction}
\alias{boxKernel}
\alias{triangleKernel}
\alias{mitchellNetravaliKernel}
\alias{mnKernel}

\title{Functions for generating and manipulating kernels}
\description{
These functions can be used to generate kernels for morphological, smoothing or resampling operations. There are two types of kernels: kernel arrays, which are used with \code{\link{morph}}, and kernel functions, which are used with \code{\link{resample}}.
}

\usage{
isKernel(object)

kernelArray(values)
isKernelArray(object)

shapeKernel(width, dim = length(width), type = c("box","disc","diamond"),
    binary = TRUE, normalised = FALSE)
gaussianKernel(sigma, dim = length(sigma), size = 6*sigma, normalised = TRUE)

kernelFunction(name = c("box","triangle","mitchell-netravali"), ...)
isKernelFunction(object)

boxKernel()
triangleKernel()
mitchellNetravaliKernel(B = 1/3, C = 1/3)
mnKernel(B = 1/3, C = 1/3)
}

\arguments{
  \item{object}{Any object.}
  \item{values}{A numeric vector or array, containing the values of the kernel array.}
  \item{width}{A numeric vector giving the width of the shape in each dimension, in array elements. Does not need to be integer-valued, or equal for all dimensions. Will be recycled to length \code{dim} if that parameter is also specified.}
  \item{sigma}{A numeric vector giving the standard deviation of the underlying Gaussian distribution in each dimension, in array elements. Does not need to be equal for all dimensions. Will be recycled to length \code{dim} if that parameter is also specified.}
  \item{dim}{An integer value giving the dimensionality of the kernel. Defaults to the length of \code{width} or \code{sigma}.}
  \item{size}{A numeric vector giving the width of the kernel in each dimension, which will be rounded up to the nearest odd integer. Defaults to four times the corresponding \code{sigma} value.}
  \item{type}{A string giving the type of shape to produce. In one dimension, these shapes are all equivalent.}
  \item{binary}{If \code{FALSE}, the value of the kernel at each point represents the proportion of the array element within the shape. If \code{TRUE}, these values are binarised to be 1 if at least half of the element is within the shape, and 0 otherwise.}
  \item{normalised}{If \code{TRUE}, the sum of non-missing elements of the kernel will be unity. Note that this is the default for \code{gaussianKernel}, but not for \code{shapeKernel}.}
  \item{name}{A string giving the name of the kernel function required.}
  \item{\dots}{Parameters for the kernel function.}
  \item{B,C}{Mitchell-Netravali coefficients, each of which must be between 0 and 1.}
}

\details{
There are two forms of kernel used by this package. Kernel arrays, otherwise known in mathematical morphology as structuring elements, are numeric arrays with class \code{kernelArray}. They are defined on a grid of odd width, and are used by \code{\link{morph}} and related functions. Kernel functions, by contrast, are represented in R as a list containing a name and, optionally, some parameters. The real implementation is in C++. They are defined everywhere within the support of the kernel, and are used by \code{\link{resample}} and friends. The key distinction is in whether the kernel will always be centred exactly on the location of an existing value in the data (for kernel arrays) or not (for kernel functions).
}

\value{
An array of the requested dimensionality, which can be used as a kernel with \code{\link{morph}} and related functions.
}

\author{Jon Clayden <code@clayden.org>}

\seealso{\code{\link{morph}} for general application of kernels to data, \code{\link{erode}} for mathematical morphology functions, and \code{\link{gaussianSmooth}} for applying smoothing to data.}

\examples{
shapeKernel(c(3,5), type="diamond")
gaussianKernel(c(0.3,0.3))
mnKernel()
}
