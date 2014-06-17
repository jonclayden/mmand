\name{sampleKernelFunction}
\alias{sampleKernelFunction}
\alias{plot.kernelFunction}
\alias{plot.kernelArray}

\title{Sampling and plotting kernels}
\description{
These functions can be used to sample and plot kernel profiles.
}

\usage{
sampleKernelFunction(kernel, values)
\method{plot}{kernelFunction}(x, y, xlim = c(-2,2), lwd = 2, col = "red", ...)
\method{plot}{kernelArray}(x, y, axis = 1, lwd = 2, col = "red", ...)
}

\arguments{
  \item{kernel}{A kernel function object.}
  \item{values}{A vector of values to sample the function at. These are in units of pixels, with zero representing the centre of the kernel.}
  \item{x}{A kernel object of the appropriate class.}
  \item{y}{Ignored.}
  \item{xlim}{The limits of the range used to profile the kernel.}
  \item{lwd}{The line width to use for the kernel profile.}
  \item{col}{The line colour to use for the kernel profile.}
  \item{axis}{The axis to profile along.}
  \item{\dots}{Additional plot parameters.}
}

\value{
For \code{sampleKernelFunction} a vector of kernel values at the locations requested. The \code{plot} methods are called for their side-effects.
}

\author{Jon Clayden <code@clayden.org>}

\seealso{\code{\link{kernels}} for kernel-generating functions.}

\examples{
sampleKernelFunction(mnKernel(), -2:2)
plot(mnKernel())
}
