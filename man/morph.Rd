\name{morph}
\alias{morph}
\alias{morph.default}

\title{Morph an array with a kernel}
\description{
The \code{morph} function applies a kernel, literally or multiplicatively, to a target array. Optionally, applying the kernel to a particular array element can be made conditional on its value, or the number of nonzero immediate neighbours that it has. The \code{morph} function is (S3) generic.
}

\usage{
morph(x, ...)

\method{morph}{default}(x, kernel, brush = TRUE, eraser = FALSE, value = NULL,
      valueNot = NULL, nNeighbours = NULL, nNeighboursNot = NULL)
}

\arguments{
  \item{x}{Any object. For the default method, this must be coercible to an array.}
  \item{kernel}{An object representing the kernel to be applied, which must be coercible to an array. It must have odd width in all dimensions, but does not have to be isotropic in size. The kernel's dimensionality may be less than that of the target array, \code{x}. See \code{\link{kernels}} for kernel-generating functions.}
  \item{brush}{If \code{TRUE}, the kernel is placed into the morphed array literally, centred at each point that meets the specified criteria (unless it is an eraser: see below). If \code{FALSE}, the value of the final morphed array at each point is the sum of the surrounding values, using the kernel as a set of weights.}
  \item{eraser}{If \code{TRUE}, the kernel is used as an eraser, so that all nonzero values in the kernel are used to indicate that the corresponding array elements should be set to zero in the morphed array. Ignored if \code{brush} is \code{FALSE}.}
  \item{value}{An optional vector of values in the target array for which to apply the kernel. Takes priority over \code{valueNot} if both are specified.}
  \item{valueNot}{An optional vector of values in the target array for which not to apply the kernel.}
  \item{nNeighbours}{An optional numeric vector giving allowable numbers of nonzero neighbours (including diagonal neighbours) for array elements where the kernel will be applied. Takes priority over \code{nNeighboursNot} if both are specified.}
  \item{nNeighboursNot}{An optional numeric vector giving nonallowable numbers of nonzero neighbours (including diagonal neighbours) for array elements where the kernel will be applied.}
}

\value{
A morphed array with the same dimensions as the original array.
}

\author{Jon Clayden <jon.clayden+mmand@gmail.com>}

\seealso{\code{\link{kernels}} for kernel-generating functions, and \code{\link{erode}} for more specific mathematical morphology functions. \code{\link{gameOfLife}} shows how this function can be used for non-morphological purposes, in that case to power a cellular automaton. See also the \code{kernel} and \code{kernapply} functions in the \code{stats} package, particularly if you want to smooth time series.}
