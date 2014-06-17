\name{filters}
\alias{filters}
\alias{medianFilter}
\alias{meanFilter}

\title{Apply a filter to an array}
\description{
These functions apply mean or median filters to an array.
}

\usage{
medianFilter(x, kernel)
meanFilter(x, kernel)
}

\arguments{
  \item{x}{An object that can be coerced to an array, or for which a \code{\link{morph}} method exists.}
  \item{kernel}{A kernel array, indicating the scope of the filter.}
}

\value{
A morphed array with the same dimensions as the original array.
}

\author{Jon Clayden <code@clayden.org>}

\seealso{\code{\link{morph}} for the function underlying this operation, and \code{\link{kernels}} for kernel-generating functions.}
