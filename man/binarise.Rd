\name{binarise}
\alias{binarise}

\title{Binarise a numeric array}
\description{
This function binarises an array, setting all nonzero elements to unity.
}

\usage{
binarise(x)
}

\arguments{
  \item{x}{An object that can be coerced to an array, or for which a \code{\link{morph}} method exists.}
}

\value{
A morphed array with the same dimensions as the original array.
}

\author{Jon Clayden <code@clayden.org>}

\seealso{\code{\link{morph}} for the function underlying this operation, and \code{\link{erode}} for mathematical morphology functions.}
