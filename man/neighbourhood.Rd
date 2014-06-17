\name{neighbourhood}
\alias{neighbourhood}

\title{Get neighbourhood information for an array}
\description{
This function provides information about the structure of a neighbourhood of a given width within a specified array.
}

\usage{
neighbourhood(x, width)
}

\arguments{
  \item{x}{An object that can be coerced to an array.}
  \item{width}{The width of neighbourhood to be considered. Only the first element will be used.}
}

\value{
A list with the following elements.
  \item{widths}{The width of the neighbourhood along each dimension. Currently all elements of this vector will be the same.}
  \item{size}{The number of pixels within the neighbourhood.}
  \item{locs}{A matrix giving the coordinates of each neighbourhood pixel, relative to the centre pixel, one per row.}
  \item{offsets}{Vector offsets of the neighbourhood values within \code{x}. These will only make sense if \code{x} is not less than \code{width} wide in any dimension.}
}

\author{Jon Clayden <code@clayden.org>}
