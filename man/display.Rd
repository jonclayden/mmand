\name{display}
\alias{display}
\alias{display.default}

\title{Display a 2D image}
\description{
This function displays a 2D greyscale image. It is a wrapper around \code{image}, with more sensible defaults for images. It is (S3) generic.
}

\usage{
display(x, ...)

\method{display}{default}(x, transpose = TRUE, useRaster = TRUE, add = FALSE, ...)
}

\arguments{
  \item{x}{An object that can be coerced to a numeric matrix.}
  \item{transpose}{Whether to transpose the matrix before display. This is usually necessary due to the conventions of \code{image}.}
  \item{useRaster}{Whether to use raster graphics if possible. This is generally preferred for speed. Passed to \code{image}.}
  \item{add}{Whether to add the image to an existing plot.}
  \item{\dots}{Additional arguments to \code{image}.}
}

\details{
Relative to the defaults for \code{image} (from the \code{graphics} package), this function transposes and then inverts the matrix along the y-direction, uses a grey colour scale, fills the entire device with the image, and tries to size the image correctly given the dot pitch of the display. Unfortunately the latter is not always possible, due to downstream limitations.
}

\value{
This function is called for its side-effect of displaying an image on a new R device.
}

\author{Jon Clayden <code@clayden.org>}
