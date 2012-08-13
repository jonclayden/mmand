\name{gameOfLife}
\alias{gameOfLife}
\alias{gosperGliderGun}

\title{Conway's Game of Life}
\description{
An implementation of Conway's Game of Life, a classical cellular automaton, using the \code{\link{morph}} function. The \code{\link{gosperGliderGun}} function provides an interesting starting configuration.
}

\usage{
gameOfLife(init, size, density = 0.3, steps = 200, viz = FALSE, tick = 0.5)

gosperGliderGun()
}

\arguments{
  \item{init}{The initial state of the automaton, a binary matrix. If missing, the initial state will be randomly generated, with a population density given by \code{density}.}
  \item{size}{The dimensions of the board. Defaults to the size of \code{init}, but must be given if that parameter is missing. If both are specified and \code{size} is larger than the dimensions of \code{init}, then the latter will be padded with zeroes.}
  \item{density}{The approximate population density of the starting state. Ignored if \code{init} is provided.}
  \item{steps}{The number of generations of the automaton to simulate.}
  \item{viz}{If \code{TRUE}, the state of the system at each generation is plotted.}
  \item{tick}{The amount of time, in seconds, to pause before plotting each successive generation. Ignored if \code{viz} is \code{FALSE}.}
}

\details{
Conway's Game of Life is a simple cellular automaton, based on a 2D matrix of ``cells''. It shows complex behaviour based on four simple rules. These are:
\enumerate{
  \item Any live cell with fewer than two live neighbours dies, as if caused by under-population.
  \item Any live cell with two or three live neighbours lives on to the next generation.
  \item Any live cell with more than three live neighbours dies, as if by overcrowding.
  \item Any dead cell with exactly three live neighbours becomes a live cell, as if by reproduction.
}
Live and dead cells are represented by 1s and 0s in the matrix, respectively.

The initial state and the rules above completely determine the behaviour of the system. The Gosper glider gun is an interesting starting configuration that generates so-called ``gliders'', which propagate across the board.

In principle the size of the board in a cellular automaton is infinite. Of course this is not easy to simulate, but this implementation adds a border of two extra cells around the board on all sides to approximate an infinite board slightly better. These are not visualised, nor returned in the final state.
}

\value{
A binary matrix representing the final state of the system after \code{steps} generations.
}

\author{Jon Clayden <jon.clayden+mmand@gmail.com>}

\seealso{The \code{\link{morph}} function, which powers this simulation.}

\examples{
\dontrun{gameOfLife(init=gosperGliderGun(), size=c(40,40), steps=50, viz=TRUE)}
}
