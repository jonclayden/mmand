#' Conway's Game of Life
#' 
#' An implementation of Conway's Game of Life, a classical cellular automaton,
#' using the \code{\link{morph}} function. The \code{\link{gosperGliderGun}}
#' function provides an interesting starting configuration.
#' 
#' Conway's Game of Life is a simple cellular automaton, based on a 2D matrix
#' of ``cells''. It shows complex behaviour based on four simple rules. These
#' are: \enumerate{ \item Any live cell with fewer than two live neighbours
#' dies, as if caused by under-population.  \item Any live cell with two or
#' three live neighbours lives on to the next generation.  \item Any live cell
#' with more than three live neighbours dies, as if by overcrowding.  \item Any
#' dead cell with exactly three live neighbours becomes a live cell, as if by
#' reproduction. } Live and dead cells are represented by 1s and 0s in the
#' matrix, respectively.
#' 
#' The initial state and the rules above completely determine the behaviour of
#' the system. The Gosper glider gun is an interesting starting configuration
#' that generates so-called ``gliders'', which propagate across the board.
#' 
#' In principle the size of the board in a cellular automaton is infinite. Of
#' course this is not easy to simulate, but this implementation adds a border
#' of two extra cells around the board on all sides to approximate an infinite
#' board slightly better. These are not visualised, nor returned in the final
#' state.
#' 
#' @aliases gameOfLife gosperGliderGun
#' @param init The initial state of the automaton, a binary matrix. If missing,
#' the initial state will be randomly generated, with a population density
#' given by \code{density}.
#' @param size The dimensions of the board. Defaults to the size of
#' \code{init}, but must be given if that parameter is missing. If both are
#' specified and \code{size} is larger than the dimensions of \code{init}, then
#' the latter will be padded with zeroes.
#' @param density The approximate population density of the starting state.
#' Ignored if \code{init} is provided.
#' @param steps The number of generations of the automaton to simulate.
#' @param viz If \code{TRUE}, the state of the system at each generation is
#' plotted.
#' @param tick The amount of time, in seconds, to pause before plotting each
#' successive generation. Ignored if \code{viz} is \code{FALSE}.
#' @return A binary matrix representing the final state of the system after
#' \code{steps} generations.
#' @author Jon Clayden <code@@clayden.org>
#' @seealso The \code{\link{morph}} function, which powers this simulation.
#' @examples
#' 
#' \dontrun{gameOfLife(init=gosperGliderGun(), size=c(40,40), steps=50, viz=TRUE)}
#' 
gameOfLife <- function (init, size, density = 0.3, steps = 200, viz = FALSE, tick = 0.5)
{
    state <- NULL
    
    if (missing(init) && missing(size))
        report(OL$Error, "At least one of \"init\" and \"size\" must be specified")
    else if (missing(init))
    {
        state <- ifelse(runif(prod(size)) < density, 1L, 0L)
        dim(state) <- size
    }
    else if (missing(size))
        state <- init
    else
    {
        if (length(size) != 2)
            report(OL$Error, "The size of the initial matrix should have length 2")
        if (!is.matrix(init))
            report(OL$Error, "Initial state should be specified as a matrix")
        state <- matrix(0L, nrow=size[1], ncol=size[2])
        state[1:nrow(init),1:ncol(init)] <- init
    }
    
    if (!is.matrix(state) || length(dim(state)) != 2)
        report(OL$Error, "Initial state is not a 2D matrix")
    
    stateWithBorder <- matrix(0L, nrow=nrow(state)+4, ncol=ncol(state)+4)
    stateWithBorder[(1:nrow(state))+2,(1:ncol(state))+2] <- state
    
    if (viz)
    {
        oldPars <- par(mai=c(0,0,0,0))
        image(state, asp=ncol(state)/nrow(state))
    }
    
    for (i in seq_len(steps))
    {
        # Rule 2 is a survival rule, so nothing changes
        rule1Diff <- morph(stateWithBorder, kernel=1L, operator="0", value=1, nNeighbours=0:1) - stateWithBorder
        rule3Diff <- morph(stateWithBorder, kernel=1L, operator="0", value=1, nNeighbours=4:8) - stateWithBorder
        rule4Diff <- morph(stateWithBorder, kernel=1L, operator="1", value=0, nNeighbours=3L) - stateWithBorder
        
        prevState <- state
        stateWithBorder <- stateWithBorder + rule1Diff + rule3Diff + rule4Diff
        state <- stateWithBorder[(1:nrow(state))+2,(1:ncol(state))+2]
        
        if (isTRUE(all.equal(prevState, state)))
        {
            report(OL$Info, "State is stable after ", i-1, " steps")
            break
        }
        
        if (viz)
        {
            Sys.sleep(tick)
            image(state, add=TRUE)
        }
    }
    
    if (viz)
        par(oldPars)
    
    invisible(state)
}

gosperGliderGun <- function ()
{
    state <- matrix(0L, nrow=11, ncol=38)
    state[c(17,18,28,29,127,128,129,137,141,147,153,158,164,172,181,185,193,194,195,205,235,236,237,246,247,248,256,260,277,278,282,283,389,390,400,401)] <- 1L
    invisible(state)
}
