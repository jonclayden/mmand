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
        rule1Diff <- morph(stateWithBorder, kernel=0L, brush=TRUE, value=1, nNeighbours=0:1) - stateWithBorder
        rule3Diff <- morph(stateWithBorder, kernel=0L, brush=TRUE, value=1, nNeighbours=4:8) - stateWithBorder
        rule4Diff <- morph(stateWithBorder, kernel=1L, brush=TRUE, value=0, nNeighbours=3L) - stateWithBorder
        
        prevState <- state
        stateWithBorder <- stateWithBorder + rule1Diff + rule3Diff + rule4Diff
        state <- stateWithBorder[(1:nrow(state))+2,(1:ncol(state))+2]
        
        if (equivalent(prevState, state))
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
