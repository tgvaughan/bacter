getPopSizes <- function(t, df, idstr) {

    Nindices <- getNindices(df, idstr)
    tindices <- getTindices(df, idstr)

    res <- list()
    res$mean <- rep(0, length(t))
    res$upper <- rep(0, length(t))
    res$lower <- rep(0, length(t))

    for (tidx in 1:length(t)) {
        thist <- t[tidx]

        N <- apply(df, 1, function (dfrow) {
                   Nindex <- max(1,findInterval(thist, dfrow[tindices]))
                   return (dfrow[Nindices[Nindex]])
        })

        q <- quantile(N, probs=c(0.025, 0.5, 0.975))
        res$lower[tidx] <- q[1]
        res$median[tidx] <- q[2]
        res$upper[tidx] <- q[3]
    }

    return(res)
}

getPopSizesLinear <- function(t, df, idstr) {

    Nindices <- getNindices(df, idstr)
    tindices <- getTindices(df, idstr)

    res <- list()
    res$mean <- rep(0, length(t))
    res$upper <- rep(0, length(t))
    res$lower <- rep(0, length(t))

    for (tidx in 1:length(t)) {
        thist <- t[tidx]

        N <- apply(df, 1, function (dfrow) {
                   Nindex <- findInterval(thist, dfrow[tindices])

                   if (Nindex==0)
                       return(dfrow[Nindices[1]])

                   if (Nindex==length(Nindices))
                       return(dfrow[Nindices[Nindex]])

                   N0 <- dfrow[Nindices[Nindex]]
                   N1 <- dfrow[Nindices[Nindex+1]]

                   t0 <- dfrow[tindices[Nindex]]
                   t1 <- dfrow[tindices[Nindex+1]]

                   return(N0 + (N1-N0)/(t1-t0)*(thist-t0))
        })

        q <- quantile(N, probs=c(0.025, 0.5, 0.975))
        res$lower[tidx] <- q[1]
        res$median[tidx] <- q[2]
        res$upper[tidx] <- q[3]
    }

    return(res)
}

getNindices <- function(df, idstr) {
    pattern <- paste("^", gsub(":", ".", idstr), ".N", sep='')

    return(which(regexpr(pattern, names(df))>0))
}

getTindices <- function(df, idstr) {
    pattern <- paste("^", gsub(":", ".", idstr), ".t", sep='')

    return(which(regexpr(pattern, names(df))>0))
}

plotBSP <- function(t, df, burnin=0.1, idstr="popModel", linear=FALSE, ...) {

    frameLen <- dim(df)[1]
    if (!linear)
        N <- getPopSizes(t, df[ceiling(burnin*frameLen):frameLen,], idstr)
    else
        N <- getPopSizesLinear(t, df[ceiling(burnin*frameLen):frameLen,], idstr)

    ellipsis <- list(...)

    if (length(ellipsis$xlab) == 0)
        ellipsis$xlab = "t"

    if (length(ellipsis$ylab) == 0)
        ellipsis$ylab = "N"

    if (length(ellipsis$ylim) == 0)
        ellipsis$ylim = c(0.9*min(N$lower), 1.1*max(N$upper))

    args <- c(list(t, N$median, 'l'), ellipsis)

    do.call(plot, args)
    polygon(c(t, rev(t)), c(N$lower, rev(N$upper)), col="grey", border=NA)
    lines(t, N$median, 'l', lwd=2)
}

redraw <- function(...) {
    #t <- 10^seq(-4,-2, length.out=100)
    t <- seq(0,0.01, length.out=100)
    df <- read.table('inferenceSimulatedDataLinear.log', header=T)

    plotBSP(t, df, ...)
    lines(t, exp(-1000*t), lty=2, lwd=2)
}
