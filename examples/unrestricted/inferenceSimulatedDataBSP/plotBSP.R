df <- read.table('inferenceSimulatedData.log', header=T)

getPopSizes <- function(t, df, idstr="popModel") {

    Nindices <- getNindices(df, idstr)
    tindices <- getTindices(df, idstr)

    res <- list()
    res$mean <- rep(0, length(t))
    res$upper <- rep(0, length(t))
    res$lower <- rep(0, length(t))

    for (tidx in 1:length(t)) {
        thist <- t[tidx]

        N <- apply(df, 1, function (dfrow) {
                   rowsub <- (thist-dfrow[tindices])>0
                   Nindex <- 1
                   if (sum(rowsub)>0)
                       Nindex <- min((1:length(Nindices))[rowsub]) + 1

                   return (dfrow[Nindices[Nindex]])
        })

        q <- quantile(N, probs=c(0.025, 0.5, 0.975))
        res$lower[tidx] <- q[1]
        res$median[tidx] <- q[2]
        res$upper[tidx] <- q[3]
    }

    return(res)
}


getNindices <- function(df, idstr) {
    pattern <- paste("^", idstr, ".N", sep='')

    return(which(regexpr(pattern, names(df))>0))
}

getTindices <- function(df, idstr) {
    pattern <- paste("^", idstr, ".t", sep='')

    return(which(regexpr(pattern, names(df))>0))
}

plotBSP <- function(t, df, idstr="popModel", ...) {

    N <- getPopSizes(t, df, idstr)

    plot(t, N$median, 'l', ylim=c(0.9*min(N$lower), 1.1*max(N$upper)), ...)
    polygon(c(t, rev(t)), c(N$lower, rev(N$upper)), col="grey", border=NA)
    lines(t, N$median, 'l', lwd=2)
}
