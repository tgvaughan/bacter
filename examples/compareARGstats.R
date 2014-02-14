# Compare two sets of ARG stats

removeBurnin <- function(df, burninFrac=0.1) {

    n <- length(df)
    burnin <- ceiling(burninFrac*n)

    return(df[-(1:burnin),])
}

compareStats <- function(df1, df2, label1="df1", label2="df2", noCF=FALSE) {

    par(mfrow=c(2,3))

    qqplot(df1$arg.nRecomb, df2$arg.nRecomb,
           xlab=label1, ylab=label2, col='blue',
           main='Recombination count')
    lines(c(0,1e10), c(0,1e10))

    qqplot(df1$arg.meanTractLength, df2$arg.meanTractLength,
           xlab=label1, ylab=label2, col='blue',
           main='Tract length')
    lines(c(0,1e10), c(0,1e10))

    qqplot(df1$arg.meanInterTractLength, df2$arg.meanInterTractLength,
           xlab=label1, ylab=label2, col='blue',
           main='Inter-tract length')
    lines(c(0,1e10), c(0,1e10))

    qqplot(df1$arg.meanDepartureHeight, df2$arg.meanDepartureHeight,
           xlab=label1, ylab=label2, col='blue',
           main='Departure height')
    lines(c(0,1e10), c(0,1e10))

    qqplot(df1$arg.meanEdgeLength, df2$arg.meanEdgeLength,
           xlab=label1, ylab=label2, col='blue',
           main='Edge length')
    lines(c(0,1e10), c(0,1e10))

    if (!noCF) {
        qqplot(df1$arg.CFheight, df2$arg.CFheight,
               xlab=label1, ylab=label2, col='blue',
               main='Clonal frame height')
        lines(c(0,1e10), c(0,1e10))
    }
}

