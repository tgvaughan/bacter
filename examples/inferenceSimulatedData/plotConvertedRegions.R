require(ggplot2)
require(reshape2)

getSiteConversionStatus <- function(convMapStr, seqLen) {
    res <- rep(0, seqLen)

    if (!is.na(convMapStr)) {
        for (regionStr in strsplit(convMapStr, ',')[[1]]) {
            bounds <- strtoi(strsplit(regionStr, ':')[[1]])
            res[(bounds[1]+1):(bounds[2]+1)] <- 1
        }
    }

    return(res)
}

getSiteConversionProb <- function(df, seqLen) {
    res <- rep(0, seqLen)
    for (i in 1:length(df$arg.converted)) {
        res <-  res + getSiteConversionStatus(df$arg.converted[[i]], seqLen)
    }

    return (res/length(df$arg.converted))
}

plotConversionProb <- function(filename, seqLen) {
    df <- read.table(filename, as.is=T, header=T)
    siteProbs <- getSiteConversionProb(df, seqLen)

    probdf <- data.frame(site=1:seqLen, prob=siteProbs)
    
    p <- ggplot(probdf, aes(x=site, y=prob)) + geom_line()
    p <- p + xlab("Site") + ylab("Posterior probability")
    p <- p + ggtitle("Per-site probability of conversion")

    return(p)
}


plotConvertedSiteTrace <- function(filename, seqLen, maxStates=400, maxSites=200) {
    df <- read.table(filename, as.is=T, header=T)

    Nstates <- length(df$arg.converted)
    NstatesDS <- min(maxStates, Nstates)
    seqLenDS <- min(maxSites, seqLen)

    sites <- seq(1, seqLen, by=floor(seqLen/seqLenDS))
    states <- seq(1, Nstates, by=floor(Nstates/NstatesDS))
    
    # Construct downsampled trace matrix
    m <- matrix(ncol=length(sites), nrow=length(states))
    for (i in 1:length(states)) {
        state <- states[i]
        row <-  getSiteConversionStatus(df$arg.converted[[state]], seqLen)
        for (j in 1:length(sites)) {
            site <- sites[j]
            m[i,j] <- row[site]
        }
    }

    # Convert to data frame
    d <- melt(m)
    d$converted <- c(F,T)[d$value + 1]
    d$states <- states[d$Var1]
    d$sites <- sites[d$Var2]

    # Plot
    p <- ggplot(d, aes(x=states, y=sites, fill=converted)) + geom_raster()
    p <- p + scale_fill_manual(name="Converted", values=c("orange","brown"))
    p <- p + xlab("MCMC steps") + ylab("Alignment sites")
    p <- p + ggtitle("Converted region trace")
    
    return(p)
}
