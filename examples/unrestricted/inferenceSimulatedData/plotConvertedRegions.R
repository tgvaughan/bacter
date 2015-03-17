require(ggplot2)
require(reshape2)
require(stringr)

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

getSiteConversionCount <- function(convMapStr, seqLen) {
    res <- rep(0, seqLen)

    if (!is.na(convMapStr)) {
        for (regionStr in strsplit(convMapStr, ',')[[1]]) {
            bounds <- strtoi(strsplit(regionStr, ':')[[1]])
            res[(bounds[1]+1):(bounds[2]+1)] <- res[(bounds[1]+1):(bounds[2]+1)] + 1
        }
    }

    return(res)
}



getTruth <- function(filename) {
    start <- NULL
    end <- NULL
    visible <- NULL
    idx <- 1
    for (line in readLines(filename)) {
        if (!str_detect(line, "conversion"))
            next

        start[idx] <- strtoi(strsplit(strsplit(line, "site1=")[[1]][2], " ")[[1]][1])
        end[idx] <- strtoi(strsplit(strsplit(line, "site2=")[[1]][2], ";")[[1]][1])

        node1 <- strtoi(strsplit(strsplit(line, "node1=")[[1]][2], " ")[[1]][1])
        node2 <- strtoi(strsplit(strsplit(line, "node2=")[[1]][2], " ")[[1]][1])
        if (node1==node2)
            visible[idx] <- FALSE
        else
            visible[idx] <- TRUE
        
        idx <- idx + 1
    }

    indices <- factor(1:length(visible))
    
    return(data.frame(start=start, end=end, visible=visible, index=indices))
}

getTrueSiteConversionCount <- function(df, seqLen, includeInvisible=TRUE) {
    res <- rep(0, seqLen)
    for (i in 1:length(df$index)) {
        if (!includeInvisible && !df$visible[i])
            next;

        range <- (df$start[i]+1):(df$end[i]+1)
        res[range] <- res[range] + 1
    }

    return(res)
}

getSiteConversionProb <- function(df, seqLen) {
    res <- rep(0, seqLen)
    for (i in 1:length(df$acg.converted)) {
        res <-  res + getSiteConversionStatus(df$acg.converted[[i]], seqLen)
    }

    return (res/length(df$acg.converted))
}

plotConversionProb <- function(filename, seqLen, truthFile=NA, burnin=0.1) {
    df <- read.table(filename, as.is=T, header=T)
    N <- dim(df)[1]
    df <- df[-(1:(round(burnin*N))),]
    
    siteProbs <- getSiteConversionProb(df, seqLen)

    probdf <- data.frame(site=1:seqLen, prob=siteProbs)

    p <- ggplot(probdf, aes(x=site, y=prob)) + geom_line()
    p <- p + xlab("Site") + ylab("Posterior probability")
    p <- p + ggtitle("Per-site probability of conversion")

    if (!is.na(truthFile)) {
        truth <- getTruth(truthFile)
        p <- p + geom_rect(data=truth, mapping=aes(x=NULL, y=NULL, xmin=start, xmax=end, ymin=0, ymax=1,
                                           fill=index), alpha=0.5)
        p <- p + scale_fill_manual(values=rep(c("blue","red","green","purple","yellow", "pink"),
                                       length.out=length(truth$index)), guide=FALSE)
    }
    
    return(p)
}


getSiteMeanConversionCount <- function(df, seqLen) {
    N <- length(df$acg.converted)

    mean <- rep(0, seqLen)
    std <- rep(0, seqLen)
    for (i in 1:N) {
        count <- getSiteConversionCount(df$acg.converted[[i]], seqLen)
        mean <-  mean + count
        std <- std + count*count
    }

    mean <- mean/N
    std <- sqrt(std/N - mean*mean)

    res <- list()
    res$mean <- mean
    res$std <- std

    return (res)
}

plotMeanConversionCount <- function(filename, seqLen, truthFile=NA, burnin=0.1) {
    df <- read.table(filename, as.is=T, header=T)
    N <- dim(df)[1]
    df <- df[-(1:(round(burnin*N))),]
    
    siteCounts <- getSiteMeanConversionCount(df, seqLen)

    countdf <- data.frame(site=1:seqLen, mean=siteCounts$mean, std=siteCounts$std)

    #return (countdf)

    p <- ggplot(countdf, aes(x=site, y=mean, ymin=mean-std, ymax=mean+std))
    p <- p + geom_ribbon(show_guide=T, alpha=0.5)
    p <- p + geom_line()
    p <- p + scale_y_discrete()
    p <- p + xlab("Site") + ylab("Count")
    p <- p + ggtitle("Per-site mean conversion count")

    if (!is.na(truthFile)) {
        truth <- getTruth(truthFile)
        siteCountsTrue <- getTrueSiteConversionCount(truth, seqLen)
        trueCountDF <-  data.frame(site=1:seqLen, count=siteCountsTrue)
        p <- p + geom_line(data=trueCountDF,
                           mapping=aes(x=site, y=count, ymin=0, ymax=0), colour="red")
    }
    
    return(p)
}

plotConversionCounts <- function(filename, seqLen, truthFile=NA, burnin=0.1, includeInvisible=TRUE, ...) {
    df <- read.table(filename, as.is=T, header=T)
    N <- dim(df)[1]
    df <- df[-(1:(round(burnin*N))),]

    for (i in seq(1,dim(df)[1], by=ceiling(dim(df)[1]/1000))) {
        count <- getSiteConversionCount(df$acg.converted[[i]], seqLen)

        if (i==1) {
            plot(1:seqLen, count, 'l', col=rgb(1,0,0,0.01),
                 ylim=c(0, max(count)),
                 xlab='Site',
                 ylab='Count',
                 main="Conversion counts per site",
                 ...)
        } else {
            lines(1:seqLen, count, 's', col=rgb(0,0,0,0.01))
        }
    }

    if (!is.na(truthFile)) {
        truth <- getTruth(truthFile)
        siteCountsTrue <- getTrueSiteConversionCount(truth, seqLen, includeInvisible)
        trueCountDF <-  data.frame(site=1:seqLen, count=siteCountsTrue)
        lines(1:seqLen, siteCountsTrue, col=rgb(0.5,0,0), lwd=3)

        legend('topleft', inset=0.05, c('Truth'), col=rgb(0.5,0,0), lty=1, lwd=3)
    }
}

plotConvertedSiteTrace <- function(filename, seqLen, maxStates=400, maxSites=200) {
    df <- read.table(filename, as.is=T, header=T)

    Nstates <- length(df$acg.converted)
    NstatesDS <- min(maxStates, Nstates)
    seqLenDS <- min(maxSites, seqLen)

    sites <- seq(1, seqLen, by=floor(seqLen/seqLenDS))
    states <- seq(1, Nstates, by=floor(Nstates/NstatesDS))
    
    # Construct downsampled trace matrix
    m <- matrix(ncol=length(sites), nrow=length(states))
    for (i in 1:length(states)) {
        state <- states[i]
        row <-  getSiteConversionStatus(df$acg.converted[[state]], seqLen)
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
    p <- p + scale_fill_manual(name="Converted",
                               values=c("orange","brown"),
                               guide=guide_legend(reverse=TRUE))
    p <- p + xlab("MCMC steps") + ylab("Site")
    p <- p + ggtitle("Converted region trace")
    
    return(p)
}

plotConversionCountTrace <- function(filename, seqLen, maxStates=400, maxSites=200) {
    df <- read.table(filename, as.is=T, header=T)

    Nstates <- length(df$acg.converted)
    NstatesDS <- min(maxStates, Nstates)
    seqLenDS <- min(maxSites, seqLen)

    sites <- seq(1, seqLen, by=floor(seqLen/seqLenDS))
    states <- seq(1, Nstates, by=floor(Nstates/NstatesDS))
    
    # Construct downsampled trace matrix
    m <- matrix(ncol=length(sites), nrow=length(states))
    for (i in 1:length(states)) {
        state <- states[i]
        row <-  getSiteConversionCount(df$acg.converted[[state]], seqLen)
        for (j in 1:length(sites)) {
            site <- sites[j]
            m[i,j] <- row[site]
        }
    }

    # Convert to data frame
    d <- melt(m)
    d$count <- d$value
    d$states <- states[d$Var1]
    d$sites <- sites[d$Var2]

    # Plot
    p <- ggplot(d, aes(x=states, y=sites, fill=count)) + geom_raster()
    #p <- p + scale_fill_gradientn(colours=rainbow(20))
    p <- p + xlab("MCMC steps") + ylab("Site")
    p <- p + ggtitle("Conversion count trace")
    
    return(p)
}
