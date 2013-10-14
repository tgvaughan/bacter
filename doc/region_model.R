# Simulate generation of regions along sequence of length L
sim <- function (rho, delta, L) {

    svec <- rep(0,L)
    svec[1] <- 1

    nrec <- 0
    tractLengths <- c()
    
    for (k in 1:(L-1)) {

        if (svec[k]==1) {
            if (runif(1)<rho/L) {
                svec[k+1] <- 2
                nrec <- nrec + 1
                tractLength <- 1
            } else
                svec[k+1] = 1
        } else {
            if (runif(1)<1/delta) {
                svec[k+1] <- 1
                tractLengths <- append(tractLengths, tractLength)
                
            } else {
                svec[k+1] <- 2
                tractLength <- tractLength + 1
            }
        }
    }

    res <- list()
    res$tractLengths <- tractLengths
    res$nrec <- nrec
    res$svec <- svec

    return(res)
}


# Plot dependence of expected number of regions against rho and delta
# (rho is in units of tree length)

expectedCount <- function(rho, delta, L) {
    return ((rho/2)/(rho*delta/(2*L) + 1))
}


pdf('regionCount.pdf', width=7, height=5)

rho <- 10^seq(-3,3.5,by=0.1)
delta <- 1e3
L <- 1.6e6

plot(rho, expectedCount(rho, delta, L), 'l', lwd=2,
     xlab=expression(paste(rho,"'",sep="")),
     ylab='Recombination count',
     ylim=c(0,1500))
lines(rho, rho/2, lty=2, lwd=2);

legend('topleft', inset=.05, c('No overlaps', 'With overlaps'),
       lty=c(1,2), )


dev.off()
