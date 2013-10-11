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
