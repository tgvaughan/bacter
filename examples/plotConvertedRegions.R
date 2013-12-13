# Take converted region sequence from log file and plot

lines <- readLines("addRemoveTest.converted")[-1]
nLines <- length(lines)

allConvStates <- c()

for (line in lines) {
    convState <- rep(0,10000)
    rangeStrings <- strsplit(strsplit(line,"\t")[[1]][2], ",")[[1]]
    for (rangeString in rangeStrings) {
        range <- strtoi(strsplit(rangeString, ":")[[1]])
        
        convState[(range[1]+1):(range[2]+1)] <- 1
    }

    allConvStates <- cbind(allConvStates,convState)
}
