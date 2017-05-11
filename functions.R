simulateT2 <- function(msCommand, path2ms = "./"){
  # Simulate T2 values (the coalescence time of two genes) using ms
  #
  # Args:
  #   msCommand: string, the corresponding ms command (must have the -T flag)
  #   path2ms: string, the route to the folder where the ms software is
  #
  # Returns:
  #   a vector of coalescence times (T2)
  msOutput <- system(paste(path2ms, msCommand, sep = ""), intern = TRUE)
  reTrees <- "1:[0-9]*\\.[0-9]*,2:[0-9]*\\.[0-9]*"
  trees <- msOutput[grepl(reTrees, msOutput)] # Get only the lines with a tree in newick format
  reT2value <- "[0-9]*\\.[0-9]*"
  T2values <- regmatches(trees, regexpr(reT2value, trees))  # Get the T2 values
  return(as.numeric(T2values))
}

computeEmpiricalIICR <- function(T2values, tVector){
  # Compute the IICR based on independent values of T2 and a given distretization
  # of time
  #
  # Args:
  #   T2values: numeric, a vector of independent values of T2
  #   tVector: numeric, the discrete time windows for computing the IICR
  #
  # Returns:
  #   a list of two vectors: the tVector and the corresponding IICR values
  
  T2values <- simulateT2(msCommand)
  tVector <- seq(0, 3, 0.01)
  if(max(T2values) > max(tVector)) tVector <- c(tVector, max(T2values))
  counts <- hist(T2values, breaks = tVector)$counts
  fx <- counts / (tVector[-1] - tVector[-length(tVector)])
  Fx <- cumsum(counts)
  IICR <- (length(T2values) - Fx) / fx
  return(list(tVector = tVector[-length(tVector)], IICR = IICR))
}

msCommand <- "ms 2 100000 -T"
msCommandNislands <- "ms 2 100000 -I 10 2 0 0 0 0 0 0 0 0 0 -T"

tVector <- seq(0, 5, 0.01)
testT2 <- simulateT2(msCommandNislands)
test <- computeEmpiricalIICR(testT2, tVector)
plot(test$tVector, test$IICR, type = 's')


