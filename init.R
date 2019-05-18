suppressPackageStartupMessages({
  if (!exists("PBM")) {
    source("pbm/load.R", chdir = TRUE)
    load("data/data.rda")
  }
  source("utils.R")
  sourceCpp("decision-funcs.cpp")
})

## common values used throughout
ixDM <- as.integer(as.factor(data$Sid))
idDM <- levels(as.factor(data$Sid))
nrDM <- length(unique(ixDM))
ixQ <- as.integer(data$Qid)
nrQ <- length(unique(ixQ))
Nr <- length(ixDM)

with(data, {
  dQ$zero <- 0
  dD$zero <- 0
  dQ$aveExpectation <- (dQ$meanA + dQ$meanB)/2
  dQ$x_at_max_p_A <- local({
    X <- as.matrix(dQ[xaNames])
    P <- as.matrix(dQ[paNames])
    sapply(1:nrQ, function(i){
      X[i, which.max(P[i, ])]
    })
  })
  dQ$x_at_max_p_B <- local({
    X <- as.matrix(dQ[xbNames])
    P <- as.matrix(dQ[pbNames])
    sapply(1:nrQ, function(i){
      X[i, which.max(P[i, ])]
    })
  })
})

