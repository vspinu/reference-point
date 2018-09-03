suppressPackageStartupMessages({
    if (!exists("PBM")) {
        library(devtools)
        library(Rcpp)
        load_all("~/dev/protoClasses/")
        sourceCpp("~/dev/pbm/src/sample.cpp")
        source("~/dev/pbm/R/hierarchy.R")
        load("../data/data.rda")
    }
    source("./utils.R")
    sourceCpp("../decision-funcs.cpp")
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

setup_for_fold <- function(cv_type = 10){
    if (cv_type == 10){
        SEED <<- 17
        FOLDS <<- split(1:70, sample(rep(1:10, times = 7)))
    } else if (cv_type == 70) {
        SEED <<- 70
        FOLDS <<- split(1:70, 1:70)
    } else {
        stop ("invalid cv_type:", cv_type)
    }
    set.seed(SEED)
    CVDIR <<- sprintf("./data/CV%s/", cv_type)
    dir.create(CVDIR, F)
    folds_file <- sprintf("%sfolds%d.rds", CVDIR, SEED)
    if (!file.exists(folds_file))
        saveRDS(FOLDS, folds_file)
}

