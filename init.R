suppressPackageStartupMessages({
  if (!exists("PBM")) {
    source("pbm/load.R", chdir = TRUE)
  }
  source("utils.R")
  Rcpp::sourceCpp("decision-funcs.cpp")
})

