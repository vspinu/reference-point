#!/usr/bin/Rscript

## In order to replicate models in the appendix run this script as:
## sec 2 model 1               : ./build.R 17 50000 T T T T pow prelec1 (same as ./build.R)
## sec 3 model 1 with exp util : ./build.R 17 50000 T T T T exp_log prelec1
## sec 4 model 1 with prelec2  : ./build.R 17 50000 T T T T pow prelec2
## sec 5 model 1 with ibeta    : ./build.R 17 50000 T T T T pow ibeta
## sec 6 model 2               : ./build.R 17 50000 F T T T pow prelec1
## sec 7 model 3               : ./build.R 17 50000 T T T F pow prelec1
## sec 8 model 4               : ./build.R 17 50000 F T T F pow prelec1

## ./build.R 17 50000 T T T T exp_log prelec1
## parallel ./build.R 17 50000 T T T T pow ibeta,prelec2
## parallel ./build.R 17 50000 T,F T T T,F pow prelec1

source("utils.R")
assignCommandArgs(seed = 17, N = 50000,
                  I_nonLSKR = T, I_LSKR = T, W_nonLSKR = T, W_LSKR = T,
                  util_type = "pow", weight_type = "prelec1")

set.seed(seed)
LS <- KR <- TRUE

model <- encode_model(seed, N,
                      I_nonLSKR = I_nonLSKR, I_LSKR = I_LSKR,
                      W_nonLSKR = W_nonLSKR, W_LSKR = W_LSKR,
                      LS = LS, KR = KR,
                      fU = util_type, fW = weight_type)
model <- paste0(model, Sys.getenv("SUFFIX"))

cat("\n------------- MODEL ", model, "-----------\n\n")

source("init.R")
source("model-RP-rule-estim.R")

modelfile <- sprintf("%s/model.rds", model)

if (file.exists(modelfile)) {
  cat("Model file exists; skipping simulation.")
} else {

  tdata <- buildData.rule_mixture_DM(data,
                                     I_nonLSKR = I_nonLSKR, I_LSKR = I_LSKR,
                                     W_nonLSKR = W_nonLSKR, W_LSKR = W_LSKR,
                                     LS = LS, KR = KR)

  M <- model_rule_mixture(N, tdata, util_type, weight_type)

  dir.create(model, FALSE, TRUE)
  cat("Saving", modelfile, "...\n")
  saveRDS(M, file = modelfile)
  result_file <- sprintf("results/%s.rds", model)
  saveRDS(M, file = result_file)
}

cat("\n\n---------- CREATING REPORT FOR MODEL ", model, "-----------\n\n")

library(knitr)
Sys.setenv("MODEL-DIR" = model)
setwd(model)
knitr::knit('../report.Rnw', output="report.tex")

cat("\n---------- COMPILING PDF ------------------------------------\n\n")
system("pdflatex -interaction=nonstopmode report.tex")
