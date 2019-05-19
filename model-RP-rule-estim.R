## Mixture of stochastic reference point models

## Each vector of RPs consists of a matrix rpX (one RP per row) and a
## probability distributions in over those in rpP mattrix of the same size as
## rpX. Now, this model takes a bunch of RPs and computes a mixture of them for
## each individual.

## - MODEL_DATA environemnt must contain xa, xb, rpX, rpP, PREFS and is usually
##   a result from the call to buildData.rule_mixture_DM. 

model_rule_mixture <- function(nr_sims = 5000, model_data,
                               util_type = "pow", weight_type = "prelec1",
                               exclude_rps = c(), nr_sims_warm = round(max(100, nr_sims *.1))){
 
  if(weight_type %in% c("prelec1", "pow1", "linlog1", "TK")){
    DM <- c(l = 1, u = 1, w = 1, xi = 10)
  } else if (weight_type %in% c("prelec2", "pow2", "linlog2", "ibeta")) {
    DM <- c(l = 1, u = 1, w1 = 1, w2 = 1, xi = 10)
  } else stop ("invalid weight_type argument ", weight_type)

  if(!util_type %in% c("pow", "exp_log")) stop("invalid util_type argument ", util_type)

  model_data$util_type <- util_type
  model_data$weight_type <- gsub("[0-9]", "", weight_type)

  M <- pbm("M", rootParentEnv = model_data)
  ## DM <- c(l = 1, u = 1, w = 1, xi = 10) 
  ## DM <- c(l = 1, u = 1, w1 = 1, w2 = 1, xi = 10) 
  ## DM <- c(l = 1, u = 1, xi = 10)
  nrRules <- model_data$nrRules

  ## control which rules to estimate for CV
  mnames <- model_data$mnames
  rp_hyper_priors <- setNames(rep.int(1, nrRules), mnames)
  rp_hyper_priors[exclude_rps] <- 1e-5
  rp_ixs <- which(!mnames %in% exclude_rps)
  stopif(length(rp_ixs) == 0)

  ## defBC and defP calls use standard R dynamic scope
  Nr <- get("Nr", envir = model_data)
  nrRP <- get("nrRP", envir = model_data)
  nrDM <- get("nrDM", envir = model_data)
  ixQ <- get("ixQ", envir = model_data)
  ixRP <- get("ixRP", envir = model_data)
  ixDM <- get("ixDM", envir = model_data)
  
  M$initCells(DM =
                defBC(prototype = "pd(LogNorm)",
                      setFields = list(
                        size = nrDM, scale=1, var = DM),
                      prDM =
                        defP("pd(LogGaNorm)",
                             cix_dim=2, size=length(DM), names = names(DM),
                             hpDM =
                               defP(prototype = "hc",
                                    cix_dim = 1, 
                                    setFields = list(
                                      size = length(DM),
                                      ## changing this prior has stronger impact
                                      ## on more complex models (prelec2 and
                                      ## ibeta) and virtually no impact on
                                      ## prelec1 models

                                      # 1. flat (most sampled values are Inf)
                                      var = c(mu0=-5, n0=.1, alpha0=.01, beta0=.01)
                                      # 2. flat with median 1, but all finite
                                      # values (pbm/distributions.R has the
                                      # code):
                                      ## $ median    : num 1.09
                                      ## $ mean      : num 1743
                                      ## $ sd        : num 481900
                                      ## c(mu0=-1, n0=10, alpha0=5, beta0=10)
                                    )))),

              ## Related RPXs (aka sampled from dirichlet population)
              RP =
                defBC("pd(Cat)",
                      N = nrRules, size = nrRP, varnames = "RP", 
                      st = sample(rp_ixs, nrRP, T), 
                      prRP = defP("pd(conj)(Dirich)",
                                  cix = 1, cix_dim = 1,
                                  varsize = nrRules, 
                                  hpRP = defP("hc", var = rp_hyper_priors)))

              ## ## Unrelated RPXs
              ## RP = defBC("pd(DiscrUnif)",
              ##     size = nrRP, varnames = "RP",
              ##     st = sample(1:nrRules, nrRP, T), 
              ##     hpRP = defP("hc", var = c(min = 1, max = nrRules), size = 1))

              ## Stochastic RP. Estimate the mixture density directly.
              ## RP = defBC("pd(Dirich)",
              ##     size = nrRP, scale = 1, varsize = nrRules, ##var = rep.int(1, nrRules),
              ##     st = dirichlet_rng(nrRP, rep.int(.1, nrRules)), 
              ##     hpRP = defP(prototype = "hc",
              ##         cix = 1:nrRules, size = nrRules, 
              ##         var = rep.int(.1, nrRules)))
              )

  M$initCells(DATA =
                defBC(prototype = "dc",
                      setFields = list(
                        st = as.matrix(ixQ),
                        foldable_objects = c("xa", "xb", "pa", "pb", "DM", "PREFS")), 
                      initForms = list(
                        set.DM = form({
                          DM <- e(PST("DM"))
                          .U <- DM[, "u", drop = F]
                          .L <- DM[, "l"]
                          {
                            if( "w" %in% colnames(DM))
                              .W <- DM[, "w", drop = F]
                            if( "w2" %in% colnames(DM))
                              .W <- DM[, c("w1", "w2")]
                          }
                        }), 
                        init.M.build.DM = form({
                          .W <- rep.int(1, Nr)
                          .W2 <- rep.int(1, Nr)
                          e(set.DM)
                        })
                        ## init.M.build.rpP = expression(.rpP <- matrix(rep.int(1, Nr)))),
                        ## init.M.build.rpP =
                        ##   form(.rpP <- cbind(rep.int(1, Nr), 0, 0, 0))
                      ),
                      setForms = list(
                        set.ll = expression({
                          e(set.DM)

                          ## decision model mixture
                          ## rpP <- PST("RP")
                          ## delta <- c_SPT_pow(xb, pb, RPXs, rpP, .U, .L) -
                          ##   c_SPT_pow(xa, pa, RPXs, rpP, .U, .L)

                          ## stochastic or rule mixture
                          ## .rpX <- matrix(RPXs[cbind(1:Nr, PST("RP"))])
                          .RP <- PST("RP")
                          .mixs <- 1:Nr + (.RP - 1)*Nr
                          .rpXa <- RPXa[.mixs, ]
                          .rpPa <- RPPa[.mixs, ]
                          .rpXb <- RPXb[.mixs, ]
                          .rpPb <- RPPb[.mixs, ]

                          .D <- ifelse(Is[.RP], meanDiff, 0)
                          .W[!Ws[.RP], ] <- 1 # ifelse(Ws[.RP], .W, 1)

                          ## meanDiff = meanA - meanB
                          delta <-
                            c_SPT(xb, pb, .rpXb, .rpPb, .L, .U, .W, util_type, weight_type) -
                            c_SPT(xa, pa, .rpXa, .rpPa, .L, .U, .W, util_type, weight_type) - .D

                          expdb <- exp(delta*DM[, "xi"])
                          prob  <- expdb/(1+expdb)
                          ll[] <- log(ifelse(PREFS, prob, 1-prob))
                          ## ll <- log(prob*PREFS + (1-prob)*(1-PREFS))
                        })),
                      defP(cell = "DM", cix = ixDM, cix_dim = 1), 
                      defP(cell = "RP", cix = ixRP, cix_dim = 1)))
  
  M$DATA$do.mc_ll <- T
  M$DM$mixin(adHST)
  M$DM$do.adapt <- T
  M$DM$do.group_scale <- T
  M$DM$do.mc_scale <- T
  M$DM$dilate_scale <- 2
  M$`*`$do.mc_ll <- T
  ndpar <- length(DM) - 1L
  ## Removing this (i.e. initializing all parameters to 1) has a minimal effect
  ## on estimated.
  ## M$DM$st <- t(replicate(nrDM, c(runif(ndpar, .1, 8), runif(1, 5, 50))))

  cat("\n")
  cat(sprintf("-- ADAPTIVE WARM UP (WITH SCALING) %d ITERATIONS ...\n", nr_sims_warm))
  update(M, nr_sims_warm)

  cat(sprintf("-- ADAPTIVE WARM UP (NO SCALING) %d ITERATIONS ...\n", 2*nr_sims_warm))
  M$DM$do.group_scale <- F
  update(M, 2L*nr_sims_warm)
  
  cat(sprintf("-- SIMULATING %d ITERATIONS ...\n", nr_sims))
  M$DM$do.adapt <- F
  update(M, nr_sims)
  M
}

buildData.rule_mixture_DM <- function(parent = data,
                                      ## I term? (consumption utility)
                                      I_nonLSKR = TRUE, I_LSKR = TRUE, 
                                      ## Use Weighting? 
                                      W_nonLSKR = TRUE, W_LSKR = TRUE, 
                                      ## fit KR, LS models?
                                      LS = TRUE, KR = TRUE,
                                      Qexclude = NULL) {

  model_data <- new.env(parent = parent)

  eval(bquote({
    W_LSKR <- .(W_LSKR)
    W_nonLSKR <- .(W_nonLSKR)
    I_LSKR <- .(I_LSKR)
    I_nonLSKR <- .(I_nonLSKR)
    KR <- .(KR)
    LS <- .(LS)
    Qexclude <- .(Qexclude)
  }), model_data)
  
  with(model_data, {

    ## subset for Q-CV
    if (!is.null(Qexclude)) {
      whichQ <- !ixQ %in% Qexclude
      ixQ <- ixQ[whichQ]
      ixDM <- ixDM[whichQ]
      xa <- xa[whichQ, ]
      xb <- xb[whichQ, ]
      pa <- pa[whichQ, ]
      pb <- pb[whichQ, ]
      PREFS <- PREFS[whichQ]
    }
    Nr <- length(ixQ)

    ## SIMPLE RULE BASED RPs (prospect specific)
    ## mnames <- c("zero", "maxMin", "minMax", "aveExpectation", "x_at_max_p")
    mnames <- c("zero", "maxMin", "minMax", "x_at_max_p")

    ## consumption utility? (I term)
    Is <- c(zero = I_nonLSKR,
            maxMin = I_nonLSKR,
            minMax = I_nonLSKR,
            x_at_max_p = I_nonLSKR)

    ## weighting function?
    Ws <- c(zero = W_nonLSKR,
            maxMin = W_nonLSKR,
            minMax = W_nonLSKR,
            x_at_max_p = W_nonLSKR)
    
    rpx <- as.matrix(dQ[mnames])
    RPXa <- RPXb <- cbind(c(rpx[ixQ, ]), NA, NA, NA) ## rule based mixture model
    RPPa <- RPPb <- cbind(rep.int(1, Nr*length(mnames)), 0, 0, 0)

    ## X_AT_MAX_P (within prospect)
    ## mnames <- c(mnames, "x_at_max_p")
    ## Is <- c(Is, x_at_max_p = 1)
    ## RPXa <- rbind(RPXa, cbind(dQ$x_at_max_p_A[ixQ], NA, NA, NA))
    ## RPXb <- rbind(RPXb, cbind(dQ$x_at_max_p_B[ixQ], NA, NA, NA))
    ## RPPa <- rbind(RPPa, cbind(rep.int(1, Nr), 0, 0, 0))
    ## RPPb <- rbind(RPPb, cbind(rep.int(1, Nr), 0, 0, 0))

    if(LS){
      ## LOOMES-SUGDEN 86
      mnames <- c(mnames, "LS")
      Is <- c(Is, LS = I_LSKR)
      Ws <- c(Ws, LS = W_LSKR)
      RPXa <- rbind(RPXa, cbind(dQ$meanA[ixQ], NA, NA, NA))
      RPXb <- rbind(RPXb, cbind(dQ$meanB[ixQ], NA, NA, NA))
      RPPa <- rbind(RPPa, cbind(rep.int(1, Nr), 0, 0, 0))
      RPPb <- rbind(RPPb, cbind(rep.int(1, Nr), 0, 0, 0))
    }

    if(KR){
      ## KOSZEGI-RABIN
      mnames <- c(mnames, "KR")
      Is <- c(Is, KR = I_LSKR)
      Ws <- c(Ws, KR = W_LSKR)
      RPXa <- rbind(RPXa, as.matrix(dQ[xaNames])[ixQ, ])
      RPXb <- rbind(RPXb, as.matrix(dQ[xbNames])[ixQ, ])
      RPPa <- rbind(RPPa, as.matrix(dQ[paNames])[ixQ, ])
      RPPb <- rbind(RPPb, as.matrix(dQ[pbNames])[ixQ, ])
    }
    
    meanDiff <- dQ$meanA[ixQ] - dQ$meanB[ixQ]
    
    nrRules <- length(mnames)

    ## Accross Decision Makers (recicled)
    nrRP <- nrDM
    ixRP <- ixDM

    ## Accross Questions
    ## nrRP <- nrQ
    ## ixRP <- ixQ

    stopifnot(identical(mnames, names(Is)))
    stopifnot(identical(mnames, names(Ws)))

  })
  model_data
}
