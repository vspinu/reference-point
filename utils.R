suppressPackageStartupMessages({
  library(methods)
  library(xtable)
  library(ggplot2)
  library(purrr)
  library(magrittr)
  library(data.table)
  library(dtplyr)
  library(dplyr)
})

### DATA BUILDERS
buildData.RDU <- function(parent=data, RP){
    envir <- new.env(hash=T, parent=parent)
    envir$RP_type <- RP
    with(envir, {
        RP <- dD[[RP_type]]
        if(is.null(RP)) stop("RP mismatch")
        xa <- xa-RP
        xb <- xb-RP
        pa[is.na(pa)] <- 0L
        pb[is.na(pb)] <- 0L
        ## xa[is.na(xa)] <- 0L
        ## xb[is.na(xb)] <- 0L
        dpa <- cbind(1-aaply(pa, 1, cumsum)+pa, 0)
        dpb <- cbind(1-aaply(pb, 1, cumsum)+pb, 0)
    })
    envir
}

buildData.PT <- function(envir = new.env(hash=T, parent=globalenv()), RP){
    envir <- new.env(parent = envir)
    buildData.PT
    envir$RP_type <- RP
    with(envir, {
        rpx <- matrix(c(dQ[[RP_type]]))
        rpp <- matrix(rep.int(1, nrQ))
        rpX <- matrix(rpx[ixQ])
        rpP <- matrix(rep.int(1, Nr))
        if(is.null(RP)) stop("RP mismatch")
        pa[is.na(pa)] <- 0L
        pb[is.na(pb)] <- 0L
        dpa <- cbind(1-t(apply(pa, 1, cumsum))+pa, 0)
        dpb <- cbind(1-t(apply(pb, 1, cumsum))+pb, 0)
        cpa <- cbind(0, t(apply(pa, 1, cumsum)))
        cpb <- cbind(0, t(apply(pb, 1, cumsum)))
        xap <- xan <- xa - RP
        dpa[, 1:4][xap < 0] <- 0L
        cpa[, 2:5][xap > 0] <- 0L
        xap[xap < 0] <- 0L
        xan[xan > 0] <- 0L
        xan <- abs(xan)
        xbp <- xbn <- xb - RP
        dpb[, 1:4][xbp < 0] <- 0L  ## replace with NA ..to see if quickens the computation!!! it goes slower:(
        cpb[, 2:5][xbp > 0] <- 0L  ## doto : xbp > 0||is.na(xbp)
        xbp[xbp < 0] <- 0L
        xbn[xbn > 0] <- 0L
        xbn <- abs(xbn)
        ## xbn[is.na(xbn)] <- 0L
        ## xan[is.na(xan)] <- 0L
        ## xbp[is.na(xbp)] <- 0L
        ## xap[is.na(xap)] <- 0L
    })
    envir
}

buildData.simple <- function(envir = globalenv()){
    envir <- new.env(hash = T, parent = envir)
    with(envir, {
        xa <- xa
        xb <- xb
        pa[is.na(pa)] <- 0
        pb[is.na(pb)] <- 0
    })
    envir
}


### PREFS
gen_prefs <- function(u = .8, l = 2.4, r = 1, beta = 10, data = buildData.SPT(data, "minMax")){
    ## Assign probs and PREFS into DATA and return DATA.
    ## DATA must contain xa, xb, rpX and rpP (usually from a call to buildData.XXX)
    Nr <- nrow(get("xb", data))
    U <- rep.int(u, Nr)
    L <- rep.int(l, Nr)
    R <- rep.int(r, Nr)
    beta <- rep.int(beta, Nr)
    ## if(RP_type == "KR"){
    ##     prefs <-
    ##         c_SPT_pow(data$xb, data$pb, data$xb, data$pb, U, L) >
    ##           c_SPT_pow(data$xb, data$pb, data$xa, data$pa, U, L)
    ##     data$PREFS <- as.integer(prefs)
    ##     data <- buildData.SPT(data, RP_type)
    ## }
    delta <-
        c_SPT_pow_prelec1(data$xb, data$pb, data$rpX, data$rpP, L, U, R) -
          c_SPT_pow_prelec1(data$xa, data$pa, data$rpX, data$rpP, L, U, R)
    expdb <- exp(delta*beta)
    data$probs  <- expdb/(1+expdb)
    data$PREFS <- as.integer(runif(Nr) < data$probs)
    ## data$PREFS <- as.integer(data$probs > .5)
    data
}

proximity <- function(probs, PREFS){
    sqrt(mean((probs - PREFS)^2, na.rm = T))
}

eval_proximity <- function(PREFS, RP_type = "minMax", u = .8, l = 2.4, beta = 10){
    ## evaluate how close PREFS are to the model, with beta == 100 it is an
    ## almost perfect DM with beta = 1. It is a random noise. proximity of .20
    ## is exceptional and .40 is quite good. Proximity of .40 is enough to
    ## recover the model.
    if(is.environment(PREFS)){
        with(PREFS, proximity(probs, PREFS))
    }else{
        probs <- gen_prefs(u, l, beta, data = buildData.SPT(data, RP_type))$probs
        proximity(probs, PREFS)
    }
}

kr_vals <- function(data, u = .8, l = 2.4, table = T){
    N <- nrow(data$pa)
    U <- rep_len(u, N)
    L <- rep_len(l, N)
    vals <- cbind(
        aa = c_SPT_pow(data$xa, data$pa, data$xa, data$pa, U, L),
        bb = c_SPT_pow(data$xa, data$pa, data$xa, data$pa, U, L),
        ba = c_SPT_pow(data$xb, data$pb, data$xa, data$pa, U, L),
        ab = c_SPT_pow(data$xa, data$pa, data$xb, data$pb, U, L))
    if(table) table(apply(vals, 1, which.max))/nrDM
    else vals
}



### OTHER
meanLNorm <- function(pobj,...){
    mc <- getMC(pobj, ..., transform_coda=F)
    out <- as.mcmc(exp(mc[,,"meanlog"] + 1/(mc[,,"taulog"]*2)))
    if(!is.array(out))
        out <- as.mcmc(array(out,  dim = c(length(out), 1L)))
    dimnames(out) <- dimnames(mc)[1:2]
    out
}

sigmaLNorm <- function(pobj,...){
    mc <-
        if(is.matrix(pobj)){
            out <- as.array(exp(pobj[,"meanlog"] + 1/(pobj[,"taulog"]*2)))
            rownames(out) <- rownames(pobj)
        }else{
            mc <- getMC(pobj, ..., transform_coda=F)
            out <- as.mcmc(sqrt((exp(1/(mc[,,"taulog"]))-1)*exp(2*mc[,,"meanlog"] + 1/(mc[,,"taulog"]))))
            if(!is.array(out))
                out <- as.mcmc(array(out,  dim = c(length(out), 1L)))
            dimnames(out) <- dimnames(mc)[1:2]
        }
    out
}

Q_and_density <- function(gr_df, formula= ~ sims|Qid, main=gr_df$gr[1], groups=rep(1L, nrow(gr_df)),...){
                                        # gr_df should have column sims and Qid
                                        # Qset should be accesable
    print(
          densityplot(formula, data=gr_df, main=paste("Question Group : ", main), subscripts=T, ylim=c(0, 30), xlim=c(-.05, 1),
                      lattice.options=lattice.options(default.args = list(as.table = TRUE)), groups=groups,...,
                      panel=function(subscripts, groups, ...){
                          Q <- QSet[[gr_df[subscripts[1], ]$Qid]]
                          y <- 20
                          ydif <- 2
                          probs_cex <- 2
                          panel.points(x=Q[[1]]@outcomes,y=y, pch=19, col="red", cex=sqrt(Q[[1]]@probs)*probs_cex)
                          panel.lines(x=Q[[1]]@outcomes,y=y, col="red")
                          panel.abline(v=Q[[1]]@outcomes, col="red", lwd=.8, lty=4)
                          panel.text(x=Q[[1]]@outcomes,y=y, adj=c(0, -1), labels=Q[[1]]@probs, col="red", cex=.6)
                          panel.points(x=Q[[2]]@outcomes,y=y-ydif, pch=19, col="blue", cex=sqrt(Q[[2]]@probs)*probs_cex)
                          panel.lines(x=Q[[2]]@outcomes,y=y-ydif, col="blue")
                          panel.abline(v=Q[[2]]@outcomes, col="blue", lwd=.8, lty=4)
                          panel.text(x=Q[[2]]@outcomes,y=y-ydif, adj=c(0, 2), labels=Q[[2]]@probs, col="blue", cex=.6)
                                        #panel.grid()
                          panel.densityplot(jitter.amount=2,pch=".", groups=groups, subscripts=subscripts, ...)
                          panel.text(.6, 28, labels=paste("Qid= ", Q@traits$Qid, " ",
                                               "Rep= ", Q@traits$rep, " ",
                                               "Shift= ", Q@traits$shift),
                                     cex=.6)
                      },
                      auto.key=list(space="right")
                      ## key = list(lines = Rows(trellis.par.get("superpose.line"), 1:2),
                      ##   text = list(lab = as.character())), title = "Model conf:")
                      ))
}

init_state_Mu <- function(M, expr = expression()){
    S_par <- M$LUR$rv_dimnames[[1]]
    if("beta"%in%S_par){
        M$prLUR$st["beta", ] <- c(runif(1, .8, 1.5), runif(1, 0, 5))
        M$LUR$st[,"beta" ] <- runif(nrDM, 10, 40)
        M$LUR$st[ , S_par[!S_par %in% "beta"]] <- runif((length(S_par) - 1)*nrDM, 0, 3)
    }else{
        M$LUR$st[ , S_par[!S_par %in% "beta"]] <- runif((length(S_par))*nrDM, 0, 3)
    }
    if("U"%in%S_par)
        M$prLUR$st["U", ] <- c(runif(1, 1, 4), runif(1, -2, 1))
    if("R"%in%S_par)
        M$prLUR$st["R", ] <- c(runif(1, 1, 4), runif(1, -2, 1))
    if("L"%in%S_par)
        M$prLUR$st["L", ] <- c(runif(1, 1, 4), runif(1, -2, 1))
    M$prLUR$ll[] <- -Inf
    M$LUR$ll[] <- -Inf
    eval(expr)
}

init_state_M2 <- function(M, expr = expression()){
    S_par <- M$Lb$rv_dimnames[[1]]
    if("beta"%in%S_par){
        M$prLb$st["beta", ] <- c(runif(1, .8, 1.5), runif(1, 0, 5))
        M$Lb$st[,"beta" ] <- runif(nrDM, 10, 40)
        M$Lb$st[ , S_par[!S_par %in% "beta"]] <- runif((length(S_par) - 1)*nrDM, 0, 3)
    }else{
        M$Lb$st[ , S_par[!S_par %in% "beta"]] <- runif((length(S_par))*nrDM, 0, 3)
    }
    if("U"%in%S_par)
        M$prLb$st["U", ] <- c(runif(1, 1, 4), runif(1, -2, 1))
    if("R"%in%S_par)
        M$prLb$st["R", ] <- c(runif(1, 1, 4), runif(1, -2, 1))
    if("L"%in%S_par)
        M$prLb$st["L", ] <- c(runif(1, 1, 4), runif(1, -2, 1))
    S_par <- M$UR$rv_dimnames[[1]]
    if("beta"%in%S_par){
        M$prUR$st["beta", ] <- c(runif(1, .8, 1.5), runif(1, 0, 5))
        M$UR$st[,"beta" ] <- runif(nrDM, 10, 40)
        M$UR$st[ , S_par[!S_par %in% "beta"]] <- runif((length(S_par) - 1)*nrDM, 0, 3)
    }else{
        M$UR$st[ , S_par[!S_par %in% "beta"]] <- runif((length(S_par))*nrDM, 0, 3)
    }
    if("U"%in%S_par)
        M$prUR$st["U", ] <- c(runif(1, 1, 4), runif(1, -2, 1))
    if("R"%in%S_par)
        M$prUR$st["R", ] <- c(runif(1, 1, 4), runif(1, -2, 1))
    if("L"%in%S_par)
        M$prUR$st["L", ] <- c(runif(1, 1, 4), runif(1, -2, 1))
    M$prUR$ll[] <- -Inf
    M$UR$ll[] <- -Inf
    M$prLb$ll[] <- -Inf
    M$Lb$ll[] <- -Inf
    eval(expr)
}

meltRP <- function(M, type=NULL, start= 1){
    rpMC <- data.frame(drop(M$RP[["mc_st"]])[-(1:start),], check.names=F)
    rp_m <- rename(stack(rpMC), c(ind="Qid", values="sims"))
                                        #rp_m$figs <- cut(as.numeric(rp_m$Qid),  8, labels=1:8)
    rp_m <- merge(dQ[c("Qid", "gr", "rb_type")], rp_m,  by="Qid")
    rp_m$type <- if(is.null(type)) .type(M, FALSE)
    else type
    rp_m
}

run_models <- function(M, nr = 10, nr_iter = 1000 , thin = 1, start = 1,
                       file = paste(deparse(substitute(M)), nr, ".pdf", sep = "_"),
                       expr = expression()){
    pdf(file, w = 14, h = 10)
    mcs <- lapply(1:nr, function(i){
        init_state_Mu(M, expr = expr)
        update(M, nr_iter = nr_iter, thin = thin, reinit = T)
        plot(meanLNorm(M$prLUR, start = start))
        meltRP(M, i, start = start)
    })
                                        #mcs <- list(meltRP(M, start=700), meltRP(tM, start = 100))
    df <- do.call(rbind, mcs)
    d_ply(df, .(gr),  function(df_cut) Q_and_density(df_cut, groups = df_cut$type))
    dev.off()
    system(paste("open ", file))
    mcs
}

run_models2 <- function(M, nr = 10, nr_iter = 1000, start = 1, thin = 1,
                        file = paste(deparse(substitute(M)), nr, ".pdf", sep = ""),
                        expr = expression()){
    pdf(file)
    mcs <- lapply(1:nr, function(i){
        init_state_M2(M, expr = expr)
        update(M, nr_iter = nr_iter, thin = thin, reinit = T)
        plot(meanLNorm(M$prUR, start = start))
        plot(meanLNorm(M$prLb, start = start))
        meltRP(M, i, start = start)
    })
                                        #mcs <- list(meltRP(M, start=700), meltRP(tM, start = 100))
    df <- do.call(rbind, mcs)
    d_ply(df, .(gr),  function(df_cut) Q_and_density(df_cut, groups = df_cut$type))
    dev.off()
    system(paste("open ", file))
    mcs
}

Qplot <- function(..., start = 1){
### plot the RP densities toghetehre with the choice loterry
### ... are the models to plot
    mcs <- lapply(list(...), meltRP, start = start)
    df <- do.call(rbind, mcs)
    d_ply(df, .(gr),  function(df_cut) Q_and_density(df_cut, groups = df_cut$type))
}


## Q testing, common_beta
crosval <- function(M_function, allFolds, testFolds, nrsim = 2000,
                    mnames = model_names, do.pc_ll = "fold",
                    expr = expression()){
    models <- list()
    predict_ll <- list()
    for (mname in mnames){
        M <- do.call(M_function, list(mname))
        models[[mname]] <- M
        M$DATA$folds <- allFolds
        M$DATA$do.pc_ll <- do.pc_ll
        eval(expr)
        M$initTestMirror(testFolds)
        cat("\n\n****** updating ", mname, substitute(M_function), "********\n")
        update(M, nr_iter=nrsim, thin=1, app=T)
        cat("\n****** predicting ", mname, "********\n")
        M$DATA$predict()
        predict_ll[[mname]] <- M$DATA$pc_ll
    }

    save(models, predict_ll, foldQ, testFolds,
         file =
         sprintf("%s_Q_%s_%s.rda",
                 paste(testFolds, collapse = ":"),
                 deparse(substitute(M_function)), nrsim))
    invisible(models)
}

update_kfold <- function(mname, crosfolds, M_function, allFolds,
                         nrsim = 2000, use_crosfolds_name = T){
    ## crosfolds <-  an integer vector of folds
    M <- M_function(mname)
    M$DATA$folds <- allFolds
    lev <- levels(allFolds)
    for(i in levels(crosfolds) ){
        cat("setting mirror ", i, "\n")
        M$initTestMirror(lev[crosfolds == i], test_mname = if(use_crosfolds_name) i)
    }
    predict_ll <- list()
    for(i in seq_along(M$mirrors_test)){
        M$mirror <- fname <- M$mirrors_test[[i]]
        cat("\n\n****** updating ", mname, substitute(M_function),
            " mirror ", i, fname, "********\n")
        M$DATA$do.pc_ll <- "fold"
        M$DATA$do.mc_ll <- T
        update(M, nr_iter=nrsim, thin=1, app=T)
        cat("\n****** predicting ", mname, " fold ", i, "(", fname, ") ********\n")
        M$DATA$predict()
    }
    nrlev <- length(levels(crosfolds))
    save(M, allFolds, crosfolds,
         file = sprintf("crosQ_%sfold_%s_%s_%s.rda", nrlev, 
           mname, deparse(substitute(M_function)), nrsim))
    M
}

crosdiffs <- function(pc_sims, start = 500){
    names <- names(pc_sims)
    delta <- list()
    for(nm1 in names)
        for(nm2 in names)
            if(nm1 != nm2){
                df <- cbind(data.frame(mod1 = nm1, mod2 = nm2),
                            as.data.frame((pc_sims[[nm1]] - pc_sims[[nm2]])[-(1:start), ]))
                delta[[paste(nm1, nm2)]] <- df
            }
    delta
}
        
    

plot_crosdiffs <- function(models, each_ll = FALSE){
    predict_ll <- lapply(names(models), function(mn){
        models[[mn]]$DATA$pc_ll
    })
    names(predict_ll) <- names(models)

    library(ggplot2)
    nm <- names(predict_ll)
    delta <- list()

    for(i in seq_along(models))
        for( j in seq_along(models)){
            if(i != j){
                ll_diff <-
                    if(each_ll) predict_ll[[i]] - predict_ll[[j]]
                    else rowSums(predict_ll[[i]]) - rowSums(predict_ll[[j]])
                if(NROW(ll_diff) > 600) ll_diff <- as.matrix(ll_diff)[-(1:500), ]
                if(NROW(ll_diff) > 1200) ll_diff <- as.matrix(ll_diff)[-(1:1000), ]
                delta[[paste(i, j)]] <-
                    data.frame(mod1 = nm[[i]],  mod2 = nm[[j]], ll_diff = ll_diff,
                               check.names = F, stringsAsFactors = F)
            }}

    df <- do.call(rbind, delta)
    df <- transform(df,
                    mod1 = relevel(as.factor(mod1), "maxMin"),
                    mod2 = relevel(as.factor(mod2), "maxMin"))

    lapply(1:NCOL(ll_diff), function(i){
        df <- df[, c(1, 2, i + 2)]
        names(df) <- c("mod1", "mod2", "ll_diff")
        ggplot(df, aes(x = ll_diff)) + geom_histogram() +
            geom_vline(xintercept = 0, color = "red") + facet_grid(mod1~mod2) +
                ggtitle(paste("Subject: ", i))
    })
}

mm_attach <- function(src)
    attach(paste0("~/data/maxmin/", src))

mm_load <- function(src)
    load(paste0("~/data/maxmin/", src))

write_folds <- function(testFolds, file){
    cat("\n\n******** ", paste0(testFolds, collapse = ":"), "********\n",
        file = file, append = T)
    cat(capture.output(print(data$QSet[testFolds])), sep = "\n", 
        file = file, append = T)
}

rp_names <- function(DM){
    mnames <- c(zero = "Status Quo",
                maxMin = "MaxMin",
                minMax = "MinMax",
                x_at_max_p = "X at Max P",
                "X with Max P" = "X at Max P", 
                LS   = "Expected Value",
                KR = "Prospect Itself")
    if(is.environment(DM)) unname(mnames[get("mnames", DM)])
    else{
        out <- unname(mnames[DM])
        out[is.na(out)] <- DM[is.na(out)]
        out
    } 
}

get_post_prob_intervals <- function(M, start = 10000, signif = .05){
    M$DM$mc_st[, 23]
}
## smarter table

stable <- function(x, ..., names = NULL){
    tbl <- table(x, ...)
    tbl <- tbl/sum(tbl)
    if(is.null(names))
        tbl
    else{
        out <- structure(c(0, 0), names = names)
        out[names(tbl)] <- tbl
        out
    }
}


### reporting

get_tables_dir <- function() {
  root <- Sys.getenv("MODEL-DIR")
  if (root == "")
    root <- "."
  dir <- sprintf("%s/table/", root)
  message("Saving: ", dir)
  dir.create(dir, FALSE, TRUE)
  dir
}

pxtable <- function(obj,  caption = "", label = "", scf = NULL, srf = identity, floating = T, ...,
                    save = getOption("maxmin.paper.save.tables", FALSE)) {
  if (save) {
    lab <- paste0("tb:", label)
    tbl <- xtable(obj, caption = caption, label = lab)
    out <- capture.output(print(tbl,
                                sanitize.colnames.function = scf,
                                sanitize.rownames.function = srf,
                                floating = floating, 
                                ...))
    cat(paste(out, collapse = "\n"),
        file = sprintf("%s/%s.tex", get_tables_dir(), label))
    cat(sprintf("\\input{table/%s}\n", label))
  } else {
    print(obj)
  }
}


densplot1 <- function(mcmc){
    ## colnames must be expressions
    mcmc <- as.mcmc(mcmc)
    op <- par(mfrow = c(2, ceiling(ncol(mcmc)/2)), c(.1, .1, .1, .1))
    on.exit(par(op))
    for(n in colnames(mcmc)){
        densplot(mcmc[, n], main = parse(text = n))
    }
}

median_summary <- function(mcmc){
    names <- colnames(mcmc)
    out <- cbind(Median = sapply(names, function(n) median(mcmc[, n])),
                 SD = sapply(names, function(n) sd(mcmc[, n])))
    ## out <- sapply(names, function(n) median(mcmc[, n]))
    rownames(out) <- names
    out
}

gen_pop_results <- function(M, label,
                            caption_beh = "",
                            caption_rp = ""){
    caption_beh <- paste("Point estimates of behavioral parameters in population.", caption_beh)
    caption_rp <- paste("Point estimates of RP rules in population.", caption_rp)

    mnames <- get("mnames", M$DM)
    
    ## PARAMETERS IN POPULATION 
    mpop <- window(meanLNorm(M$prDM)[, par_names], start = START, thin = THIN)
    spop <- window(sigmaLNorm(M$prDM)[, par_names], start = START, thin = THIN)

    ## colnames(mpop) <- gsub("^\\$", "$\\\\mu_", par_latex)
    ## colnames(spop) <- gsub("^\\$", "$\\\\sigma_", par_latex)

    sum <- cbind("$\\mu$" = median_summary(mpop)[, "Median"],
                 "$\\sigma$" = median_summary(spop)[, "Median"])
    rownames(sum) <- par_latex
    pxtable(sum, 
            caption = caption_beh, 
            label = sprintf("prLUW_summaries_%s", label),
            scf = identity)

    ## RP POSTERIORS
    mn <- window(M$prRP$mc_st, start = START)
    colnames(mn) <- mnames

    pxtable(median_summary(mn),
            caption = caption_rp,
            label = sprintf("RP_tbl_%s", label), 
            srf = NULL)
}

encode_model <- function(seed, N, I_nonLSKR = T, I_LSKR = T, W_nonLSKR = T, W_LSKR = T,
                         LS = T, KR = T, fU = "pow", fW = "prelec1"){
  sprintf("md%s_%d_I%d%d_W%d%d_LS%d_KR%d_%s_%s",
          seed, N, I_nonLSKR, I_LSKR, W_nonLSKR, W_LSKR, LS, KR, fU, fW)
}

assignCommandArgs <- function(..., envir = .GlobalEnv, collapse = T){
  dots <- list(...)
  if(any(sapply(dots, is.null)))
    stop("default arguments cannot be NULL. Need to infer the type.")
  stopifnot(all(nzchar(names(dots))))
  args <- commandArgs(T)
  for(i in seq_along(args))
    dots[[i]] <- as(args[[i]], class(dots[[i]]))
  for(nm in names(dots)){
    assign(nm, dots[[nm]], envir = envir)
  }
  out <-
    if(collapse){
      paste(lapply(dots, as.character), collapse = "_")
    } else
      dots
  invisible(out)
}

stopif <- function (...) {
  n <- length(ll <- list(...))
  if (n == 0L) 
    return(invisible())
  mc <- match.call()
  for (i in 1L:n)
    if (!is.logical(r <- ll[[i]]) || anyNA(r) || any(r)) {
      ch <- deparse(mc[[i + 1]], width.cutoff = 60L)
      if (length(ch) > 1L)
        ch <- paste(ch[1L], "....")
      stop(sprintf("%s is TRUE", ch), call. = FALSE, domain = NA)
    }
  invisible()
}
