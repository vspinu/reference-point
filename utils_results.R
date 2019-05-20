
### DEFAULT GLOBALS
## !! All these globals are reset in report.Rnw

IND_IXS <- c(17, 50, 100)
MAJ_THRESHOLD<-50
START <- 1000
THIN <- 1

PAR_NAMES <- c("u", "w", "l", "xi")
PAR_LATEX <- c("$\\alpha$", "$\\gamma$", "$\\lambda$", "$\\xi$")
PAR_EXPR <- expression(alpha, gamma, lambda, xi)
PAR_RANGE <- expression(alpha[1] - alpha[139], gamma[1] - gamma[139],
                        lambda[1] - lambda[139], xi[1] - xi[139])

mc_mean_population <- function(M) {
  mpop <- window(meanLNorm(M$prDM)[, PAR_NAMES], start = START, thin = THIN)
  colnames(mpop) <- paste0("mu[", PAR_EXPR, "]")
  mpop
}

mc_sd_population <- function(M) {
  spop <- window(sigmaLNorm(M$prDM)[, PAR_NAMES], start = START, thin = THIN)
  colnames(spop) <- paste0("sigma[", PAR_EXPR, "]")
  spop
}

mc_ind <- function(M) {
  window(M$DM$mc_st, start = START, thin = THIN)  
}

mc_RP_pop <- function(M) {
  mn <- window(M$prRP$mc_st, start = START)
  colnames(mn) <- MNAMES
  mn
}

mc_RP_ind <- function(M) {
  window(M$RP$mc_st, start = START)
}

prLUW <- function(M, prefix = "tmp"){
  ## Posterior densities of behavioral parameters in population
  mpop <- mc_mean_population(M)
  spop <- mc_sd_population(M)
  densplot1(cbind(mpop, spop))
  ## pop_stack <- stack(as.data.frame(pop))
  ## ggplot(pop_stack) + geom_density(aes(x = values)) +
  ## facet_grid(ind~., labeller = label_parsed, scales = "free")

  sum <- cbind("$\\mu$" = median_summary(mpop)[, "Median"],
               "$\\sigma$" = median_summary(spop)[, "Median"])
  rownames(sum) <- PAR_LATEX

  pxtable(sum,
          caption = "Posterior point estimates of behavioral parameters in population: $\\alpha$ - exponent of power utility,  $\\gamma$ - parameter of Prelec weighting function, $\\lambda$ - loss aversion.", 
          label = paste0(prefix, "_prLUW_summaries"),
          scf = identity)
}

stats_ind <- function(M, prefix = "tmp", subject = 17){
  ## Posterior densities of behavioral parameters for subject 17 ($B_{17}$)
  ## mn <- window(M$prLUR, start = 10000, thin = 1)
  mn <- mc_ind(M)
  mn <- mn[, paste0(PAR_NAMES, subject)]
  colnames(mn) <- PAR_EXPR
  densplot1(mn)
  colnames(mn) <- PAR_LATEX
  pxtable(median_summary(mn),
          caption = "Posterior summaries for subject 17.", 
          label = paste0(prefix, "_summaries_ind"))
}

ind_spec_all <- function(M, prefix = "tmp"){
  ## Histograms of point estimates of behavioral parameters of all 139 subjects.
  ## beh <- colMeans(M$DM[["mc_st"]])[, PAR_NAMES]
  DM <- M$DM[["mc_st"]]
  beh <- apply(DM[START:dim(DM)[[1]], , ], c(2, 3), median)[, PAR_NAMES]
  colnames(beh) <- PAR_RANGE

  beh_stack <- stack(as.data.frame(beh))
  names(beh_stack)[[2]] <- "luw"

  print(ggplot(data = beh_stack) +
          geom_histogram(aes(x = values)) + 
          ## geom_density(aes(x = values, color = LUW)) + 
          facet_grid(~luw, scales = "free", labeller = label_parsed))
  ## scale_fill_discrete(name=expression(B[17]),
  ##                     labels = )

  colnames(beh) <- PAR_LATEX
  pxtable(summary(as.mcmc(beh))$quantiles,
          caption = "Quantiles of point estimates of behavioral parameters of all subjects.",
          label = paste0(prefix, "_ind_spec_all"))
}

RP_dens <- function(M, prefix = "tmp") {
  ## Posterior densities for RP rules in population
  mn <- mc_RP_pop(M)
  out <- stack(as.data.frame(mn))
  out$ind <- ordered(out$ind, levels = MNAMES)

  print(ggplot(aes(x = values), data = out) +
          geom_density(aes(fill = ind), color = "white") +
          scale_fill_brewer(palette="Set2", name = "Reference Point Rules") + 
          scale_y_continuous("density", limits = c(-1, 25), oob = scales::squish) +
          scale_x_continuous("probability") +
          facet_grid(ind ~ .))
  
  ## print(ggplot(aes(x = values, color = ind), data = out,) +
  ##       geom_density() +
  ##       scale_color_discrete(name = "Reference Point Rule") + 
  ##       scale_y_continuous("density", limits = c(-1, 25), oob = scales::squish) +
  ##       scale_x_continuous("probability"))

  pxtable(median_summary(mn), 
          caption = "Point estimates of RP mixture in population.",
          label = paste0(prefix, "_RP_tbl"), 
          srf = NULL)
}

RP_ind <- function(M, prefix = "tmp"){
  ## Posterior probability of a selection of subjects using a particular model
  
  mn <- mc_RP_ind(M)[, IND_IXS]

  models <- data.table(stack(as.data.frame(mn)))
  osubjects <- unique(models$ind)
  subjects <- gsub("RP", "", osubjects)
  subjects <- paste("Subject ", subjects[order(as.numeric(subjects))])
  models[, subject := factor(ind, levels = osubjects, labels = subjects, ordered = T)]
  models[, model := factor(MNAMES[values], levels = MNAMES, ordered = T)]

  models_aggr <-  models[, N := .N,  by =  .(subject)][, .(prop = .N/N[1]),  by = .(subject, model)]

  print(ggplot(models_aggr, aes(x = model, fill = model)) + 
          geom_bar(aes(y = prop), stat = "identity") +
          scale_fill_brewer(palette="Set2", name = "Reference Point Rules") + 
          facet_wrap( ~ subject) + 
          scale_x_discrete("", drop=FALSE) +
          scale_y_continuous("posterior probability") + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1)),
                strip.text = element_text(size = rel(1))))

}

get_rp_props_tbl <- function(M, start = START) {
  K <- length(MNAMES)
  tbl <- aaply(window(M$RP$mc_st, start), 2, tabulate, nbins = K)
  tbl <- round(tbl/rowSums(tbl), 4)
  colnames(tbl) <- MNAMES
  ## tbl[order(tbl[, 1], decreasing = T), ]*100
  tbl*100
}

RP_maj_support <- function(M, prefix = "tmp"){
  tbl <- get_rp_props_tbl(M, start = START)
  sharps <- unlist(alply(tbl, 1,
                         function(x){
                           if (any(x > MAJ_THRESHOLD))
                             MNAMES[which.max(x)]
                           else "OTHERS"
                         }), use.names = F)
  out <- structure(rep.int(0L, length(MNAMES)), names = MNAMES)
  tbl <- table(sharps[sharps != "OTHERS"])
  out[names(tbl)] <- tbl
  out <- data.frame(Nr = out)
  pxtable(out, label = paste0(prefix, "_maj_support"),
          caption = "Sharp predictions. \\small Number of subjects for which the posterior probability of the respective reference point is higher than 0.5.", 
          srf = NULL, scf = NULL)
}

sharp_params <- function(M, prefix = "tmp") {
  
  tbl <- get_rp_props_tbl(M, start = START)

  sharps <- unlist(alply(tbl, 1,
                         function(x){
                           if (any(x > MAJ_THRESHOLD))
                             MNAMES[which.max(x)]
                           else "OTHERS"
                         }), use.names = F)

  ## beh <- colMeans(M$DM[["mc_st"]])[, PAR_NAMES]
  DM <- M$DM[["mc_st"]]
  beh <- apply(DM[START:nrow(DM),, ], c(2, 3), median)[, PAR_NAMES]

  colnames(beh) <- PAR_EXPR
  beh1 <- data.table(cbind(as.data.table(beh), model = sharps))
  m <- function(x) round(median(x), 2)

  out <-
    rbind(
      beh1[, {
        c(map(.SD, m), Nr = .N)
      }, by = "model"],
      as.data.table(c(model = "ALL",
                      map(as.data.frame(beh), m),
                      Nr = nrow(beh))),
      use.names = T)
  
  ## mods <- model_names(as.character(out$model))

  setkey(out, model)
  allnms <- c(MNAMES, "OTHERS", "ALL")
  out <- out[allnms[allnms %in% model]]
  
  colnames(out) <- c("model", PAR_LATEX, "Nr")
  pxtable(out, label = paste0(prefix, "_sharp-params"),
          caption = "Median individual level parameters for each sharply clasified group.",
          scf = identity, srf = NULL)
}


est_ref_prob <- function(M, prefix = "tmp"){
  ## Level plot of posterior classification probabilities. Intensity is proportional to probability.
  tbl <- get_rp_props_tbl(M, start = START)
  rownames(tbl) <- gsub("RP", "", rownames(tbl))
  names(dimnames(tbl)) <- NULL

  pxtable(tbl, label = paste0(prefix, "_est_ref_prob"),
          caption = "Estimated parameters of the reference point mixture for each subject.", 
          tabular.environment = 'longtable', floating = FALSE)
  library(lattice)
  levelplot(tbl, xlab = "Subject", ylab = "Type",
            scales = list(x = list(rot = 90, cex = .5)))
}

## RElied on global which_SQ
## which_SQ <- which(tbl[, "Status Quo"] > MAJ_THRESHOLD)
## NOTE: doesn't work for 2 parameter weighting function
EU_tests <- function(M, subj_ixs = NULL, sign_level = .05, prefix = "tmp"){

  library(Hmisc)

  DM <- M$DM[["mc_st"]]
  DM <- DM[START:nrow(DM),, ]

  if (is.null(subj_ixs))
    subj_ixs <- 1:dim(DM)[[2]]
  
  ## one sided tests
  pil <- data.table(apply(DM[, subj_ixs, , drop = F], c(2, 3), quantile, probs = sign_level))
  pih <- data.table(apply(DM[, subj_ixs, , drop = F], c(2, 3), quantile, probs = 1-sign_level))

  alpha <- factor(ifelse(pil$u > 1, "$>1$", ifelse(pih$u < 1,  "$<1$", "$1$")), levels = c("$<1$", "$1$", "$>1$"))
  gamma <- factor(ifelse(pil$w > 1, "$>1$", ifelse(pih$w < 1,  "$<1$", "$1$")), levels = c("$<1$", "$1$", "$>1$"))

  out <- table(alpha, gamma)

  label <- paste0(prefix, "_EU_", sign_level)
  file <- sprintf("%s/%s.tex", get_tables_dir(), label)
  
  latex(addmargins(out),
        label = label,
        caption = paste0("EU test for individuals sharply classified as using Status Quo rule. ",
                         "Confidence level $", sign_level, "$\\%"),
        file = file, 
        title = "",
        cgroup=c("$\\gamma$"),
        n.cgroup=c(4),
        rgroup = c("$\\alpha$"), 
        n.rgroup=rep(4))

  cat(sprintf("\\input{table/%s}\n", label))
}
