# tools for running standard error simulation studies

condNumLimit <- 1e7

calcCondNum <- function(hess) {
  d <- try(svd(hess, nu=0, nv=0)$d)
  if (is(d, "try-error")) return(1e16)
  if (all(d > 0)) {
    max(d)/min(d)
  } else {
    1e16
  }
}

MCphase <- function(modelGen, reps=500, verbose=TRUE, maxCondNum) {
  emcycles <- rep(NA, reps)
  condnum <- rep(NA, reps)
  est <- matrix()
  for (rep in 1:reps) {
    set.seed(rep)
    model <- modelGen()
    em <- model$compute
    getCondNum <- list(mxComputeOnce('fitfunction', 'information', 'meat'),
                       mxComputeReportDeriv())
    plan <- mxComputeSequence(c(em, getCondNum))
    model$compute <- plan
    fit <- try(mxRun(model, silent=TRUE, suppressWarnings = TRUE), silent=TRUE)
    if (inherits(fit, "try-error")) {
      print(fit)
      condnum[rep] <- 1e16
      next
    } else if (fit$output$status$code != 0) {
      print(paste("status code", fit$output$status$code))
      next
    }
    emstat <- fit$compute$steps[[1]]$output
    emcycles[rep] <- emstat$EMcycles
    condnum[rep] <- calcCondNum(fit$output$hessian)
    par <- omxGetParameters(fit)
    if (any(is.na(par))) {
      print(par)
      condnum[rep] <- 1e16
      next
    }
    if (verbose) print(paste(c(rep, emstat, round(condnum[rep])), collapse=" "))
    if (all(dim(est) == 1)) {
      est <- matrix(NA, length(par), reps)
      rownames(est) <- names(par)
    }
    est[,rep] <- par
  }
  list(condnum=condnum, est=est)
}

getMCdata <- function(name, modelGen, correct, recompute=FALSE, reps=500,
                      envir=parent.frame(), maxCondNum) {
  if (missing(maxCondNum)) stop("Provide a maxCondNum")
  
  correct <- c(correct)
  rda <- paste("data/", name, ".rda", sep="")

  if (!recompute) {
    if (file.exists(rda)) {
      load(rda, envir=envir)
    } else if (file.exists(paste("models/enormous/", rda, sep=""))) {
      load(paste("models/enormous/", rda, sep=""), envir=envir)
    } else {
      recompute <- TRUE
    }
  }
  
  if (recompute) {
    got <- MCphase(modelGen, reps, maxCondNum=maxCondNum)
    mcMask <- rep(TRUE, reps)
    if (!is.na(maxCondNum)) {
      mcMask <- !is.na(got$condnum) & got$condnum < maxCondNum
    }
    est <- got$est
    mcEst <- apply(est[,mcMask], 1, mean)
    bias <- mcEst - correct
    if (reps < length(correct)) stop("Not enough replications to estimate the Hessian")
    mcCov <- cov(t(est))
    #max(abs(apply(est, 1, sd) - sqrt(diag(mcCov))))
    mcHessian <- solve(mcCov/2)
    mcBias <- bias
    mcSE <- sqrt(2*diag(solve(mcHessian)))
    save(mcMask, mcBias, mcSE, mcHessian, file=rda)
    if (!is.na(maxCondNum)) {
      cat(paste("Note:", sum(!mcMask), "excluded due to condition number\n"))
    }
    cat("Monte-Carlo study complete. Proceeding with accuracy benchmark.\n")
#    stop("stopped here")
    load(rda, envir=envir)  # copy to parent environment
  }
}

# This was written in terms of the information matrix, but
# I switched it to check the parameter covariance matrix (invert everything).
mvn_KL_D <- function(invH, H) {
  pcm <- solve(mcHessian)
  .5*(tr(invH %*% pcm) - nrow(H) - log(det(pcm)/det(H)))
}

summarizeInfo1 <- function(condnum, emstat=list(EMcycles=NA, semProbeCount=NA),
                           H, standardErrors, cputime, method) {
  numReturn <- 6
  
  if (!is.na(condnum) && condnum > condNumLimit) return(rep(NA, numReturn))
  
  normH <- NA
  if (!is.null(H) && all(eigen(H, only.values =TRUE)$val > 0)) {
#    normH <- norm(H - mcHessian, "2")
#    normH <- norm(H %*% solve(mcHessian) - diag(nrow(H)), "2")
    # D_KL(H | mcHessian)   -- backwards
#    normH <- .5*(tr(solve(mcHessian) %*% H) - nrow(H) - log(det(H)/det(mcHessian)))
    # D_KL(mcHessian | H)
    iH <- try(solve(H), silent=TRUE)
    if (is(iH, "try-error")) return(rep(NA, numReturn))
    normH <- mvn_KL_D(H, iH)
  }
  
  normRd <- NA
  rd <- (standardErrors - mcSE) / mcSE
  if (!is.na(condnum)) {
    if (all(is.finite(rd))) {
      normRd <- norm(rd, "2")
    } else {
      print(paste("Method", method,"condition number", condnum, "but some SEs are NA"))
      condnum <- NA
    }
  }
  
  got <- c(cputime, emstat$EMcycles, emstat$semProbeCount, condnum, normH, normRd)

  if (length(got) != numReturn) {
    print('wrong length')
    print(got)
    return(rep(NA, numReturn))
  } else {
    return(got)
  }
}

summarizeInfo <- function(fitModel, method) {
  emstat <- list(EMcycles=NA, semProbeCount=NA)
  if (length(intersect(c('mr', 'tian', 'agile'), method))) {
    emstat <- fitModel$compute$steps[[1]]$output
  }
  
  if (fitModel$output$status$code != 0) {
    summarizeInfo1(NA, emstat, NULL, NULL,
                   fitModel$output$cpuTime, method)
    return()
  }

  H <- fitModel$output$hessian
  if (is.null(H)) H <- fitModel$output$ihessian
  condnum <- calcCondNum(H)

  H <- NULL
  if (!is.na(condnum) && condnum < 1e12) {
    if (!is.null(fitModel$output[['hessian']])) {
      H <- fitModel$output[['hessian']]
    }
    if (is.null(H) && !is.null(fitModel$output[['ihessian']])) {
      H <- solve(fitModel$output[['ihessian']])
    }
  }
  
  summarizeInfo1(condnum, emstat, H, fitModel$output$standardErrors,
                 fitModel$output$cpuTime, method)
}

summarizeDetail <- function(detail, maxCondNum=NA) {
  mask <- rep(TRUE, dim(detail)[3])
  if (!is.na(maxCondNum)) {
    mask <- apply(is.na(detail['condnum',,]) | detail['condnum',,] < maxCondNum, 2, all)
    detail <- detail[,,mask]
  }
  excluded <- 0
  if (dim(detail)[3] > 1) {
    excluded <- apply(detail['condnum',,], 1, function (c) sum(is.na(c)))
  }
  print(round(rbind(excluded, apply(detail, 1:2, mean, na.rm=TRUE)), 4))
  cat(paste("  N=", sum(mask), "\n", sep=""))
}

testPhase <- function(modelGen, reps = 500, verbose=TRUE, methods=c('agile', 'meat')) {
  rec <- c('cputime', 'emcycles', 'probes', 'condnum', 'hNorm', 'rdNorm')
  detail <- array(NA, dim=c(length(rec), length(methods), reps),
                  dimnames=list(rec, methods, NULL))
  
  for (rep in 1:reps) {
    warnings()
    set.seed(rep)
    model <- modelGen()
    em <- model$compute
    fit <- NULL    # needed for MLE
  
    fitfun <- c()
    if (is(em$mstep, "MxComputeSequence")) {
      fitfun <- sapply(em$mstep$steps, function(step) step$fitfunction)
    } else {
      fitfun <- em$mstep$fitfunction
    }
    
    sem <- intersect(c('mr', 'tian'), methods)
    if (length(sem)) {
      em$accel <- ""
      em$tolerance <- 1e-11
      em$maxIter <- 750L
      for (semType in sem) {
        em$information <- "mr1991"
        em$infoArgs <- list(fitfunction=fitfun, semMethod=semType, semTolerance=sqrt(1e-6))
        plan <- mxComputeSequence(list(
          em,
          mxComputeStandardError(),
          mxComputeReportDeriv()
        ))
        model$compute <- plan
        fit <- try(mxRun(model, silent=TRUE, suppressWarnings=TRUE), silent=TRUE)
        if (inherits(fit, "try-error")) {
          print(paste("error in", semType))
          print(fit)
          next
        } else if (fit$output$status$code != 0) {
          print(paste("status code", fit$output$status$code, "without acceleration"))
          break
        } else {
          detail[,semType,rep] <- summarizeInfo(fit, semType)
        }
      }
    }
    
    # need the MLE
    if (is.null(fit) || inherits(fit, "try-error")) {
      em$tolerance <- 1e-11
      model$compute <- em
      fit <- try(mxRun(model, silent=TRUE, suppressWarnings = TRUE), silent=TRUE)
      if (inherits(fit, "try-error")) {
        print(paste("error finding MLE"))
        print(fit)
        next
      } else if (fit$output$status$code != 0) {
        print(paste("status code", fit$output$status$code))
        next
      }
    }
    
    if (length(intersect(methods, "agile"))) {
      em$accel <- 'ramsay1975'
      em$tolerance <- 1e-11
      em$information <- "mr1991"
      em$infoArgs <- list(fitfunction=fitfun, semMethod="agile")
      plan <- mxComputeSequence(list(
        em,
        mxComputeStandardError(),
        mxComputeReportDeriv()
      ))
      if (is.null(fit)) fit <- model
      fit$compute <- plan
      # reuse the MLE, if possible
      fit <- try(mxRun(fit, silent=TRUE, suppressWarnings = TRUE), silent=TRUE)
      if (inherits(fit, "try-error")) {
        print(paste("error in agile"))
        print(fit)
        next
      } else if (fit$output$status$code != 0) {
        print(paste("status code", fit$output$status$code, "in agile"))
        next
      } else {
        detail[,"agile",rep] <- summarizeInfo(fit, "agile")
      }
    }
    
    if (length(intersect(methods, "meat"))) {
      meat <- mxModel(fit,
                      mxComputeSequence(steps=list(
                        mxComputeOnce('fitfunction', 'information', "meat"),
                        mxComputeStandardError(),
                        mxComputeReportDeriv())))
      meat <- mxRun(meat, silent=TRUE)
      detail[,"meat",rep] <- summarizeInfo(meat, "meat")
    }
    
    if (length(intersect(methods, "sandwich"))) {
      sandwich <- mxModel(fit,
                      mxComputeSequence(steps=list(
                        mxComputeOnce('fitfunction', 'information', "sandwich"),
                        mxComputeStandardError(),
                        mxComputeReportDeriv())))
      sandwich <- mxRun(sandwich, silent=TRUE)
      detail[,"sandwich",rep] <- summarizeInfo(sandwich, "sandwich")
    }

    if (length(intersect(methods, c("oakes")))) {
      em$information <- "oakes1999"
      em$infoArgs <- list(fitfunction=fitfun)
      plan <- mxComputeSequence(list(
        em,
        mxComputeStandardError(),
        mxComputeReportDeriv()
      ))
      fit$compute <- plan
      # reuse the MLE
      fit <- try(mxRun(fit, silent=TRUE, suppressWarnings = TRUE), silent=TRUE)
      if (inherits(fit, "try-error")) {
        print(paste("error in agile"))
        print(fit)
        next
      } else if (fit$output$status$code != 0) {
        print(paste("status code",fit$output$status$code,"in agile"))
        next
      } else {
        detail[,"oakes",rep] <- summarizeInfo(fit, "oakes")
      }
    }
    
    if (length(intersect(methods, "estepH"))) {  # should be Mstep, oops
      estepH <- mxModel(fit,
                        mxComputeSequence(steps=list(
                          mxComputeOnce(em$expectation, 'scores'),
                          mxComputeOnce(fitfun, 'information', "hessian"),
                          mxComputeStandardError(),
                          mxComputeReportDeriv())))
      estepH <- mxRun(estepH, silent=TRUE)
      detail[,"estepH",rep] <- summarizeInfo(estepH, "estepH")
    }
    
    if (length(intersect(methods, "re"))) {
      re <- mxModel(fit,
                    mxComputeSequence(steps=list(
                      mxComputeNumericDeriv(stepSize = 1e-3, iterations = 2),
                      mxComputeStandardError(),
                      mxComputeReportDeriv())))
      re <- mxRun(re, silent=TRUE)
      detail[,"re",rep] <- summarizeInfo(re, "re")
    }
    
    if (verbose) {
      summarizeDetail(detail)
    }
  }
  detail
}

quantifyAsymmetry <- function(info) {
  sym1 <- (info + t(info))/2
  sym2 <- try(chol(solve(sym1)), silent=TRUE)
  if (inherits(sym2, "try-error")) return(NA)
  asymV <- (info - t(info))/2
  norm(sym2 %*% asymV %*% sym2, type="2")
}

summarizeAgile <- function(fit) {
  numReturn <- 4
  
  condnum <- calcCondNum(fit$output$ihessian)
  if (is.null(condnum)) condnum <- NA
  if (is.na(condnum) || (!is.na(condnum) && condnum > condNumLimit)) return(rep(NA, numReturn))

  H <- fit$compute$steps[[1]]$debug$outputInfo
  if (is.null(H)) return(rep(NA, numReturn))
  
  # Jamshidian (2000) defined this in terms of the inverse Hessian
  # even though it seems to work regardless of the inverse.
  asym <- quantifyAsymmetry(solve(H))
  
#  max(abs((H + t(H))/2 - solve(fit$output[['ihessian']])))  # == 0
  H <- (H + t(H))/2
  
  normH <- NA
  if (!is.null(H) && all(eigen(H, only.values =TRUE)$val > 0)) {
    #    normH <- norm(H - mcHessian, "2")
    #    normH <- norm(H %*% solve(mcHessian) - diag(nrow(H)), "2")
    # D_KL(H | mcHessian)   -- backwards
    #    normH <- .5*(tr(solve(mcHessian) %*% H) - nrow(H) - log(det(H)/det(mcHessian)))
    # D_KL(mcHessian | H)
    iH <- try(solve(H), silent=TRUE)
    if (is(iH, "try-error")) return(rep(NA, numReturn))
    normH <- mvn_KL_D(H, iH)
  }
  
  normRd <- NA
  rd <- (fit$output$standardErrors - mcSE) / mcSE
  if (all(is.finite(rd))) {
    normRd <- norm(rd, "2")
  }
  
  c(condnum, asym, normH, normRd)
}

summarizeASEM <- function(detail) {
  excluded <- apply(detail[,'condnum',], 1, function (c) sum(is.na(c)))
  grid <- cbind(excluded,
                apply(detail, 1:2, mean, na.rm=TRUE),
                apply(detail, 1:2, var, na.rm=TRUE))
  cperm <- c(1, 2,6, 3,7, 4,8, 5,9)
  print(round(grid[,cperm], 4))
}

studyASEM <- function(modelGen, reps = 100, verbose=TRUE) {
  targets=c(seq(-8.1, -3.9, .2), seq(-5.8, -4.4, .2))
  # should not see any order effects, but just to check
  targets <- targets[order(runif(length(targets)))]
  rec <- c('condnum', 'asym', 'hNorm', 'rdNorm')
  detail <- array(NA, dim=c(length(targets), length(rec), reps),
                  dimnames=list(targets, rec, NULL))
  
  for (rep in 1:reps) {
    set.seed(rep)
    model <- modelGen()
    em <- model$compute
    em$tolerance <- 1e-10   # this only affects the MLE, not the individual trials
    em$information <- "mr1991"
    fitfun <- c()
    if (is(em$mstep, "MxComputeSequence")) {
      fitfun <- sapply(em$mstep$steps, function(step) step$fitfunction)
    } else {
      fitfun <- em$mstep$fitfunction
    }
    fit <- NULL
    
    for (tx in 1:length(targets)) {
      if (is.null(fit) || inherits(fit, "try-error")) fit <- model
      em$infoArgs=list(fitfunction=fitfun, semMethod="agile", semDebug=TRUE,
                       noiseTarget=exp(targets[tx]), semFixSymmetry=TRUE)
      plan <- mxComputeSequence(list(
        em,
        mxComputeStandardError(),
        mxComputeReportDeriv()
      ))
      fit$compute <- plan
      fit <- try(mxRun(fit, silent=TRUE), silent=TRUE)
      if (inherits(fit, "try-error")) {
#        print(paste("error in agile"))
#        print(fit)
        next
      } else {
        detail[tx,,rep] <- summarizeAgile(fit)
      }
    }
    if (verbose) summarizeASEM(detail)
  }
  detail
}

checkSmoothness <- function(mkmodel, probePoints=50) {
  set.seed(which(mcMask)[1])  # any positive definite model
  model <- mkmodel()
  em <- model$compute
  fitfun <- c()
  if (is(em$mstep, "MxComputeSequence")) {
    fitfun <- sapply(em$mstep$steps, function(step) step$fitfunction)
  } else {
    fitfun <- em$mstep$fitfunction
  }
  
  em$information <- "mr1991"
  em$tolerance <- 1e-9
  em$infoArgs <- list(fitfunction='fitfunction', semDebug=TRUE,
                      semMethod=seq(.0005, .01, length.out=probePoints))
  model$compute <- em
  model <- mxRun(model, silent=TRUE)
  em <- model$compute

  phl <- em$debug$paramHistLen
  probeOffset <- em$debug$probeOffset
  semDiff <- em$debug$semDiff
  
  upper <- 20
  modelfit <- list()
  result <- data.frame()
  for (vx in 1:length(model$output$estimate)) {
    len <- phl[vx]
    offset <- probeOffset[1:len, vx]
    dd <- semDiff[1:(len-1), vx]
    mid <- offset[1:(len-1)] + diff(offset)/2
    mask <- abs(diff(offset)) < .01 & dd < upper
    df <- data.frame(mid=mid[mask], diff=dd[mask])
    m1 <- lm(diff ~ 0 + I(1/mid^2), data=df)
    modelfit[[vx]] <- m1
    df$model <- predict(m1)
    result <- rbind(result, cbind(vx=vx, vname=names(model$output$estimate)[vx], df))
  }
  list(result=result, fits=modelfit, modelfit=sapply(modelfit, function(m) summary(m)$r.squ))
}

if (0) {
  # the worst fitting
  ggplot(subset(smooth$result, vx %in% order(smooth$modelfit)[1:4])) +
    geom_point(aes(mid, diff), size=2) + geom_line(aes(mid, model), color="green") +
    facet_wrap(~vname) + labs(x="x midpoint")
}


# Slow version:
#
# if (length(intersect(methods, "oaks"))) {
#   require(numDeriv)
#   stopifnot(sum(fit$item$free) == length(fit$output[['estimate']]))
#   
#   ejacob <- jacobian(function(eitemNew) {
#     eitem <- fit$item$values
#     eitem[fit$item$free] <- eitemNew
#     pm <- fit
#     pm$expectation$EstepItem <- eitem
#     pm$compute <- mxComputeSequence(list(
#       mxComputeOnce('expectation', 'scores'),
#       mxComputeOnce('fitfunction', 'gradient'),
#       mxComputeReportDeriv()))
#     pm <- mxRun(pm, silent=TRUE)
#     grad <- pm$output$gradient
#     # print(grad)
#     grad
#   }, fit$item$values[fit$item$free],
#   method="simple", method.args=list(eps=1e-3)
#   #  method="Richardson", method.args=list(r=2, d=1e-3, eps=1e-3)
#   )
#   
#   ejacob <- (ejacob + t(ejacob))/2  # usually almost symmetric anyway
#   
#   estepH <- mxModel(fit,
#                     mxComputeSequence(steps=list(
#                       mxComputeOnce(em$expectation, 'scores'),
#                       mxComputeOnce(fitfun, 'information', "hessian"),
#                       mxComputeReportDeriv())))
#   estepH <- mxRun(estepH, silent=TRUE)
#   H <- estepH$output$hessian + ejacob
#   se <- sqrt(2*diag(solve(H)))
#   ev <- eigen(H, TRUE, only.values=TRUE)$val
#   
#   detail[,"oaks",rep] <- 
#     summarizeInfo1(condnum=max(ev)/min(ev), H=H, standardErrors=se, cputime=NA, method="oaks")
# }
