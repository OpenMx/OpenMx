# Maes, H. H., Neale, M. C., Kendler, K. S., Hewitt, J. K., Silberg, J. L., Foley, D. L., ... & Eaves, L. J. (1998). Assortative mating for major psychiatric diagnoses in two population-based samples. Psychological medicine, 28(6), 1389-1401.
#
# https://vipbg.vcu.edu/vipbg/Articles/psychological-assortative-1998.pdf

library(OpenMx)
library(testthat)

diagno <- toupper(c('alc', 'gad', 'mdd', 'pan', 'pho'))
parent <- c('wife', 'husband')

c1 <- apply(expand.grid(diagno, parent)[,c(2,1)], 1, paste, collapse='')
wife <- c1[1:5]
husband <- c1[6:10]

c2 <- matrix(0, length(c1), length(c1),
       dimnames=list(c1,c1))

# Table 2 ABD upper triangle
c2['wifeALC', 'wifeGAD'] <- .399
c2['wifeALC', 'wifeMDD'] <- .341
c2['wifeALC', 'wifePAN'] <- .276
c2['wifeALC', 'wifePHO'] <- .284
c2['wifeGAD', 'wifeMDD'] <- .582
c2['wifeGAD', 'wifePAN'] <- .372
c2['wifeGAD', 'wifePHO'] <- .394
c2['wifeMDD', 'wifePAN'] <- .31
c2['wifeMDD', 'wifePHO'] <- .304
c2['wifePAN', 'wifePHO'] <- .371

# Table 2 ABD lower triangle
c2['husbandGAD', 'husbandALC'] <- .295
c2['husbandMDD', 'husbandALC'] <- .341
c2['husbandMDD', 'husbandGAD'] <- .576
c2['husbandPAN', 'husbandALC'] <- .283
c2['husbandPAN', 'husbandGAD'] <- .314
c2['husbandPAN', 'husbandMDD'] <- .576
c2['husbandPHO', 'husbandALC'] <- .249
c2['husbandPHO', 'husbandGAD'] <- .220
c2['husbandPHO', 'husbandMDD'] <- .181
c2['husbandPHO', 'husbandPAN'] <- .276

# Table 2 ABD
c2[wife,wife] + c2[husband, husband]

# Table 3 ABD
c2['wifeALC', husband] <- c(.119, .064, .19, .029, .005)
c2['wifeGAD', husband] <- c(.209, .208, .207, .119, .133)
c2['wifeMDD', husband] <- c(.132, .088, .162, .112, .105)
c2['wifePAN', husband] <- c(.064, .092, .151, .219, .102)
c2['wifePHO', husband] <- c(-.079, .13, .221, -.016, .064)

# Repair symmetry
c2[husband,wife] <- t(c2[wife, husband])
c2[wife, wife][lower.tri(diag(length(diagno)))] <-
  t(c2[wife, wife])[lower.tri(diag(length(diagno)))]
c2[husband, husband][upper.tri(diag(length(diagno)))] <-
  t(c2[husband, husband])[upper.tri(diag(length(diagno)))]
diag(c2) <- 1

c3 <- matrix(0, length(c1), length(c1),
             dimnames=list(c1,c1))

# Table 2 AFT upper triangle
c3['wifeALC', 'wifeGAD'] <- .378
c3['wifeALC', 'wifeMDD'] <- .313
c3['wifeALC', 'wifePAN'] <- .076
c3['wifeALC', 'wifePHO'] <- .095
c3['wifeGAD', 'wifeMDD'] <- .709
c3['wifeGAD', 'wifePAN'] <- .471
c3['wifeGAD', 'wifePHO'] <- .228
c3['wifeMDD', 'wifePAN'] <- .367
c3['wifeMDD', 'wifePHO'] <- .227
c3['wifePAN', 'wifePHO'] <- .371

# Table 2 AFT lower triangle
c3['husbandGAD', 'husbandALC'] <- .125
c3['husbandMDD', 'husbandALC'] <- .185
c3['husbandMDD', 'husbandGAD'] <- .64
c3['husbandPAN', 'husbandALC'] <- .267
c3['husbandPAN', 'husbandGAD'] <- .483
c3['husbandPAN', 'husbandMDD'] <- .467
c3['husbandPHO', 'husbandALC'] <- .194
c3['husbandPHO', 'husbandGAD'] <- .312
c3['husbandPHO', 'husbandMDD'] <- .149
c3['husbandPHO', 'husbandPAN'] <- .344

# Table 2 AFT
c3[wife,wife] + c3[husband, husband]

# Table 3 AFT
c3['wifeALC', husband] <- c(.068, .213, .267, -.047, -.078)
c3['wifeGAD', husband] <- c(.107, -.014, .169, .116, -.006)
c3['wifeMDD', husband] <- c(.251, .014, .121, -.196, .074)
c3['wifePAN', husband] <- c(.002, .114, .015, .138, .062)
c3['wifePHO', husband] <- c(.003, .14, .033, .028, .071)

# Repair symmetry
c3[husband,wife] <- t(c3[wife, husband])
c3[wife, wife][lower.tri(diag(length(diagno)))] <-
  t(c3[wife, wife])[lower.tri(diag(length(diagno)))]
c3[husband, husband][upper.tri(diag(length(diagno)))] <-
  t(c3[husband, husband])[upper.tri(diag(length(diagno)))]
diag(c3) <- 1

diag2d <- matrix(apply(expand.grid(diagno, diagno)[,c(2,1)],
                       1, paste, collapse=""),
                 length(diagno), length(diagno))

paths <- c(
  mxPath(wife, arrows=2, connect="unique.bivariate",
         lbound=-.99, ubound=.99,
         labels=paste0('w', diag2d[lower.tri(diag2d)])),
  mxPath(husband, arrows=2, connect="unique.bivariate",
         lbound=-.99, ubound=.99,
         labels=paste0('h', diag2d[lower.tri(diag2d)])),
  mxPath(wife, husband, arrows=0, connect="unique.pairs",
         lbound=-.99, ubound=.99,
         labels=paste0('d', c(diag2d))))

abd <- mxModel(
  'abd', type="RAM",
  manifestVars = c1,
  mxData(c2, type = "cov", numObs = 854),
  mxPath(c1, arrows=2, values=1, free=FALSE), paths)

aft <- mxModel(
  'aft', type="RAM",
  manifestVars = c1,
  mxData(c3, type = "cov", numObs = 568),
  mxPath(c1, arrows=2, values=1, free=FALSE), paths)

fig6 <- mxModel("fig6", abd, aft,
                mxFitFunctionMultigroup(c("abd", "aft")))

if (0) {
  fig6 <- mxOption(fig6,"Always Checkpoint","Yes")
  fig6 <- mxOption(fig6,"Checkpoint Units","evaluations")
  fig6 <- mxOption(fig6,"Checkpoint Count",1)
  fig6 <- mxOption(fig6, "Checkpoint Fullpath", "/dev/fd/2")
}

fig6 <- mxRun(fig6)
fig6r <- mxRefModels(fig6, run=TRUE)
fig6s <- summary(fig6, refModels=fig6r)

expect_equal(fig6s$ChiDoF, 65)
expect_equal(fig6s$Chi, 544.0544, 1e-3)
expect_equal(fig6s$CFI, .8756, 1e-4)
