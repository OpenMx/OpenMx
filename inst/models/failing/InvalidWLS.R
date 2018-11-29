library(OpenMx)

a <- factor(sample(c('a', 'b', 'c'), size=100, replace=T))
b <- factor(a, levels=c('a', 'b', 'c', 'd')) #create factor is unused level 'd'
ma <- mxFactor(a, levels=levels(a))
mb <- mxFactor(b, levels=levels(b))
#mb[mb %in% 'c'] <- 'd' #make 'c' the unused level instead of 'd'
ds <- data.frame(a=ma, b=mb)

omxCheckError(wd <- mxDataWLS(ds), "Variable 'b' has a zero frequency category 'd'.
Eliminate this level in your mxFactor() or combine categories in some other way.
Do not pass go. Do not collect $200.")

Bollen1 <- Bollen
Bollen1[1,'y1'] <- NA
omxCheckError(mxDataWLS(Bollen1[, 1:8]),
              "All continuous data with missingness cannot be handled in the WLS framework. Use na.omit(yourDa...")

badJointData <- jointData
badJointData[,c(2,4,5)] <- data.frame(mapply(factor, jointData[,c(2,4,5)],
                                            levels=list(c(0,1), c(0, 1, 2, 3), c(0, 1, 2)),
                 SIMPLIFY=FALSE), check.names = FALSE, row.names=rownames(jointData))
omxCheckError(mxDataWLS(badJointData),
              "Factors 'z2', 'z4', and 'z5' must be ordered and are not")
