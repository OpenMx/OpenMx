library(OpenMx)

a <- factor(sample(c('a', 'b', 'c'), size=100, replace=T))
b <- factor(a, levels=c('a', 'b', 'c', 'd')) #create factor is unused level 'd'
ma <- mxFactor(a, levels=levels(a))
mb <- mxFactor(b, levels=levels(b))
#mb[mb %in% 'c'] <- 'd' #make 'c' the unused level instead of 'd'
ds <- data.frame(a=ma, b=mb)

omxCheckError(mxDataWLS(ds, compute=TRUE), "fake.data: variable 'b' has a zero frequency category 'd'.
Eliminate this level in your mxFactor() or combine categories in some other way.
Do not pass go. Do not collect $200.")

Bollen1 <- Bollen
Bollen1[1,'y1'] <- NA
omxCheckError(mxDataWLS(Bollen1[, 1:8], compute=TRUE),
              "fake.data: all continuous data with missingness (column 'y1') cannot be handled using the cumulants method. Use na.omit(yourDataFrame) to remove rows with missing values or use allContinuousMethod='marginals' or use maximum likelihood")

ds <- data.frame(a=a, b=ma)
omxCheckWarning(
  omxCheckError(mxDataWLS(ds, compute=TRUE),
                "fake.data: variable 'a' must be an ordered factor but is of type unordered factor"),
  "In data 'fake.data', column 'a' must be an ordered factor. Please use mxFactor()")
