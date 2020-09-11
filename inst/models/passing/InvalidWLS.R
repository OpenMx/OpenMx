library(testthat)
library(OpenMx)

a <- factor(sample(c('a', 'b', 'c'), size=100, replace=T))
b <- factor(a, levels=c('a', 'b', 'c', 'd')) #create factor is unused level 'd'
ma <- mxFactor(a, levels=levels(a))
mb <- mxFactor(b, levels=levels(b))
#mb[mb %in% 'c'] <- 'd' #make 'c' the unused level instead of 'd'
ds <- data.frame(a=ma, b=mb)

expect_error(omxAugmentDataWithWLSSummary(mxData(ds, 'raw')),
"fake.data: variable 'b' has a zero frequency outcome 'd'.")

Bollen1 <- Bollen
Bollen1[1,'y1'] <- NA
omxCheckError(omxAugmentDataWithWLSSummary(mxData(Bollen1[, 1:8], 'raw')),
              "fake.data: all continuous data with missingness (column 'y1') cannot be handled using the cumulants method. Use na.omit(yourDataFrame) to remove rows with missing values or use allContinuousMethod='marginals' or use maximum likelihood")

ds <- data.frame(a=a, b=ma)
omxCheckError(omxAugmentDataWithWLSSummary(mxData(ds, 'raw')),
              "Don't know how to interpret unordered factor 'a' as numeric")
