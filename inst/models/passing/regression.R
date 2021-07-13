library(OpenMx)
library(testthat)

set.seed(1)

N <- 10

x1 <- rnorm(N, mean=2, sd = 2)
x2 <- rnorm(N, mean=3, sd = 2)

df <- data.frame(y= 2 * x1 + 3 * x2 + rnorm(N, mean = 5, sd = 2),
                 x1=x1, x2=x2)

r1 <- lm(y~1, data = df)

m1 <- mxModel("lm", type="RAM",
              manifestVars = 'y',
              mxData(df, type="raw"),
              mxPath('y', arrows=2, values=5),
              mxPath('one', 'y'),
              mxFitFunctionWLS(allContinuousMethod = 'marginals'))
m1 <- mxRun(m1)
expect_equivalent(coef(r1), mean(df$y))
expect_equivalent(m1$data$observedStats$means, mean(df$y))
expect_equivalent(vcov(r1), m1$data$observedStats$y.vcov)

# ----------------------

r2 <- lm(y~x1+x2, data = df)

m2 <- mxModel("lm", type="RAM",
        manifestVars = 'y',
        latentVars = c('x1','x2'),
        mxData(df, type="raw"),
        mxPath('y', arrows=2, values=5),
        mxPath('one', 'y'),
        mxPath('one', c('x1','x2'), free=FALSE,
               labels=paste0("data.",c('x1','x2'))),
        mxPath(c('x1','x2'), 'y'))

m3 <- mxModel(m2,
              mxFitFunctionWLS(allContinuousMethod = 'marginals'))

m2 <- mxRun(m2)
m3 <- mxRun(m3)

expect_equivalent(m3$data$observedStats$y.vcov, vcov(r2))

expect_equivalent(coef(r2), coef(m2)[c(4,1:2)], tolerance=1e-4)
expect_equal(coef(m2), coef(m3), tolerance=1e-4)
