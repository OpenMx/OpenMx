library(OpenMx)
library(rpf)
library(numDeriv)
#options(error = utils::recover)

unpackHession <- function(deriv, np) {
  hess <- matrix(NA, nrow=np, ncol=np)
  dx <- np+1
  for (hr in 1:np) {
    hess[1:hr,hr] <- hess[hr,1:hr] <- deriv[dx:(dx+hr-1)]
    dx <- dx + hr
  }
  hess
}

#myseed <- as.integer(runif(1) * 1e7)
#print(paste("set.seed =",myseed))
#set.seed(myseed)
set.seed(1)

numItems <- 3
items <- list()
items[[1]] <- rpf.drm(factors=2)
items[[2]] <- rpf.grm(outcomes=5, factors=2)
items[[3]] <- rpf.nrm(outcomes=4, factors=2,
                      T.a="random", T.c="random")

data <- rpf.sample(5000, items)

starting <- list(c(1.4, 1, 0, .1, .9),
                 c(1.4, 1, seq(2,-2, length.out=4)),
                 c(1.4,  1,  rep(0,6)))
starting.len <- max(vapply(starting, length, 0))

ip.mat <- mxMatrix(name="itemParam", nrow=starting.len, ncol=numItems,
                   values=0, free=FALSE)

for (sx in 1:length(starting)) {
  v <- starting[[sx]]
  ip.mat@values[1:length(v),sx] <- v
  ip.mat@free[1:length(v),sx] <- TRUE
}
starting.free <- ip.mat@free
starting.values <- ip.mat@values

m.mat <- mxMatrix(name="mean", nrow=1, ncol=2, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=2, ncol=2, values=diag(2), free=FALSE)
m2 <- mxModel(model="drm1", ip.mat, m.mat, cov.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(mean="mean", cov="cov",
                                ItemSpec=items,
                                ItemParam="itemParam",
                                EItemParam=starting.values,
                qwidth=5, qpoints=21),
              mxFitFunctionML())

if (0) {   # enable to generate answer file
  samples.per.item <- 10
  ans <- list()
  for (spi in 1:samples.per.item) {
    for (ii in 1:numItems) {
      m2@matrices$itemParam@values <- starting.values
      m2@matrices$itemParam@free[,] <- FALSE
      
      spoint <- rpf.rparam(items[[ii]])
      np <- length(spoint)
      m2@matrices$itemParam@values[1:np,ii] <- spoint
      
      deriv <- genD(function(param) {
        np <- length(param)
        m2@matrices$itemParam@values[1:np,ii] <- param
        lModel <- mxModel(m2,
                          mxComputeSequence(steps=list(
                            mxComputeOnce('expectation', context='EM'),
                            mxComputeOnce('fitfunction', fit=TRUE)
                          )))
        fit <- mxRun(lModel, silent=TRUE)
        fit@output$minimum
      }, spoint)

      print(c(ii, spoint))
      ans[[length(ans)+1]] <- c(ii, spoint, deriv$D)
    }
  }
  
  ans.len <- max(sapply(ans, length))
  ans.padded <- lapply(ans, function (elem) elem <- c(elem, rep(NA, ans.len-length(elem))))
  write.table(t(simplify2array(ans.padded)), file="data/dLL.csv", row.names=FALSE, col.names=FALSE)
}

if (1) {
  ans <- suppressWarnings(try(read.table("models/passing/data/dLL.csv"), silent=TRUE))
  if (is(ans, "try-error")) ans <- read.table("data/dLL.csv")
  sqerror <- c()

  m2 <- mxModel(m2,
                mxComputeSequence(steps=list(
                  mxComputeOnce('expectation', context='EM'),
                  mxComputeOnce('fitfunction', gradient=TRUE, hessian=TRUE)
                )))

  for (tx in 1:dim(ans)[1]) {
    m2@matrices$itemParam@values <- starting.values
    m2@matrices$itemParam@free[,] <- FALSE
    
    ii <- ans[tx,1]
    np <- rpf.numParam(items[[ii]])
    spoint <- ans[tx,2:(np+1)]
    m2@matrices$itemParam@values[1:np,ii] <- simplify2array(spoint)
    m2@matrices$itemParam@free[,ii] <- starting.free[,ii]
    
    m2 <- mxRun(m2, silent=TRUE)
    
    grad1 <- m2@output$gradient
    names(grad1) <- NULL
    hess <- m2@output$hessian
    
    if (0) {
      print(paste("Item", ii))
      print(grad1)
      print(deriv$D[1:np])
      cat("T.a=",deparse(T.a),"\n")
      cat("T.c=",deparse(T.c),"\n")
      cat("an=",deparse(hess),"\n")
      cat("emp=",deparse(emp.hess),"\n")
      print(round(hess - emp.hess, 2))
    }
    
    emp.grad <- simplify2array(ans[tx,(2+np):(1+2*np)])
    emp.hess <- unpackHession(simplify2array(ans[tx, -1:-(1+np)]), np)
#    omxCheckCloseEnough(emp.grad, grad1, 1e-4)
#    omxCheckCloseEnough(emp.hess, hess, 1)
    diff <- abs(c(emp.grad - grad1, emp.hess - hess))
    sqerror <- c(sqerror, diff^2)
  }
  rms <- sqrt(sum(sqerror) / length(sqerror))
  print(rms)
  omxCheckTrue(rms < 3.4)
}
