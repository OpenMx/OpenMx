#options(digits=20)
library(OpenMx)
library(rpf)
library(numDeriv)
#options(error = utils::recover)

mxOption(NULL, 'loglikelihoodScale', -1)

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

items <- list()
items[[1]] <- rpf.drm(factors=2)
items[[2]] <- rpf.grm(outcomes=5, factors=2)
items[[3]] <- rpf.nrm(outcomes=4, factors=2,
                      T.a=diag(3), T.c=diag(3))
numItems <- length(items)

params <- lapply(items, rpf.rparam, version=1)
data <- rpf.sample(5000, items, params)

starting <- list(c(1.4, 1, 0, logit(.1), logit(.9)),
                 c(1.4, 1, seq(2,-2, length.out=4)),
                 c(1.4,  1,  rep(0,6)))
starting.len <- max(vapply(starting, length, 0))

ip.mat <- mxMatrix(name="item", nrow=starting.len, ncol=numItems,
                   values=0, free=FALSE)
colnames(ip.mat) <- colnames(data)
rownames(ip.mat) <- c(paste('f', 1:2, sep=""), rep('p', nrow(ip.mat)-2))

for (sx in 1:length(starting)) {
  v <- starting[[sx]]
  ip.mat$values[1:length(v),sx] <- v
  ip.mat$free[1:length(v),sx] <- TRUE
}
starting.free <- ip.mat$free
starting.values <- ip.mat$values

m2 <- mxModel(model="drm1", ip.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(ItemSpec=items,
                                EstepItem=starting.values,
                qwidth=5, qpoints=21),
              mxFitFunctionML())

if (0) {   # enable to generate answer file
  samples.per.item <- 10
  ans <- list()
  for (ii in 1:numItems) {  # 
    spi <- 0
    while (spi < samples.per.item) {
      m2$item$values <- starting.values
      m2$item$free[,] <- FALSE
      
      spoint <- rpf.rparam(items[[ii]], version=1)
      # exclude GRM with close adjacent intercepts, too much curvature for numDeriv
      if (ii==2 && min(abs(diff(spoint[3:6]/max(spoint[1:2])))) < .3) next

      np <- length(spoint)
      m2$item$values[1:np,ii] <- spoint
      
      deriv <- genD(function(param) {
        np <- length(param)
        m2$item$values[1:np,ii] <- param
        lModel <- mxModel(m2,
                          mxComputeSequence(steps=list(
                            mxComputeOnce('expectation', 'scores'),
                            mxComputeOnce('fitfunction', 'fit')
                          )))
        fit <- mxRun(lModel, silent=TRUE)
        fit$output$fit
      }, spoint, method.args=list(d=.01, r=2))
      
      if (any(is.na(deriv$D))) next

      print(c(ii, spoint))
      ans[[length(ans)+1]] <- c(ii, spoint, deriv$D)
      spi <- spi + 1
    }
  }
  
  ans.len <- max(sapply(ans, length))
  ans.padded <- lapply(ans, function (elem) elem <- c(elem, rep(NA, ans.len-length(elem))))
  write.table(t(simplify2array(ans.padded)), file="data/dLL.csv", row.names=FALSE, col.names=FALSE)
}

ans <- suppressWarnings(try(read.table("models/nightly/data/dLL.csv"), silent=TRUE))
if (is(ans, "try-error")) ans <- read.table("data/dLL.csv")

m2 <- mxModel(m2,
              mxComputeSequence(steps=list(
                mxComputeOnce('expectation', 'scores'),
                mxComputeOnce('fitfunction', c('gradient', 'hessian')),
                mxComputeReportDeriv()
              )))

if (1) {  # enable to examine the RMSE by item model
  for (ix in 1:numItems) {
    np <- rpf.numParam(items[[ix]])
    diff.grad <- rep(0, np)
    diff.hess <- matrix(0, np, np)
    diff.count <- 0
    
    for (tx in 1:dim(ans)[1]) {
      ii <- ans[tx,1]
      if (ii != ix) next
      
      m2$item$values <- starting.values
      m2$item$free[,] <- FALSE
      
      spoint <- ans[tx,2:(np+1)]
      m2$item$values[1:np,ii] <- simplify2array(spoint)
      m2$item$free[,ii] <- starting.free[,ii]
      
      m2 <- mxRun(m2, silent=TRUE)
      
      grad1 <- m2$output$gradient
      names(grad1) <- NULL
      hess <- m2$output$hessian
      
      emp.grad <- simplify2array(ans[tx,(2+np):(1+2*np)])
      emp.hess <- unpackHession(simplify2array(ans[tx, -1:-(1+np)]), np)
      
      diff.grad <- diff.grad + (emp.grad - grad1)^2
      if (0) {
        print(tx)
        print(grad1)
        print(emp.grad)
      }
      diff.hess <- diff.hess + (emp.hess - hess)^2
      if (0) {
        print(tx)
        print(hess)
        print(emp.hess)
      }
      diff.count <- diff.count+1
    }

    if (diff.count > 0) {
      diff.grad <- sqrt(diff.grad / diff.count)
      diff.hess <- sqrt(diff.hess / diff.count)
      if (0) {
        print(round(diff.grad,3))
        print(max(diff.grad))
        print(round(diff.hess,3))
        print(max(diff.hess))
      }
    }
#     print(max(diff.grad))
#     print(max(diff.hess))
    # The poor accuracy here is probably due to numDeriv, not the
    # math for analytic derivs.
    omxCheckCloseEnough(diff.grad, rep(0,length(diff.grad)), 1e-5)
    omxCheckCloseEnough(c(diff.hess), rep(0,length(diff.hess)), 1e-3)
  }
}

if (0) {
  kat <- c()
  badest <-  c(6, 15, 40, 55, 67, 82, 84, 85, 87, 94)
 # badest <- c(107, 114, 121, 126, 132, 138, 139, 164, 172, 177, 182, 199)
      for (tx in badest) {
    ii <- ans[tx,1]
#    print(params[[ii]])
    np <- rpf.numParam(items[[ii]])
    
    m2$item$values <- starting.values
    m2$item$free[,] <- FALSE
    
    spoint <- simplify2array(ans[tx,2:(np+1)])
    print(spoint)
      next;
  
    m2$item$values[1:np,ii] <- spoint
    m2$item$free[,ii] <- starting.free[,ii]
      
    m2 <- mxRun(m2, silent=TRUE)
      
    grad1 <- m2$output$gradient
    names(grad1) <- NULL
    hess <- m2$output$hessian
    
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

#    print(grad1)
#    print(grad1 - emp.grad)
#    print(emp.hess)
#    print(hess)
#    print(hess - emp.hess)
    
    evalLL <- function(param) {
      np <- length(param)
      m2$item$values[1:np,ii] <- param
      lModel <- mxModel(m2,
                        mxComputeSequence(steps=list(
                          mxComputeOnce('expectation', 'scores'),
                          mxComputeOnce('fitfunction', 'fit')
                        )))
      fit <- mxRun(lModel, silent=TRUE)
      ll <- fit$output$minimum
      ll
    }
    deriv <- genD(evalLL, spoint, method.args=list(d=.1, r=3))
    print(round(hess - unpackHession(deriv$D, np), 2))
    
    if (0) {
      grid <- expand.grid(x=seq(.07,-.05,-.01))
      for (pp in grid$x) {
        spoint[5] <- pp
        grid$LL <- evalLL(spoint)
      }
    }
  }
}
