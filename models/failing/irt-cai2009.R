# This data is from an email:
#
# Date: Wed, 06 Feb 2013 19:49:24 -0800
# From: Li Cai <lcai@ucla.edu>
# To: Joshua N Pritikin <jpritikin@pobox.com>
# Subject: Re: how did you control item bias in Cai (2010, p. 592) ?

options(error = utils::recover)
library(OpenMx)
library(rpf)
library(ggplot2)
library(stringr)

data.raw <- read.csv("/opt/OpenMx/models/failing/data/cai-2009.csv")
data.g1 <- as.data.frame(data.raw[data.raw$G==1, 2:13])
data.g2 <- as.data.frame(data.raw[data.raw$G==2, 2:17])
if (0) {
  write.table(data.g1, "cai2009-g1.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(data.g2, "cai2009-g2.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)  
}
for (col in colnames(data.g1)) data.g1[[col]] <- ordered(data.g1[[col]])
for (col in colnames(data.g2)) data.g2[[col]] <- ordered(data.g2[[col]])

mk.model <- function(model.name, data, latent.free) {
#   name <- "g1"
#   data <- data.g1
#   latent.free <- TRUE
#   
  numItems <- dim(data)[2]
  numPersons <- dim(data)[1]
  spec <- list()
  spec[1:numItems] <- rpf.grm(factors = 2)
  
  dims <- (1 + numItems/4)
  design <- matrix(c(rep(1,numItems),
                     kronecker(2:dims,rep(1,4))), byrow=TRUE, ncol=numItems)
  
  ispec <- mxMatrix(name="ItemSpec", nrow=3, ncol=numItems,
           values=sapply(spec, function(m) slot(m,'spec')),
           free=FALSE)

  ip.mat <- mxMatrix(name="ItemParam", nrow=3, ncol=numItems,
                     values=c(1.4,1,0),
                     free=c(TRUE,TRUE,TRUE),
                     lbound=c(1e-5, 1e-5, NA))
  ip.mat@free.group <- 'param'
  
  for (ix in 1:numItems) {
    for (px in 1:3) {
      name <- paste(c('p',ix,',',px), collapse='')
      ip.mat@labels[px,ix] <- name
    }
  }
  eip.mat <- mxAlgebra(ItemParam, name="EItemParam", fixed=TRUE)

  m.mat <- mxMatrix(name="mean", nrow=1, ncol=dims, values=0, free=latent.free)
  cov.mat <- mxMatrix(name="cov", nrow=dims, ncol=dims, values=diag(dims),
                      free=latent.free)
  
  m1 <- mxModel(model=model.name, ip.mat, eip.mat, ispec, m.mat, cov.mat,
                mxMatrix(name="Design", nrow=dim(design)[1], ncol=numItems, values=design),
                mxData(observed=data, type="raw"),
                mxExpectationBA81(
                  ItemSpec="ItemSpec",
                  Design="Design",
                  EItemParam="EItemParam",
                  mean="mean", cov="cov",
                  qpoints=21, qwidth=5, scores="full"),
                mxFitFunctionBA81(ItemParam="ItemParam", rescale=FALSE))
  m1
}

g1 <- mk.model("g1", data.g1, TRUE)
g2 <- mk.model("g2", data.g2, FALSE)

groups <- paste("g", 1:2, sep="")

if (1) {
  fm <- read.flexmirt("~/irt/cai2009/cai2009-prm.txt")  # use relative path TODO
  cg1 <- g1
  cg1@matrices$ItemParam@values <-
    rbind(fm$G1$param[1,], apply(fm$G1$param[2:4,], 2, sum), fm$G1$param[5,])
  cg1@matrices$mean@values <- t(fm$G1$mean)
  cg1@matrices$cov@values <- fm$G1$cov
  cg2 <- g2
  cg2@matrices$ItemParam@values <-
    rbind(fm$G2$param[1,], apply(fm$G2$param[2:5,], 2, sum), fm$G2$param[6,])
  cModel <- mxModel(model="cModel", cg1, cg2,
                    mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
                    mxComputeSequence(steps=list(
                      mxComputeOnce(paste(groups, 'expectation', sep="."), context='M'),
                      mxComputeOnce('fitfunction'))))
  cModel <- mxRun(cModel)
  # latent distribution is totally wrong:
#  cModel@submodels$g1@matrices$mean@values - fm$G1$mean
  print("correct:")
  print(fm$G1$mean)
  print(fm$G1$cov)
#  cModel@submodels$g1@matrices$cov@values - fm$G1$cov
#  cModel@output$minimum
}

if(0) {
grpModel <- mxModel(model="groupModel", g1, g2,
                    mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
                    mxComputeIterate(steps=list(
                      mxComputeOnce(paste(groups, "EItemParam", sep=".")),
                      mxComputeOnce(paste(groups, 'expectation', sep='.'), context='E'),
#                      mxComputeGradientDescent(free.group='param'),
                      mxComputeNewtonRaphson(free.group='param'),
                      mxComputeOnce(paste(groups, 'expectation', sep="."), context='M'),
                      mxComputeOnce('fitfunction')
                      ), verbose=TRUE))

#grpModel <- mxOption(grpModel, "Number of Threads", 1)

# NPSOL options:
grpModel <- mxOption(grpModel, "Analytic Gradients", 'Yes')
grpModel <- mxOption(grpModel, "Verify level", '-1')
grpModel <- mxOption(grpModel, "Function precision", '1.0E-7')

grpModel <- mxRun(grpModel)

calc.bias <- function (bank) {
  bias <- matrix(0, nrow=dim(correct)[1], ncol=dim(correct)[2])
  for (sx in 1:length(bank)) {
    bias <- bias + bank[[sx]]$param
  }
  bias <- (bias / length(bank)) - correct
  bias
}

#abs(calc.bias(bank[filter])[ip.mat@free])

if(0) {
  bygh <- bank[sapply(bank, function (b) b$seed == 1)]
  df <- as.data.frame(t(sapply(bygh, function(b) c(points=b$ghp, LL=b$LL))))
  ggplot(df, aes(x=points, y=LL)) + geom_line()
}

#qplot(c(0, 3), stat="function", fun=function (x) dlnorm(x, sdlog=.25), geom="line") + ylim(0,1.5)
if (0) {
  var.type <- apply(ip.mat@free,1,sum)
  bias <- calc.bias(bank[sapply(bank, function (b) b$ghp == 40)])
  df <- data.frame(x=t(correct)[t(ip.mat@free)],
                   var=c(rep('a', var.type[1]),
                         rep('b', var.type[2])),
                   bias=t(bias)[t(ip.mat@free)])
  df$var <- factor(df$var)
  ggplot(df, aes(x, bias, color=var)) + geom_point() + xlab("true parameter value")
}
}
