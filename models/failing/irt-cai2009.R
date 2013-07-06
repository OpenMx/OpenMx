options(error = utils::recover)
library(OpenMx)
library(rpf)
library(ggplot2)
library(stringr)

data.raw <- read.csv("/opt/OpenMx/models/failing/data/cai-2009.csv")
data.g1 <- as.data.frame(data.raw[data.raw$G==1, 2:13])
data.g2 <- as.data.frame(data.raw[data.raw$G==2, 2:17])
for (col in colnames(data.g1)) data.g1[[col]] <- ordered(data.g1[[col]])
for (col in colnames(data.g2)) data.g2[[col]] <- ordered(data.g2[[col]])

correct.g1 <- matrix(c(0.99, 1.42, 1.77, 2.19,  1.38, 1.8, 2.16, 1.18, 1.8, 2.61, 1.02, 1.69,
                       0.65,1.26, 1.21, 0.85, 1.06, 0.81, 1.58, 1.56, 1.08, 1.24, 0.73, 1.39,
                       0.88, 0.08, -0.35, -0.98, 0.99, 0.21, -0.42, -1.24, 0.81, 0.06, -0.29, -1.14,
                       rep(0,12)), byrow=TRUE, nrow=4)
correct.g2 <- matrix(c(0.99, 1.42, 1.77, 2.19, 1.38, 1.8,  2.16, 1.18, 1.8,  2.61, 1.02, 1.69, 1.76, 1.26, 1.45, 1.9,
                       0.65, 1.26, 1.21, 0.85, 1.06, 0.81, 1.58, 1.56, 1.08, 1.24, 0.73, 1.39, 1.21, 1.25, 0.99, 0.85,
                       0.88, 0.08,-0.35,-0.98, 0.99, 0.21,-0.42,-1.24, 0.81, 0.06,-0.29,-1.14, 0.88, 0.2, -0.35,-1.09,
                       rep(0,16)), byrow=TRUE, nrow=4)

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
  cov.mat <- mxMatrix(name="cov", nrow=dims, ncol=dims, values=diag(dims), free=latent.free)
  
  m1 <- mxModel(model=model.name, ip.mat, eip.mat, ispec, m.mat, cov.mat,
                mxMatrix(name="Design", nrow=dim(design)[1], ncol=numItems, values=design),
                mxData(observed=data, type="raw"),
                mxExpectationBA81(
                  ItemSpec="ItemSpec",
                  Design="Design",
                  EItemParam="EItemParam",
                  mean="mean", cov="cov",
                  qpoints=21),
                mxFitFunctionBA81(ItemParam="ItemParam"))
  m1
}

g1 <- mk.model("g1", data.g1, TRUE)
g2 <- mk.model("g2", data.g2, FALSE)

groups <- paste("g", 1:2, sep="")
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

