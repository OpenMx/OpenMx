options(error = utils::recover)
library(OpenMx)
library(rpf)
library(ggplot2)
library(stringr)

qwidth <- 5
quad <- imxEqualIntervalQuadratureData(21, qwidth)

data.raw <- read.csv("/opt/OpenMx/models/failing/data/cai-2009.csv")
data.g1 <- data.raw[data.raw$G==1, 2:13]
data.g2 <- data.raw[data.raw$G==2, 2:17]
colnames(data.g2) <- str_replace(colnames(data.g2), 'I', 'J')
data <- as.data.frame(cbind(data.g1, data.g2))
for (col in colnames(data)) data[[col]] <- ordered(data[[col]])

correct.g1 <- matrix(c(0.99, 1.42, 1.77, 2.19,  1.38, 1.8, 2.16, 1.18, 1.8, 2.61, 1.02, 1.69,
                       0.65,1.26, 1.21, 0.85, 1.06, 0.81, 1.58, 1.56, 1.08, 1.24, 0.73, 1.39,
                       0.88, 0.08, -0.35, -0.98, 0.99, 0.21, -0.42, -1.24, 0.81, 0.06, -0.29, -1.14,
                       rep(0,12)), byrow=TRUE, nrow=4)
correct.g2 <- matrix(c(0.99, 1.42, 1.77, 2.19, 1.38, 1.8,  2.16, 1.18, 1.8,  2.61, 1.02, 1.69, 1.76, 1.26, 1.45, 1.9,
                       0.65, 1.26, 1.21, 0.85, 1.06, 0.81, 1.58, 1.56, 1.08, 1.24, 0.73, 1.39, 1.21, 1.25, 0.99, 0.85,
                       0.88, 0.08,-0.35,-0.98, 0.99, 0.21,-0.42,-1.24, 0.81, 0.06,-0.29,-1.14, 0.88, 0.2, -0.35,-1.09,
                       rep(0,16)), byrow=TRUE, nrow=4)

# 
numItems <- dim(data)[2]
numPersons <- dim(data)[1]

spec <- list()
spec[1:numItems] <- rpf.drm(dimensions = 2)

fit.g1 <- rpf.ot2000.chisq(spec[1:12], correct.g1, correct.g1!=0, data[1:12], quad)
fit.g2 <- rpf.ot2000.chisq(spec[13:numItems], correct.g2, correct.g2!=0, data[13:numItems], quad)
print(sum(sapply(fit.g1, function (f) f$statistic), sapply(fit.g2, function (f) f$statistic)))

design <- matrix(c(rep(c(1,2),4),
                   rep(c(1,3),4),
                   rep(c(1,4),4),
                   rep(c(1,5),4),
                   rep(c(1,6),4),
                   rep(c(1,7),4),
                   rep(c(1,8),4)), ncol=12+16)

ip.mat <- mxMatrix(name="ItemParam", nrow=4, ncol=numItems,
                   values=c(1.4,1,0,0),
                   free=c(TRUE,TRUE,TRUE,FALSE),
                   lbound=c(1e-5, 1e-5, -qwidth, 0),
                   ubound=c(10, 10, qwidth, 0))

for (ix in 1:12) {
  for (px in 1:3) {
    name <- paste(c('p',px,',',ix), collapse='')
    ip.mat@labels[px,ix] <- name
    ip.mat@labels[px,12+ix] <- name
  }
}

fit1 <- function() {  
  m1 <- mxModel(model="drm1", ip.mat,
              mxMatrix(name="ItemSpec", nrow=6, ncol=numItems,
                       values=sapply(spec, function(m) slot(m,'spec')),
                       free=FALSE, byrow=TRUE),
              mxMatrix(name="Design", nrow=dim(design)[1], ncol=numItems, values=design),
              mxData(observed=data, type="raw"),
              mxExpectationBA81(
                ItemSpec="ItemSpec",
                ItemParam="ItemParam",
                Design="Design",
            		quadrature=quad),
              mxFitFunctionBA81())

  if (1) {
    m1 <- mxOption(m1, "Analytic Gradients", 'yes')
    m1 <- mxOption(m1, "Verify level", '-1')
  } else {
    m1 <- mxOption(m1, "Analytic Gradients", 'no')
  }
  m1 <- mxOption(m1, "Function precision", '1.0E-5')
  m1 <- mxOption(m1, "Calculate Hessian", "No")
  m1 <- mxOption(m1, "Standard Errors", "No")
  m1 <- mxRun(m1)
  m1
}

calc.bias <- function (bank) {
  bias <- matrix(0, nrow=dim(correct)[1], ncol=dim(correct)[2])
  for (sx in 1:length(bank)) {
    bias <- bias + bank[[sx]]$param
  }
  bias <- (bias / length(bank)) - correct
  bias
}

rda <- "cai2009-fit.rda"
if(0) {
  fit <- fit1()
  save(fit, file=rda)
} else {
  load(rda)
}


openmx.itemparam <- fit@matrices$ItemParam@values
openmx.g1 <- openmx.itemparam[,1:12]
openmx.g2 <- openmx.itemparam[,13:28]

fit.g1 <- rpf.ot2000.chisq(spec[1:12], openmx.g1, openmx.g1!=0, data[1:12], quad)
fit.g2 <- rpf.ot2000.chisq(spec[13:numItems], openmx.g2, openmx.g2!=0, data[13:numItems], quad)
print(sum(sapply(fit.g1, function (f) f$statistic), sapply(fit.g2, function (f) f$statistic)))
# 881

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
