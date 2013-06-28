#options(error = utils::recover)
library(OpenMx)
library(rpf)
library(ggplot2)

numItems <- 20
numPersons <- 500

i1 <- rpf.drm(numChoices=4, a.prior.sdlog=.5, multidimensional=TRUE)
#i1 <- rpf.drm(numChoices=4, a.prior.sdlog=.5)
items <- list()
items[1:numItems] <- i1

if (0) {
  correct <- sapply(items, rpf.rparam)
  correct[2,] <- scale(correct[2,])
  correct[3,] <- 0
  cat(deparse(correct))
} else {
  correct <- structure(c(1.54250842425733, 1.29865606045677, 0, 0.582738899753006,  0.439074405155082, 0, 0.479997481961601, 1.6288972150727, 0,  1.69444870619279, 0.415803213429429, 0, 0.617838919096744, 0.0339637225500039,  0, 0.89104263140429, -0.310656987834084, 0, 0.856862341680302,  -1.39777061130952, 0, 0.399785394294146, 0.452539004701956, 0,  1.13750982818114, -1.12551515354692, 0, 0.551920087175635, 0.234962883698946,  0, 1.317228172309584, -0.767585189231164, 0, 1.06680521649664,  -0.17744848799824, 0, 0.820181809806619, 0.697978914209079, 0,  1.2446674744156, 0.579329622962913, 0, 3.17834154396065, 0.0422351712197786,  0, 1.29268733199003, -1.50421669038342, 0, 0.620813738032935,  -1.72835375813844, 0, 0.57191713624839, 0.748756461527565, 0,  1.44045131430693, 1.69557125401783, 0, 0.462792889649201, -0.956221050560265,  0), .Dim = c(3L, 20L), .Dimnames = list(c("a", "b", "c"), NULL))
}

ip.mat <- mxMatrix(name="ItemParam", nrow=3, ncol=numItems,
                   values=c(1,0,0),
                   free=c(TRUE,TRUE,FALSE),
                   lbound=c(1e-5, -1e6, 0))

fit1 <- function(seed=5, ghp=11) {
  result <- list(seed=seed, ghp=ghp, sdlog=i1@a.prior.sdlog, mdim=(class(i1) == "rpf.mdim.drm"))
  
  set.seed(seed)
  ability <- rnorm(numPersons)
  gen.param <- correct
  if (class(i1) == "rpf.mdim.drm") {
    gen.param[2,] <- gen.param[2,] * -gen.param[1,]
  }
  data <- rpf.sample(ability, items, gen.param)

  m1 <- mxModel(model="drm1", ip.mat,
              mxMatrix(name="ItemSpec", nrow=6, ncol=numItems,
                       values=sapply(items, function(m) slot(m,'spec')),
                       free=FALSE, byrow=TRUE),
              mxData(observed=data, type="raw"),
              mxExpectationBA81(
                ItemSpec="ItemSpec",
                ItemParam="ItemParam",
            		GHpoints=ghp),
              mxFitFunctionBA81())

  if (1) {
    result$gradient <- 1
    m1 <- mxOption(m1, "Analytic Gradients", 'yes')
    m1 <- mxOption(m1, "Verify level", '-1')
  } else {
    result$gradient <- 0
    m1 <- mxOption(m1, "Analytic Gradients", 'no')
  }
  m1 <- mxOption(m1, "Function precision", '1.0E-5')
  m1 <- mxOption(m1, "Calculate Hessian", "No")
  m1 <- mxOption(m1, "Standard Errors", "No")
  m1 <- mxRun(m1)

  result$cpuTime <- m1@output$cpuTime
  result$LL <- m1@output$Minus2LogLikelihood
  result$param <- m1@matrices$ItemParam@values
  if (class(i1) == "rpf.mdim.drm") {
    result$param[2,] <- result$param[2,] / -result$param[1,]
  }
  result
}

calc.bias <- function (bank) {
  bias <- matrix(0, nrow=dim(correct)[1], ncol=dim(correct)[2])
  for (sx in 1:length(bank)) {
    bias <- bias + bank[[sx]]$param
  }
  bias <- (bias / length(bank)) - correct
  bias
}

if (0) {
  bank <- list()
}
setwd("/opt/OpenMx")
rda <- "irt-drm1-bias.rda"
load(rda)
for (cnt in 1:500) {
  if (length(bank)) {
    filter <- sapply(bank, function (b) b$mdim==(class(i1) == "rpf.mdim.drm") & b$sdlog == .5 & b$ghp==17)
    if (any(sapply(bank, function (b) b$seed==cnt) & filter)) next;
  }
  bi <- length(bank) + 1
  bank[[bi]] <- fit1(ghp=17, seed=cnt)
  save(bank, file=rda)
  
  print(cor(bank[[bi]]$param[ip.mat@free], correct[ip.mat@free]))
#  print(cor(c(calc.bias(bank[filter])[ip.mat@free]), c(correct[ip.mat@free])))
}

for (cnt in 7:50) {
  filter <- sapply(bank, function (b) b$mdim==(class(i1) == "rpf.mdim.drm") & b$sdlog == .5)
  if (any(sapply(bank, function (b) b$ghp==cnt) & filter)) next;
  bi <- length(bank) + 1
  bank[[bi]] <- fit1(ghp=cnt, seed=1)
  save(bank, file=rda)
  
  print(cor(bank[[bi]]$param[ip.mat@free], correct[ip.mat@free]))
  #  print(cor(c(calc.bias(bank[filter])[ip.mat@free]), c(correct[ip.mat@free])))
}

#abs(calc.bias(bank[filter])[ip.mat@free])

if(1) {
  df <- list()
  for (mdim in c(TRUE,FALSE)) {
    x <- correct[ip.mat@free]
    bias <- calc.bias(bank[sapply(bank, function (b) b$mdim == mdim & b$ghp==17)])[ip.mat@free]
    df <- rbind(df, data.frame(mdim=mdim, bias=bias, x=x))
  }
  df$mdim <- factor(df$mdim)
  ggplot(df, aes(x=x, y=bias, color=mdim)) + geom_point()
}

if(1) {
  bygh <- bank[sapply(bank, function (b) b$seed == 1)]
  df <- as.data.frame(t(sapply(bygh, function(b) c(points=b$ghp, LL=b$LL, mdim=b$mdim))))
  df$mdim <- factor(df$mdim)
  ggplot(df, aes(x=points, y=LL, color=mdim)) + geom_line()
}

if(1) {
  bygh <- bank[sapply(bank, function (b) b$seed == 1)]
  df <- as.data.frame(t(sapply(bygh, function(b) c(mdim=b$mdim, points=b$ghp,
                                                   time=as.numeric(b$cpuTime)))))
  df$mdim <- factor(df$mdim)
  ggplot(df, aes(x=points, y=time, color=mdim)) + geom_line()
}

#qplot(c(0, 3), stat="function", fun=function (x) dlnorm(x, sdlog=.25), geom="line") + ylim(0,1.5)
if (0) {
  plot(correct[ip.mat@free],
       calc.bias(bank[sapply(bank, function (b) b$mdim==TRUE & b$sdlog==.5)])[ip.mat@free])
}
