library(xxm)

ex96 <- read.table("ex9.6.dat")
ex96$V8 <- as.integer(ex96$V8)
bData <- ex96[!duplicated(ex96$V8), c('V7', 'V8')]
colnames(bData) <- c('w', 'clusterID')
bData <- bData[,c('clusterID', colnames(bData)[-2])]

wData <- ex96[,-match(c('V7'), colnames(ex96))]
colnames(wData) <- c(paste0('y', 1:4), paste0('x', 1:2), 'clusterID')
wData$within <- 1:nrow(wData)
wData <- wData[,c('within', 'clusterID', colnames(wData)[c(-7,-8)])]

m1 <- xxmModel(levels = c("within", "clusterID"))

wModel <- xxmSubmodel(model = m1, level = "within", parents = c("clusterID"), 
                     ys = paste0('y', 1:4),
                     xs = paste0('x', 1:2),
                     etas = "fw",
                     data = wData)

cModel <- xxmSubmodel(model = m1, level = "clusterID", parents = ,
                      ys = , xs = 'w', 
                     etas = c("fb"), data = bData)

# -------------------------- clusterID model

# observed to latent
cModel <- xxmWithinMatrix(model = cModel, level = "clusterID", "gamma",
                          pattern=matrix(1, 1, 1),
                          value=matrix(.5, 1, 1))

# latent covariance
cModel <- xxmWithinMatrix(model = cModel, level = "clusterID", "psi",
  pattern=diag(c(1)),
  value=diag(c(1)))

# -------------------------- within model
# loadings
wModel <- xxmWithinMatrix(model = wModel, level = "within", "lambda",
                          pattern = matrix(c(0, 1, 1, 1), 4, 1),
                          value = matrix(rep(1,4), 4, 1))

# residuals
wModel <- xxmWithinMatrix(model = wModel, level="within", "theta",
                          pattern=diag(4), value=diag(4))

# observed intercepts
wModel <- xxmWithinMatrix(model = wModel, level="within", "nu",
                          pattern=matrix(1,4,1), value=matrix(0,4,1))

# latent intercepts
wModel <- xxmWithinMatrix(model = wModel, level="within", "alpha",
                          pattern=matrix(c(0),1,1),
                          value=matrix(0,1,1))

# latent covariance
wModel <- xxmWithinMatrix(model = wModel, level="within", "psi",
                          pattern=diag(c(1)),
                          value=diag(c(1)))

# observed to latent
wModel <- xxmWithinMatrix(model = wModel, level="within", "gamma",
                          pattern=matrix(1, 1, 2),
                          value=matrix(.5, 1, 2))

m1 <- xxmBetweenMatrix(model = m1, parent = "clusterID", child = "within", 
                        type = "lambda",
                       pattern = matrix(c(0, 1, 1, 1), 4, 1),
                       value = matrix(rep(1,4), 4, 1))

m1 <- xxmRun(m1)
