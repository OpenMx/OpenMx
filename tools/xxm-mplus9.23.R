library(xxm)

ex923 <- read.table("ex9.23.dat")
colnames(ex923) <- c(paste0('y',1:4), 'x', 'w', 'z', 'level2', 'level3')
ex923$level2 <- as.integer(ex923$level2)
ex923$level3 <- as.integer(ex923$level3)

data1 <- cbind(data.frame(within=1:nrow(ex923)),
               ex923[,c('level2', paste0('y',1:4), 'x')])

data2 <- ex923[!duplicated(ex923$level2), c('level2', 'level3', 'w')]

data3 <- ex923[!duplicated(ex923$level3), c('level3', 'z')]

m1 <- xxmModel(levels = c("within", "level2", 'level3'))

l1 <- xxmSubmodel(model = m1, level = "within", parents = c("level2"), 
                      ys = paste0('y', 1:4),
                      xs = 'x',
                      etas = c("iw","sw"),
                      data = data1)

l2 <- xxmSubmodel(model = m1, level = "level2", parents = c('level3'),
                      ys = , xs = 'w',
                      etas = c("ib2",'sb2',paste0('y', 1:4)),
                  data = data2)

l3 <- xxmSubmodel(model = m1, level = "level3",
                  ys = , xs = 'z',
                  etas = c("ib3",'sb3'),
                  data = data3)

# ----------------------

# exo observed to latent
l1 <- xxmWithinMatrix(model = l1, level = "within", "gamma",
                          pattern=matrix(1, 2, 1),
                          value=matrix(.5, 2, 1))
l2 <- xxmWithinMatrix(model = l2, level = "level2", "gamma",
                      pattern=matrix(c(1,1,0,0,0,0), 6, 1),
                      value=matrix(0, 6, 1))
l3 <- xxmWithinMatrix(model = l3, level = "level3", "gamma",
                      pattern=matrix(c(1,1), 2, 1),
                      value=matrix(0, 2, 1))

# latent covariance
l1 <- xxmWithinMatrix(model = l1, level = "within", "psi",
                          pattern=matrix(1, 2, 2),
                          value=diag(2))

psiUpper2 <- diag(6)
psiUpper2[1:2,1:2] <- 1
l2 <- xxmWithinMatrix(model = l2, level = "level2", "psi",
                      pattern=psiUpper2,
                      value=diag(6))
l3 <- xxmWithinMatrix(model = l3, level = "level3", "psi",
                      pattern=matrix(1,2,2),
                      value=diag(c(1,1)))

# loadings
l1 <- xxmWithinMatrix(model = l1, level = "within", "lambda",
                      pattern = matrix(0,4,2),
                      value=matrix(c(1, 1, 1, 1, 0:3), 4, 2))

# loadings
l1 <- xxmWithinMatrix(model = l1, level = "within", "theta",
                      pattern = diag(4),
                      value=diag(4))

# upper level loadings
upperBeta <- matrix(0.0,6,6)
upperBeta[3:6,1] <- 1.0
upperBeta[3:6,2] <- 0:3
l2 <- xxmWithinMatrix(model = l2, level = "level2", "beta",
                      pattern = matrix(0,6,6),
                      value=upperBeta)

# latent intercepts
l3 <- xxmWithinMatrix(model = l3, level="level3", "alpha",
                          pattern=matrix(1,2,1),
                          value=matrix(0,2,1))

lambdaBet <- matrix(0,4,6)
lambdaBet[,3:6] <- diag(4)
m1 <- xxmBetweenMatrix(model = m1, parent = "level2", child = "within", 
                       type = "lambda",
                       pattern = matrix(0, 4, 6),
                       value = lambdaBet)
m1 <- xxmBetweenMatrix(model = m1, parent = "level3", child = "level2", 
                       type = "beta",
                       pattern = matrix(0, 6, 2),
                       value = matrix(c(0,0,rep(1,4),
                                        0,0,0:3), 6,2))

m1 <- xxmRun(m1)
