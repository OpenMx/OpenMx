
source("funcKalman2.R")

require(mvtnorm)
require(Matrix)

xdim <- 3
udim <- 2
ydim <- 9
tdim <- 500
set.seed(948)
tA <- matrix(c(-.4, 0, 0, 0, -.9, .1, 0, -.1, -.9), xdim, xdim)
tB <- matrix(c(0), xdim, udim)
tC <- matrix(c(runif(3, .4, 1), rep(0, ydim), runif(3, .4, 1), rep(0, ydim), runif(3, .4, 1)), ydim, xdim)
tD <- matrix(c(0), ydim, udim)
tQ <- matrix(c(0), xdim, xdim); diag(tQ) <- runif(xdim, 1, 2) # 0, 1
tR <- matrix(c(0), ydim, ydim); diag(tR) <- runif(ydim)

x0 <- matrix(c(rnorm(xdim)), xdim, 1)
P0 <- diag(c(runif(xdim)))
tdx <- matrix(0, xdim, tdim+1)
tx <- matrix(0, xdim, tdim+1)
tu <- matrix(0, udim, tdim)
ty <- matrix(0, ydim, tdim)

tT <- matrix(0:tdim, nrow=1, ncol=tdim+1)#/100

tI <- diag(1, nrow=xdim)

tx[,1] <- x0
for(i in 2:(tdim+1)){
	q <- t(rmvnorm(1, rep(0, xdim), tQ))
	tdx[,i] <- tA %*% tx[,i-1] + tB %*% tu[,i-1] + q
	expA <- as.matrix(expm(tA * (tT[,i]-tT[,i-1])))
	intA <- solve(tA) %*% (expA - tI)
	tx[,i] <- expA %*% tx[, i-1] + intA %*% tB %*% tu[,i-1] + intA %*% q
	#ty[,i-1] <- tC %*% tx[,i-1] + tD %*% tu[,i-1] + t(rmvnorm(1, rep(0, ydim), tR))
	ty[,i-1] <- tC %*% tx[,i] + tD %*% tu[,i-1] + t(rmvnorm(1, rep(0, ydim), tR))
}

#plot(tx[1,], type='l')

rownames(ty) <- paste('y', 1:ydim, sep='')
rownames(tx) <- paste('x', 1:xdim, sep='')


#------------------------------------------------------------------------------
# Source most of funcKalman2.R

tlist2 <- list(A=tA, B=tB, C=tC, D=tD, Q=tQ, R=tR, x=x0, P=P0, Y=ty, U=tu, t=tT)
res2 <- KalmanPED.Hybrid(tlist2, predict.method="rk4")
res2$m2ll

res2a <- KalmanPED.Hybrid(tlist2, predict.method="expm")
res2a$m2ll

sqrt(mean((res2$X - res2a$X)^2))

mbuildfun <- function(x){
	mA <- matrix(c(x[1], 0, 0, 0, x[2], x[3], 0, -x[3], x[2]), xdim, xdim)
	mC <- matrix(c(x[7:9], rep(0, ydim), x[10:12], rep(0, ydim), x[13:15]), ydim, xdim)
	mR <- diag(x[16:24])
	return( list(A=mA, B=tB, C=mC, D=tD, Q=tQ, R=mR, x=x0, P=P0, Y=ty, U=tu, t=tT) )
}


tinit <- c(-.4, -.9, .1, diag(tQ), tC[tC!=0], diag(tR))
# mbuildfun(tinit)

dlmBegin <- Sys.time()
optfit <- optim(
	par=tinit,
	fn=KalmanOptim.Hybrid,
	predict.method="expm",
	buildfun=mbuildfun,
	method="L-BFGS-B",
	lower=c(rep(NA, 3), rep(0.00001, 3), rep(0.00001, 9), rep(0.00001, 9)),
	upper=c(rep(1.5, 3), rep(0.001, 3), rep(5, 9), rep(5, 9)),
	#method="BFGS",
	control=list(maxit=50)
)
dlmEnd <- Sys.time()
dlmEnd-dlmBegin

optfit


mbuildfun(optfit$par)$A; tA
mbuildfun(optfit$par)$C[tC!=0]; tC[tC!=0]
diag(mbuildfun(optfit$par)$R); diag(tR)


#res3 <- KalmanPED.Hybrid(mbuildfun(tinit))
#res3 <- res2
res3 <- KalmanPED.Hybrid(mbuildfun(optfit$par))


plot(res3$t[-1], tx[1,-1], type='l')
lines(res3$t[-1], res3$X[1,], lty=2, col='red')
cor(res3$X[1,], tx[1,-1]) #Not well correlated. Why?  N.B. Much better correlated when using underlying continuous model to generate data.
cor(res3$X[2,], tx[2,-1])
cor(res3$X[3,], tx[3,-1])



rms <- function(x, y){sqrt(mean((x-y)^2))}
rms(res3$X[1,], tx[1,-1])
rms(res3$X[2,], tx[2,-1])
rms(res3$X[3,], tx[3,-1])


plot(res3$t[-1], tx[2,-1], type='l')
lines(res3$t[-1], res3$X[2,], lty=2, col='red')


plot(res3$t[-1], tx[3,-1], type='l')
lines(res3$t[-1], res3$X[3,], lty=2, col='red')



ccf(tx[1,-1], res3$X[1,])
ccf(tx[2,-1], res3$X[2,])
ccf(tx[3,-1], res3$X[3,])

