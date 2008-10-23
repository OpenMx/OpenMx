# given the permutation number (model number) desired
# model will be generated and run

library(OpenMx)
dyn.load("/scratch/projects/tg/SIDGrid/usr/bin/singleperm.so")

## should be like R_SWIFT_ARGS="matrices/gestspeech.cov 4 6298"
allinputs   <- Sys.getenv("R_SWIFT_ARGS")
#inputfile <- Sys.getenv("R_INPUT")

print(allinputs)

covfile <- noquote(strsplit(allinputs," ")[[1]][1])
numcol      <- noquote(strsplit(allinputs," ")[[1]][2])
permnum      <- noquote(strsplit(allinputs," ")[[1]][3])

# need to set default here .75 if not given by user
#initweight    <- noquote(strsplit(allinputs," ")[[1]][3])
initweight <- 0.75
rowcol = as.integer(numcol)
size = rowcol*rowcol

print(permnum)

retval = .C("run_perm", size=as.integer(size), ncol=as.integer(rowcol), total=as.integer(0), permnum=as.integer(permnum))

total_perms = retval$total

covMatrix = as.matrix(read.table(covfile))

# read in the data files to generate models and write as objects
	model = MxModel()
	model$A <- FullMatrix(rowcol, rowcol, free=TRUE)
	model$S <- DiagMatrix(rowcol, rowcol, free=TRUE) 
	model$F <- IdenMatrix(rowcol, rowcol) # modifiable in future, for now no filter
	model$I <- IdenMatrix(rowcol, rowcol)
	model$cov <- MxAlgebra(model$F %&% (solve(model$I - model$A) %&% model$S))

	xindex = array(c(1:rowcol))
	yindex = array(c(1:rowcol))

	abits = as.matrix(read.table(paste(permnum,".adat",sep="")))
	sbits = as.matrix(read.table(paste(permnum,".sdat",sep="")))

	parcnt = 1.0
	scnt = rowcol+1.0

	for(x in xindex)
	{
	 for(y in yindex)
	 {
			if (sbits[x,y] == 1){
	        	model$S$parameters[x,y] = parcnt
			model$S$values[x,y] = initweight
			parcnt = parcnt + 1
			}
		else if (abits[x,y] == 1){
			model$A$parameters[x,y] = parcnt
			model$A$values[x,y] = initweight
			parcnt = parcnt + 1
			}
		}
	}

# now run simplecovariance 


objective <- CovarianceObjective(model$cov, covMatrix)
job <- MxJob(model, objective)
jobClosure <- createMxClosure(job, use_R=TRUE)

try((fmin=jobClosure()), FALSE)

if (exists("fmin")){
write.table(fmin$minimum, file=paste("results/",permnum,".min",sep=""), row.names=FALSE, col.names=FALSE)
}
if (!exists("fminobj")){
system(paste("touch results/",permnum,".min",sep=""))
}

save.Object(model,file=paste("results/",permnum,".rdata",sep=""))
