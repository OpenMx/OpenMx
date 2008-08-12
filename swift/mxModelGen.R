# generates all models for given number of connections
# and writes them out as mxmodel objects

library(OpenMx)
dyn.load("permall.so")

## should be like R_SWIFT_ARGS=3 3
allinputs   <- Sys.getenv("R_SWIFT_ARGS")
connections <- noquote(strsplit(allinputs," ")[[1]][1])
numcol      <- noquote(strsplit(allinputs," ")[[1]][2])
# need to set default here .75 if not given by user
#initweight    <- noquote(strsplit(allinputs," ")[[1]][3])
initweight <- 0.75
rowcol = as.integer(numcol)
size = rowcol*rowcol
retval = .C("run_perm", connections=as.integer(connections), size=as.integer(size), ncol=as.integer(rowcol), total=as.integer(0))

total_perms = retval$total

# read in the data files to generate models and write as objects
count = 0
fnums = array(c(0:total_perms))
for(i in fnums)
{
print("i is :")
print(i)
	model = MxModel()
	model$A <- FullMatrix(rowcol, rowcol, free=TRUE)
	model$S <- DiagMatrix(rowcol, rowcol, free=TRUE) 
	model$F <- IdenMatrix(rowcol, rowcol) # modifiable in future, for now no filter
	model$I <- IdenMatrix(rowcol, rowcol)
	model$cov <- MxAlgebra(model$F %&% (solve(model$I - model$A) %&% model$S))

	xindex = array(c(1:rowcol))
	yindex = array(c(1:rowcol))

	abits = as.matrix(read.table(paste(i,".adat",sep="")))
	sbits = as.matrix(read.table(paste(i,".sdat",sep="")))

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
	save.Object(model,file=paste(count,"_model.rdata", sep=""))
	count = count + 1
}
