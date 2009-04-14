# generates and runs a batch of models
# in parallel using swift.
# the workflow runs a check for independence
# before processing though in this example 
# all the models are generated as independent

library(OpenMx)
library(R.oo)
dyn.load("/scratch/projects/tg/SIDGrid/usr/bin/singleperm.so")

## set env with matrix info for generator
## should be like R_SWIFT_ARGS="matrices/gestspeech.cov 4 10"
allinputs   <- Sys.getenv("R_SWIFT_ARGS")

print(allinputs)

covfile <- noquote(strsplit(allinputs," ")[[1]][1])
numcol      <- noquote(strsplit(allinputs," ")[[1]][2])
max_permnum      <- noquote(strsplit(allinputs," ")[[1]][3])

# default initialization weight

initweight <- 0.75
rowcol = as.integer(numcol)
size = rowcol*rowcol
print(max_permnum)

covMatrix = as.matrix(read.table(covfile))

modnums = array(0:max_permnum)
### generation of models (not parallelized)

for(permnum in modnums)
{
	retval = .C("run_perm", size=as.integer(size), ncol=as.integer(rowcol), total=as.integer(0), permnum=as.integer(permnum))
	# read in the data files to generate models and write as objects
	model = mxModel()
	abits = as.matrix(read.table(paste(permnum,".adat",sep="")))
	sbits = as.matrix(read.table(paste(permnum,".sdat",sep="")))
	model <- mxModel(model, mxMatrix("Full", values = NA, name = "A", nrow = rowcol, ncol = rowcol))
	model <- mxModel(model, mxMatrix("Full", values = NA, name="S", nrow=rowcol, ncol=rowcol, free=TRUE))
	#hard-coded for now, will be identity matrix or other later
	model <- mxModel(model, mxMatrix("Full", values = NA, name="F", nrow=rowcol, ncol=rowcol))

	xindex = array(c(1:rowcol))
	yindex = array(c(1:rowcol))


	parcnt = 1.0
	scnt = rowcol+1.0

	for(x in xindex)
	{
	 for(y in yindex)
	 {
			if (sbits[x,y] == 1){
	        	model[["S"]]@spec[x,y] = parcnt
			model[["S"]]@values[x,y] = initweight
			parcnt = parcnt + 1
			}
		else if (abits[x,y] == 1){
			model[["A"]]@spec[x,y] = parcnt
			model[["A"]]@values[x,y] = initweight
			parcnt = parcnt + 1
			}
		}
}
	objective <- mxRAMObjective("objective")
	model <- mxModel(model, objective, covMatrix)
	model@independent = TRUE
	if(model@independent){
	save.Object(model,file=paste("models/",permnum,".rdata",sep=""))
	}
	retval = 0

}

# cleanup output from generator

system("rm *adat *sdat")

# call swift to run on models dir

mxGenSwift("./tc.data", "./sites.xml", "./optim.swift")



