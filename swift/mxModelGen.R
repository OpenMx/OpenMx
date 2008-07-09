# generates all models for given number of connections
# and writes them out as mxmodel objects

library(OpenMx)
dyn.load("/home/ac/skenny/perm.so")

## should be like R_SWIFT_ARGS=3 3

## get number of connections and array size (will be converted to nxn)

allinputs <- Sys.getenv("R_SWIFT_ARGS")

print(allinputs)

connections <- noquote(strsplit(allinputs," ")[[1]][1])

print(connections)

numcol        <- noquote(strsplit(allinputs," ")[[1]][2])

print(numcol)

rowcol = as.integer(numcol)
size = rowcol*rowcol


# creates the .adat files 

.C("run_perm", connections=as.integer(connections), size=as.integer(size))

# read files and generate associated ap.adat, s.dat and sp.dat files
# will do this in the c code

# read in the data files to generate models and write as objects

# dot is wildcard here, will need to adjust to it gets all 3 dat's
count = 0
for (i in list.files(pattern=".adat")){
	model = MxModel()
#	model$A <- FullMatrix(rowcol, rowcol, free=TRUE)
#	model$S <- DiagMatrix(rowcol, rowcol) 
	model$F <- IdenMatrix(rowcol, rowcol) # modifiable in future, for now no filter
	model$I <- IdenMatrix(rowcol, rowcol)
	model$cov <- MxAlgebra(model$F %&% (solve(model$I - model$A) %&% model$S))

# hard-coding for now
a = double(length = size)
apar = double(length = size)
xindex = array(c(1:rowcol))
yindex = array(c(1:rowcol))
	apar = matrix(c(read.table(i)), nrow = rowcol, ncol = rowcol)
	cnt = rowcol+1.0
for(x in xindex)
{
 for(y in yindex)
   {
	if (apar[x,y] == 1){
        apar[x,y] = cnt
	cnt = cnt + 1}

   }
}
print(apar)
aval = matrix(c(read.table(i)), nrow = rowcol, ncol = rowcol)

for(x in xindex)
{
 for(y in yindex)
   {
	if (aval[x,y] == 1)
        aval[x,y] = 0.75

   }
}
print(aval)

	
	s = matrix(c(read.table("matrices/s.dat")), nrow = rowcol, ncol = rowcol)
	spar = matrix(c(read.table("matrices/sp.dat")), nrow = rowcol, ncol = rowcol)
	
	model$A$values <- aval
	model$A$parameters <- apar
	model$S$values <- s
	model$S$parameters <- spar
	count = count+1
print(model$S)
	save.Object(model,file=paste(count,"_model.rdata", sep=""))
}
system("tar -cf models.tar *rdata")
