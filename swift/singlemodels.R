# generates matrix based on given perm number and size

# determine number of connections for given permutation

mxGetPerm <- function(rowcol,perm){
	perms = 0
	size = rowcol*rowcol
	print(paste("getting model for permutation ", perm, " of ",2^size))
	for(k in c(0:size))
		{
		perms = perms+choose(size,k)
		if(perm <= perms)
			{
			nconn = k
			break
			}
		}
	print(paste("number of connections: ",nconn))
	range = 0
	connslots = vector(length=nconn, mode="list")
	for(i in c(0:(nconn-1)))
		{
		range = range+choose(size,i)
		}

	perm = perm-range
	print(paste("secondary perm is ", perm))
	nslotranges = list()
	n = 1
	nslotranges[[n]] = vector(length=(size-(nconn-1)), mode="list")
	for(r in c(1:(size-(nconn-1))))
   		{
		      nslotranges[[n]][r] = choose((size-r),nconn-n)

        	}
	range = 0
	for(s in c(1:length(nslotranges[[n]])))
	{
		curr_range = as.numeric(nslotranges[[n]][s])
		range = range+curr_range
		if(perm <= range)
			{
			connslots[[n]] = s # changed to 1 indexing
			break
			}
	}
	rangegroup = connslots[[n]]
	if(rangegroup>1)
	{
	   for(y in 1:(rangegroup-1))
	   {
		perm = perm-as.numeric(nslotranges[[n]][y])
	   }
	}

# initialize matrix, determine new permutation
# index based on range for number of connections
# and first non-zero index within the matrix

idata = list()
istart = connslots[[1]]
istop = istart+(nconn-1)
for(x in 1:size)
	{
	if(x<istart)
	   idata[[x]] = 0
	if(x>=istart && x<=istop)
	   idata[[x]] = 1
	if(x>istop)
	   idata[[x]] = 0
	}
print(paste("perm is now ",perm))
print(paste("bit1 slot is ",connslots[[n]]))
perm_array = mxPerm(perm,nconn,idata)
return(perm_array)
}

# ------------------------------------------

# iterate thru matrices in subset 
# until matching index is found
# that is the matrix to be optimized

find_n <- function(n,idata){
	cnt = 0
	for(j in c(1:length(idata)))
	{
	   if(idata[[j]] == 1)
		cnt = cnt+1
	   if(cnt == n)
		return(j)
	}
	return(j)
}

find_shift_elmt <- function(connections,idata){
	x = connections-1
	while(x>=1)
	{
	   curr_n = find_n(x,idata)
	   if(idata[[curr_n+1]] == 0)
		return(curr_n)
	x = x-1
	}
	return(0)
}

insert_zero <- function(pos,connections,alldata) {
	size = length(alldata[[1]])

	if( pos+1>size || alldata[[1]][pos+1] != 0)
	{
	   selmt_pos = find_shift_elmt(connections, alldata[[1]])
    	   pos = selmt_pos;
	}
  	bitcount = 0;
        alldata[[2]] = alldata[[1]]
	i = size;
  	# set position of interest 
  	alldata[[1]][pos] = 0;

  # set those to left which should 
  # be static but keep track of them
if(pos>1){
	for(k in c(1:(pos-1)))
	{
   	alldata[[1]][k] = alldata[[2]][k];
	if(alldata[[1]][k] == 1){ bitcount = bitcount+1; }
	}
}
  # now remaining bits not set on left
  # are set immediately to the right 

	rbits = connections-bitcount
	rcnt = 0	
	newpos = pos+1

	while(rcnt<rbits)
	{
	alldata[[1]][newpos] = 1;
	rcnt = rcnt+1
	newpos = newpos+1
    	}
	if(newpos<=size)
	{
	   for(p in newpos:size)
	   {
           alldata[[1]][p] = 0;
	   }
	}
return(alldata)
}

mxPerm <- function(permnum,connections,idata) {
	recentdata = idata
	alldata = list()
	alldata[[1]] = idata
	alldata[[2]] = recentdata
	size = length(idata)
	rowcol = sqrt(size)
  	k = 0;
  	asymm = 0;
  	last_shift_pos = size - connections;
	permcount = 0

	while(k<=permnum)
	{
	   alldata = insert_zero(find_n(connections,alldata[[1]]), connections, alldata)
	   permcount = permcount+1
	   k = k+1
	}
return(alldata[[1]])
}

# --------------- main ------------------------

library(OpenMx)
library(R.oo)

# get parameters from environment variable R_SWIFT_ARGS: 
# this allows variables to be easily passed by swift and/or
# other R wrappers
# example R_SWIFT_ARGS:
# 	R_SWIFT_ARGS="net1_gesture.cov 57134 .5 gesture 1 8"

allinputs <- Sys.getenv("R_SWIFT_ARGS")
covfile <- noquote(strsplit(allinputs," ")[[1]][1])
perm <- as.numeric(noquote(strsplit(allinputs," ")[[1]][2]))
initweight <- as.numeric(noquote(strsplit(allinputs," ")[[1]][3]))
cond <- noquote(strsplit(allinputs," ")[[1]][4])
net <- noquote(strsplit(allinputs," ")[[1]][5])
num_obs <- as.numeric(noquote(strsplit(allinputs," ")[[1]][6]))
covMatrix = as.matrix(read.table(covfile))
rowcol = dim(covMatrix)[[1]]

# make sure index is in set

if(perm <= 2^(rowcol*rowcol))
{

permutation = mxGetPerm(rowcol,perm)
pmatrix = matrix(c(permutation),nrow=rowcol)

model <- mxModel(name = paste("perm",perm, sep=""), type = "RAM", independent = TRUE)
model <- mxModel(model, mxMatrix("Full", values = NA, name="A", nrow = rowcol, ncol = rowcol, lbound = NA, ubound = NA))
model <- mxModel(model, mxMatrix("Symm", values = NA, name="S", nrow=rowcol, ncol=rowcol, free=TRUE, lbound = NA, ubound = NA))
model <- mxModel(model, mxMatrix("Iden", values = NA, name="F", nrow=rowcol, ncol=rowcol, free=FALSE))


parcnt = 1
for(x in 1:rowcol)
{
  for(y in 1:rowcol)
   {
 # initialize 0 matrix
    model[["S"]]@free[x,y] = FALSE
    model[["S"]]@labels[x,y] = 0
    model[["S"]]@values[x,y] = 0
    model[["A"]]@free[x,y] = FALSE
    model[["A"]]@values[x,y] = 0
    model[["A"]]@labels[x,y] = 0

    if(pmatrix[x,y] == 1 && x == y){
    model[["S"]]@labels[x,y] = parcnt
    model[["S"]]@values[x,y] = initweight
    model[["S"]]@free[x,y] = TRUE
    parcnt = parcnt + 1
    }
    else if (pmatrix[x,y] == 1){
    model[["A"]]@labels[x,y] = parcnt
    model[["A"]]@values[x,y] = initweight
    model[["A"]]@free[x,y] = TRUE
    parcnt = parcnt + 1
    }

   }
}

model@manifestVars = "fMRI"

# run model

freeparams = c(model[["A"]]@labels,model[["S"]]@labels)
data <- mxData(covMatrix, 'cov', numObs = num_obs)
objective <- mxRAMObjective("A", "S", "F")
model <- mxModel(model, objective, data)
model <- mxRun(model)

# print degrees of freedom along with fit statistic
# for given model so results can be summarized

deg_of_freedom = (rowcol*(rowcol+1)/2)-parcnt

perm = formatC(perm, format="f", digits=0)
if(is.na(model@output[[1]]) || model@output[[1]] == -Inf){
		system(paste("touch results/",net,"_",cond,"_",perm,".stat", sep=""))
				     } else {
write(paste(deg_of_freedom,model@output[[1]],sep="\t"), file=paste("results/",net,"_",cond,"_",perm,".stat", sep=""))
					    }
	print(model@output)
}else {
	system(paste("touch results/",net,"_",cond,"_",perm,".stat", sep=""))
      }

