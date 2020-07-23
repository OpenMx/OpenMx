#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


#------------------------------------------------------------------------------
# Authors: Robert M. Kirkpatrick, Jiang Yang
# Date: 2019-02-15
# Filename: omxReadGRMBin.R
# Purpose: Read GCTA-format GRMs
#------------------------------------------------------------------------------

#This function is adapted from syntax written by Jiang Yang for the GCTA User Manual:
omxReadGRMBin <- function(prefix, AllN=FALSE, size=4, returnList=FALSE){
	sum_i = function(i){
		return(sum(1:i))
	}	
	BinFileName = paste(prefix,".grm.bin",sep="")
	NFileName = paste(prefix,".grm.N.bin",sep="")
	IDFileName = paste(prefix,".grm.id",sep="")
	id = read.table(IDFileName)
	n = dim(id)[1]
	BinFile = file(BinFileName, "rb")
	grm = readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
	NFile = file(NFileName, "rb")
	if(AllN==T){
		N = readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
	}
	else{
		N = readBin(NFile, n=1, what=numeric(0), size=size)
	}
	i = sapply(1:n, sum_i)
	if(returnList){
		return(list(diag=grm[i], off=grm[-i], id=id, N=N))
	}
	else{
		GRM <- matrix(0.0, nrow=n, ncol=n, dimnames=list(id$V1+id$V2, id$V1+id$V2))
		GRM[!lower.tri(GRM,diag=T)] <- grm[-i]
		GRM <- GRM + t(GRM)
		diag(GRM) <- grm[i]
		return(GRM)
	}
}

