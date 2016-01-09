#
#   Copyright 2007-2016 The OpenMx Project
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

######################################################################## 
#RUNMX
#R function by Matthew C Keller & Mike Neale
#Many thanks to Eric Schmitt, Steve Boker, Gary Xie, and Sara Medland
# 5/30/2007


##' imxOriginalMx
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param mx.filename mx.filename
##' @param output.directory output.directory
imxOriginalMx <- function(mx.filename, output.directory) {
	original.directory <- getwd()
	result <- tryCatch(originalMxHelper(mx.filename, output.directory),
		finally=setwd(original.directory))
	return(result)
}

originalMxHelper <- function(mx.filename, output.directory) {
  
	if (is.na(file.info(mx.filename)$size)) {
		stop(paste("Cannot find the file named", mx.filename))
	}
  
	############ 
	#Grab the mx.filename, place "!@front" at the top, and write it back out to the working directory
	mx.bottom <- suppressWarnings(readLines(mx.filename))
	paths <- strsplit(mx.filename, '[/\\]')[[1]]
	mx.filename <- paths[[length(paths)]]
	mx.tot <- c("!@front",mx.bottom)
	mxfile <- paste(output.directory, "/", mx.filename, sep="")
	write(mx.tot,file=mxfile)
	mxo.filename <- paste(gsub(".mx","",mx.filename),".mxo",sep="") 
	############ 

	setwd(output.directory)

	############ 
	#Run MX file
	if (.Platform$OS.type == "windows") {
		command <- substitute(system(paste("mx", x, y)), list(x = mx.filename, y = mxo.filename))
		eval(command)
	}
	
	if (.Platform$OS.type == "unix") {
		command <- substitute(system(paste("mx < ", x, "> ", y)), list(x = mx.filename, y = mxo.filename))
		eval(command)
	}
	############ 


	############ 
	#Read in the MXO file
	mxo <- suppressWarnings(readLines(mxo.filename))

	#Check for Mx Errors
	if (length(grep("Error: file not found",mxo)) != 0) {
	  stop(paste("It appears that Mx cannot find the datafile it needs.",
	        "Make sure the datafile(s) are placed in your working directory"))
	}

	if (length(grep("!@error",mxo)) != 0) {
	  stop(paste("Your Mx script did not complete running for some reason. Check over",
	        "your *.mxo file to diagnose the problem and try, try again."))
	}

	#Find the start of the fit functions
	if (length(grep("!@machine; GROUP_FIT",mxo)) != 0) {  #grepping for how windows seems to write the .mxo file
   		fits <- grep("!@machine; GROUP_FIT",mxo)
 	} else {
	   machines <- grep("!@machine;",mxo)                 #grepping for how unix seems to write the .mxo file
	   grpfits <- grep("GROUP_FIT",mxo)
	   tots <- c(machines,grpfits)
	   tots <- tots[order(tots)]
	   fits <- tots[c(FALSE,diff(tots)==1)]
	}

	#Find the start of the matrices
	if (length(grep("!@machine; SPECS",mxo)) != 0) {      #grepping for how windows seems to write the .mxo file
	   mats <- grep("!@machine; SPECS",mxo)
	} else {
	   machines <- grep("!@machine;",mxo)                 #grepping for how unix seems to write the .mxo file
	   grpspecs <- grep("SPECS",mxo)
	   tots <- c(machines,grpspecs)
	   tots <- tots[order(tots)]
	   mats <- tots[c(FALSE,diff(tots)==1)]
	}

############ 



############ 
#Number of jobs, starts & ends
	number.jobs <- length(grep("!@SUBJOB",mxo))+1
	x.starts <- mats
	x.ends <- fits[(1:number.jobs)*2]
	y.starts <- fits[((1:number.jobs)*2)-1]
	y.ends <- mats

	#Start the list
	matrices <- list()

	#Outer Loop - looping through different jobs
	for (k in 1:number.jobs){
		x <- mxo[x.starts[k]:x.ends[k]]  #x = matrices
 		y <- mxo[y.starts[k]:y.ends[k]]  #y = fit functions
		j <- 0
		ms <- grep("VALUE",x)
		me <- grep("^ ;",x)
		me <- me[seq(2,length(me),by=2)]

		#Remove the parts of "x" that break matrix elements into 100 element 
		#chunks & recreate "x" s.t. matrices are contiguous
		to.remove.ends <- (as.numeric(x[ms+4])*as.numeric(x[ms+5])) - as.numeric(x[ms+6]) > 100
		start.remove <- (ms+108)*to.remove.ends
		start.remove <- start.remove[start.remove>0]
		end.remove <- (ms+116)*to.remove.ends
		end.remove <- end.remove[end.remove>0]

		if (length(start.remove)>0) {
			removes <- vector()
  			for (u in 1:length(start.remove)) {
		    	removes <- c(removes,start.remove[u]:end.remove[u])
			}
			x <- x[-removes]
		}

		#Find starts and ends of new matrix "x" 
		mat.starts <- grep("VALUE",x)
 		mat.ends <- grep("^ ;",x)
 		mat.ends <- mat.ends[seq(2,length(mat.ends),by=2)]
		############ 

 
 
		############ 
		#Finding the fit function information if "User defined function value" is in the mxo file
  		if (length(grep("User defined function value =",y)) != 0) {
   			User.fit <- as.numeric(strsplit(y[grep("User defined function value =",y)],"=")[[1]][2])
			DF <- as.numeric(strsplit(y[grep("Degrees of freedom",y)],">")[[1]][17])

   			#Attaching fit functions onto the matrices list
		   cat(paste("matrices$User.fit.",k,  " <- ",User.fit,sep=""),file="temp");eval(parse(file="temp"))
			cat(paste("matrices$DF.",k,  " <- ",DF,sep=""),file="temp");eval(parse(file="temp"))
		}
############ 

 
# XXX CODE SO BROKEN.  Sometimes log-likelihood appears but not the other statistics  
############ 
#Finding the fit function information if "log-likelihood" is in the mxo file
#		if  (length(grep("-2 times log-likelihood of data >",y)) != 0) {
#			LL <- as.numeric(strsplit(y[grep("-2 times log-likelihood of data >",y)],">")[[1]][4])
#			DF <- as.numeric(strsplit(y[grep("Degrees of freedom >",y)],">")[[1]][17])
#			AIC <- as.numeric(strsplit(y[grep("Akaike's Information Criterion >",y)],">")[[1]][5])
#			BIC <- as.numeric(strsplit(y[grep("Bayesian Information Criterion >",y)],">")[[1]][5])
#			Adj.BIC <- as.numeric(strsplit(y[grep("Sample size Adjusted BIC       >",y)],">")[[1]][5])
#			DIC <- as.numeric(strsplit(y[grep("Deviance Information Criterion >",y)],">")[[1]][5])
 
			#Attaching fit functions onto the matrices list
#			cat(paste("matrices$LL.",k,  " <- ",LL,sep=""),file="temp");eval(parse(file="temp"))
#			cat(paste("matrices$DF.",k,  " <- ",DF,sep=""),file="temp");eval(parse(file="temp"))
#			cat(paste("matrices$AIC.",k,  " <- ",AIC,sep=""),file="temp");eval(parse(file="temp"))
#			cat(paste("matrices$BIC.",k,  " <- ",BIC,sep=""),file="temp");eval(parse(file="temp"))
#			cat(paste("matrices$Adj.BIC.",k,  " <- ",Adj.BIC,sep=""),file="temp");eval(parse(file="temp"))
#			cat(paste("matrices$DIC.",k,  " <- ",DIC,sep=""),file="temp");eval(parse(file="temp"))
#		}
############ 

 

############ 
 		#Inner Loop - looping through different matrices within groups
		for (i in mat.starts) {
			j <- j+1
			Grp <- as.numeric(x[(mat.starts[j]+1)])
			Mat <- substr(x[(mat.starts[j]+2)],2,2)
			Type <- as.numeric(x[(mat.starts[j]+3)])
			Dims <- as.numeric(x[(mat.starts[j]+4):(mat.starts[j]+5)])
			num.elements <- as.numeric(x[mat.starts[j]+7])
			if (length(start.remove)>0 & as.numeric(x[mat.starts[j]+7])==100 & (Dims[1]*Dims[2])> 100) {
	        	num.elements <- Dims[1]*Dims[2]
			}
    		if (Type !=1 & Type !=2 & Type !=3 & Type !=4 & Type != 11) {
        		elements <- as.numeric(x[(mat.starts[j]+8):(mat.starts[j]+7+num.elements)])
			}
    
    		#Calculated, Constructed, or Full Matrices
    		if (Type <= -1|Type==9){
        		cat(paste(Mat,Grp,".",k," <- matrix(elements,nrow=Dims[1],ncol=Dims[2],byrow=TRUE)",sep=""),file="temp")
			}
    		#Zero Matrix
    		if (Type==1) {
        		cat(paste(Mat,Grp,".",k," <- matrix(0,nrow=Dims[1],ncol=Dims[2])",sep=""),file="temp")
			}
		    #Identity Matrix
			if (Type==2) {
				cat(paste(Mat,Grp,".",k," <- diag(Dims[1])",sep=""),file="temp")
			}
			#Identity-Zero Matrix
			if (Type==3) {
				ident <- diag(Dims[1])
				zeros <- matrix(0,nrow=Dims[1],ncol=(Dims[2]-Dims[1]))
				cat(paste(Mat,Grp,".",k," <- cbind(ident,zeros)",sep=""),file="temp")
			}
		    #Zero-Identity Matrix
   			if (Type==4) {
				ident <- diag(Dims[1])
				zeros <- matrix(0,nrow=Dims[1],ncol=(Dims[2]-Dims[1]))
				cat(paste(Mat,Grp,".",k," <- cbind(zeros,ident)",sep=""),file="temp")
			}
			#Diagonal Matrix
			if (Type==5) {
				cat(paste(Mat,Grp,".",k," <- diag(elements)",sep=""),file="temp")
			}
			#Subdiagonal Matrix
			if (Type==6){
				tempmat <- matrix(0,nrow=Dims[1],ncol=Dims[2])
				tempmat[upper.tri(tempmat)] <- elements
				tempmat <- t(tempmat)
				cat(paste(Mat,Grp,".",k," <- tempmat",sep=""),file="temp")
			}
			#Standardized Matrix
			if (Type==7) {
				tempmat <- matrix(1,nrow=Dims[1],ncol=Dims[2])
				tempmat[upper.tri(tempmat)] <- elements
				tempmat[lower.tri(tempmat)] <- elements
				tempmat <- t(tempmat)
				cat(paste(Mat,Grp,".",k," <- tempmat",sep=""),file="temp")
			}
			#Symmetric Matrix
			if (Type==8) {
				tempmat <- matrix(0,nrow=Dims[1],ncol=Dims[2])
				tempmat[upper.tri(tempmat,diag=TRUE)] <- elements
				tempmat <- t(tempmat)
 				tempmat[upper.tri(tempmat)] <- tempmat[lower.tri(tempmat)]
				cat(paste(Mat,Grp,".",k," <- tempmat",sep=""),file="temp")
			}
		    #Lower Matrix
			if (Type==10){
				tempmat <- matrix(0,nrow=Dims[1],ncol=Dims[2])
				tempmat[upper.tri(tempmat,diag=TRUE)] <- elements
				tempmat <- t(tempmat)
				cat(paste(Mat,Grp,".",k," <- tempmat",sep=""),file="temp")
			}
			#Unit Matrix
			if (Type==11) {
		        cat(paste(Mat,Grp,".",k," <- matrix(1,nrow=Dims[1],ncol=Dims[2])",sep=""),file="temp")
			}

		    #Make the new matrix & latch it onto the "matrices" list			eval(parse(file="temp"))
			cat(paste("matrices$",Mat,Grp,".",k,  " <- ",Mat,Grp,".",k,sep=""),file="temp")
			eval(parse(file="temp"))
############ 
		} # Inner loop - matrices
	} # Outer loop - jobs

	return(matrices)
}

######################################################################## 
 


