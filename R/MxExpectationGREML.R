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

setClass(Class = "MxExpectationGREML",
         slots = c(
           V = "MxCharOrNumber",
           yvars = "character",
           Xvars = "list",
           addOnes="logical", 
           blockByPheno="logical",
           staggerZeroes="logical",
           dataset.is.yX="logical",
           X="matrix",
           y="MxData",
           yXcolnames="character",
           casesToDrop="integer",
           b="matrix",
           bcov="matrix",
           numFixEff = "integer",
           dims = "character",
           numStats = "numeric",
           dataColumns = "numeric",
           name = "character"),
         contains = "MxBaseExpectation")

#Since there is a small chance that future developers or sophisticated users might create MxExpectationGREML
#objects with new() instead of mxExpectationGREML(), the class constructor should provide defaults for all 
#the slots...
setMethod("initialize", "MxExpectationGREML",
          function(.Object, V=character(0), yvars=character(0), Xvars=list(), addOnes=TRUE, 
                   blockByPheno=TRUE, staggerZeroes=TRUE, dataset.is.yX=FALSE, casesToDrop=integer(0),
                   data = as.integer(NA), name = 'expectation') {
            .Object@name <- name
            .Object@V <- V
            .Object@yvars <- yvars
            .Object@Xvars <- Xvars
            .Object@addOnes <- addOnes
            .Object@blockByPheno <- blockByPheno
            .Object@staggerZeroes <- staggerZeroes
            .Object@dataset.is.yX <- dataset.is.yX
            .Object@numFixEff <- integer(0)
            .Object@casesToDrop <- casesToDrop
            .Object@data <- data
            .Object@X <- matrix(as.numeric(NA),1,1)
            .Object@dims <- "foo"
            return(.Object)
          }
)


setMethod("qualifyNames", signature("MxExpectationGREML"), 
          function(.Object, modelname, namespace) {
            .Object@name <- imxIdentifier(modelname, .Object@name)
            .Object@V <- imxConvertIdentifier(.Object@V, modelname, namespace)
            .Object@data <- imxConvertIdentifier(.Object@data, modelname, namespace)
            return(.Object)
          })


setMethod("genericExpDependencies", signature("MxExpectationGREML"),
          function(.Object, dependencies) {
            sources <- c(.Object@V)
            sources <- sources[!is.na(sources)]
            dependencies <- imxAddDependency(sources, .Object@name, dependencies)
            return(dependencies)
          })


setMethod("genericExpAddEntities", "MxExpectationGREML",
          function(.Object, job, flatJob, labelsData) {return(job)}
)


#setMethod("genericExpConvertEntities", "MxExpectationGREML",


setMethod("genericExpRename", signature("MxExpectationGREML"),
          function(.Object, oldname, newname) {
            .Object@V <- renameReference(.Object@V, oldname, newname)
            .Object@data <- renameReference(.Object@data, oldname, newname)
            return(.Object)
          })

setMethod("genericGetExpected", signature("MxExpectationGREML"),
	  function(.Object, model, what, defvar.row=1) {
		  if ('covariance' %in% what) {
			  covname <- .Object@V
			  cov <- mxEvalByName(covname, model, compute=TRUE, defvar.row=defvar.row)
			  ret[['covariance']] <- cov
		  }
		  if ('means' %in% what) {
			  ret[['means']] <- NA
		  }
		  if ('thresholds' %in% what) {
			  ret[['thresholds']] <- NA
		  }
		  ret
	  })

mxExpectationGREML <- function(V, yvars=character(0), Xvars=list(), addOnes=TRUE, blockByPheno=TRUE, 
                               staggerZeroes=TRUE, dataset.is.yX=FALSE, casesToDropFromV=integer(0)){
  blockByPheno <- as.logical(blockByPheno)[1]
  staggerZeroes <- as.logical(staggerZeroes)[1]
  addOnes <- as.logical(addOnes)[1]
  dataset.is.yX <- as.logical(dataset.is.yX)[1]
  if (missing(V) || typeof(V) != "character") {
    stop("argument 'V' is not of type 'character' (the name of the expected covariance matrix)")
  }
  if(!dataset.is.yX){
    casesToDropFromV <- integer(0) #<--Ignore casesToDropFromV unless dataset.is.yX is true
    if ( missing(yvars) || typeof(yvars) != "character" )  {
      stop("argument 'yvars' is not of type 'character' (the data column names of the phenotypes)")
    }
    if(!length(yvars)){
      stop("you must specify at least one phenotype in argument 'yvars'")
    }
    if( !is.list(Xvars) ){
      if(length(yvars)==1){Xvars <- list(Xvars)}
      else{stop("argument 'Xvars' must be provided as a list when argument 'yvars' is of length greater than 1")}
    }
    if(length(Xvars)){
      if( !all(sapply(Xvars,is.character)) ){
        stop("elements of argument 'Xvars' must be of type 'character' (the data column names of the covariates)")
      }
      if( !staggerZeroes && !all(sapply(Xvars,length)==sapply(Xvars,length)[1]) ){
        stop("all phenotypes must have the same number of covariates when staggerZeroes=FALSE")
      }
      if(length(Xvars)!=length(yvars)){
        #In the polyphenotype case, the same covariates will often be used for all phenotypes:
        if(length(Xvars)<length(yvars) && length(Xvars)==1){Xvars <- rep(Xvars,length.out=length(yvars))}
        else{stop("conflicting number of phenotypes specified by arguments 'Xvars' and 'yvars'")}
      }
  }}
  return(new("MxExpectationGREML", V, yvars, Xvars, addOnes, blockByPheno, staggerZeroes, dataset.is.yX, 
             casesToDrop=casesToDropFromV))
}


setMethod("genericExpFunConvert", "MxExpectationGREML", 
          function(.Object, flatModel, model, labelsData, dependencies) {
            modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
            name <- .Object@name
            #There just needs to be something in the data slot, since the backend expects it:
            if(is.na(.Object@data)){
              msg <- paste("the GREML expectation function",
                           "does not have a dataset associated with it in model",
                           omxQuotes(modelname))
              stop(msg, call. = FALSE)
            }
            mxDataObject <- flatModel@datasets[[.Object@data]]
            checkNumericData(mxDataObject)
            if (mxDataObject@type != "raw") {
              stop("GREML expectation only compatible with raw data",call.=F)
            }
            if(sum(apply(mxDataObject@observed, 2, is.factor))>0){
              stop("GREML expectation not compatible with ordinal data", call.=F)
            }
            if(!length(colnames(mxDataObject@observed))){
              msg <- paste("dataset does not have column names in model",omxQuotes(modelname))
              stop(msg, call. = FALSE)
            }
            if(.Object@dataset.is.yX){
              .Object@y <- mxData(
                observed=matrix(mxDataObject@observed[,1], nrow=1,
                                dimnames=list(NULL,paste("y",1:nrow(mxDataObject@observed),sep=""))),
                type="raw",sort=FALSE)
              .Object@X <- as.matrix(mxDataObject@observed[,-1])
              .Object@yXcolnames <- colnames(mxDataObject@observed)
              .Object@numFixEff <- as.integer(ncol(mxDataObject@observed)-1)
              .Object@dataColumns <- as.double(0:(nrow(mxDataObject@observed)-1))
            }
            else{
              if(length(.Object@Xvars)){
                .Object@numFixEff <- as.integer(sum(sapply(.Object@Xvars, length)) + (length(.Object@yvars) * .Object@addOnes))
              } #If Xvars is length 0, then the only covariates will be 1s (for the intercepts)
              else{.Object@numFixEff <- as.integer(length(.Object@yvars))}
              if( !all(.Object@yvars %in% colnames(mxDataObject@observed)) ){
                badname <- (.Object@yvars[!(.Object@yvars %in% colnames(mxDataObject@observed))])[1]
                msg <- paste("'",badname,"' is not among the data column names",sep="")
                stop(msg)
              }
              if( length(.Object@Xvars) && !all(unlist(.Object@Xvars) %in% colnames(mxDataObject@observed)) ){
                badname <- (unlist(.Object@Xvars)[!(unlist(.Object@Xvars) %in% colnames(mxDataObject@observed))])[1]
                msg <- paste("'",badname,"' is not among the data column names",sep="")
                stop(msg)
              }
              mm <- mxGREMLDataHandler(data=mxDataObject@observed, yvars=.Object@yvars, Xvars=.Object@Xvars, 
                                     addOnes=.Object@addOnes, blockByPheno=.Object@blockByPheno, 
                                     staggerZeroes=.Object@staggerZeroes)
              .Object@y <- mxData(
                observed=matrix(mm$yX[,1], nrow=1,
                                dimnames=list(NULL,paste("y",1:length(mm$yX[,1]),sep=""))),
                type="raw",sort=FALSE)
              .Object@X <- as.matrix(mm$yX[,-1])
              .Object@yXcolnames <- colnames(mm$yX)
              .Object@casesToDrop <- mm$casesToDrop
              .Object@numFixEff <- ncol(.Object@X)
              .Object@dataColumns <- as.double(0:(nrow(.Object@X)-1))
            }
            #Get number of observed statistics BEFORE call to backend, so summary() can use it:
            .Object@numStats <- nrow(.Object@X)
            dataName <- .Object@data
            .Object@data <- imxLocateIndex(flatModel, .Object@data, name)
            .Object@V <- imxLocateIndex(flatModel, .Object@V, name)
            return(.Object)
          })


mxGREMLDataHandler <- function(data, yvars=character(0), Xvars=list(), addOnes=TRUE, blockByPheno=TRUE, 
                               staggerZeroes=TRUE){
  
  #Input checking:
  blockByPheno <- as.logical(blockByPheno)[1]
  staggerZeroes <- as.logical(staggerZeroes)[1]
  addOnes <- as.logical(addOnes)[1]
  if( !is.matrix(data) && !is.data.frame(data) ){
    stop("argument 'data' must be either a matrix or dataframe")
  }
  if(!length(colnames(data))){stop("data must have column names")}
  if ( missing(yvars) || typeof(yvars) != "character" )  {
    stop("argument 'yvars' is not of type 'character' (the data column names of the phenotypes)")
  }
  if(!length(yvars)){
    stop("you must specify at least one phenotype in argument 'yvars'")
  }
  if( !is.list(Xvars) ){
    if(length(yvars)==1){Xvars <- list(Xvars)}
    else{stop("argument 'Xvars' must be provided as a list when argument 'yvars' is of length greater than 1")}
  }
  if( !all(yvars %in% colnames(data)) ){
    badname <- (yvars[!(yvars %in% colnames(data))])[1]
    msg <- paste("'",badname,"' in argument 'yvars' is not among the data column names",sep="")
    stop(msg)
  }
  if(length(Xvars)){
    if( !all(sapply(Xvars,is.character)) ){
      stop("elements of argument 'Xvars' must be of type 'character' (the data column names of the covariates)")
    }
    if( !all(unlist(Xvars) %in% colnames(data)) ){
      badname <- (unlist(Xvars)[!(unlist(Xvars) %in% colnames(data))])[1]
      msg <- paste("'",badname,"' in argument 'Xvars' is not among the data column names",sep="")
      stop(msg)
    }
    if( !staggerZeroes && !all(sapply(Xvars,length)==sapply(Xvars,length)[1]) ){
      stop("all phenotypes must have the same number of covariates when staggerZeroes=FALSE")
    }
    if(length(Xvars)!=length(yvars)){
      #In the polyphenotype case, the same covariates will often be used for all phenotypes:
      if(length(Xvars)<length(yvars) && length(Xvars)==1){Xvars <- rep(Xvars,length.out=length(yvars))}
      else{stop("conflicting number of phenotypes specified by arguments 'Xvars' and 'yvars'")}
    }
  }
  
  #Handle phenotypes:
  y <- NULL
  i <- 1
  if(blockByPheno){ #Stack phenotypes
    while(i <= length(yvars)){
      y <- rbind(y,as.matrix(data[,yvars[i]]))
      i <- i+1
  }}
  else{ #Stack cases
    for(i in 1:nrow(data)){
      y <- rbind(y,matrix(data[i,yvars]))
  }}
  if(length(yvars)==1){colnames(y) <- yvars}
  else{colnames(y) <- "y"}
  
  #Assemble matrix of covariates:
  X <- NULL
  i <- 1
  if(!length(Xvars)){
    X <- diag(length(yvars)) %x% matrix(1,nrow=nrow(data),ncol=1)
    colnames(X) <- paste("x",1:length(yvars),sep="")
  }
  else{
    if(staggerZeroes){
      if(blockByPheno){
        while(i <= length(yvars)){
          ncolprev <- ncol(X)
          Xcurr <- as.matrix(data[ ,Xvars[[i]] ])
          if(addOnes){Xcurr <- cbind(1,Xcurr)}
          if(i==1){X <- Xcurr}
          else{
            X <- rbind(
              cbind( X, matrix(0,nrow(X),ncol(Xcurr)) ),
              cbind( matrix(0,nrow(Xcurr),ncol(X)), Xcurr )
            )
          }
          if(length(yvars)==1){
            if(addOnes){colnames(X) <- c("Intrcpt",Xvars[[1]])}
            else{colnames(X) <- Xvars[[1]]}
          }
          else{
            if(i==1){
              if(addOnes){colnames(X) <- paste(yvars[1], c("1",Xvars[[1]]), sep="_")}
              else{colnames(X) <- paste(yvars[1], c(Xvars[[1]]), sep="_")}
            }
            else{
              if(addOnes){colnames(X)[(ncolprev+1):ncol(X)] <- paste(yvars[i], c("1",Xvars[[i]]), sep="_")}
              else{colnames(X)[(ncolprev+1):ncol(X)] <- paste(yvars[i], c(Xvars[[i]]), sep="_")}
            }}
          i <- i+1
        }}
      else{ #if !blockByPheno
        while(i <= length(yvars)){
          k <- i
          Xcurr <- matrix(0,nrow(data)*length(yvars),ncol=length(Xvars[[i]])+addOnes)
          if(length(yvars)==1){
            if(addOnes){colnames(Xcurr) <- c("Intrcpt",Xvars[[1]])}
            else{colnames(Xcurr) <- Xvars[[1]]}
          }
          else{
            if(addOnes){colnames(Xcurr) <- paste(yvars[i], c("1",Xvars[[i]]), sep="_")}
            else{colnames(Xcurr) <- paste(yvars[i], c(Xvars[[i]]), sep="_")}
          }
          for(j in 1:nrow(data)){
            if(addOnes){Xcurr[k,] <- c(1,as.vector(data[j,Xvars[[i]]]))}
            else{Xcurr[k,] <- as.vector(data[j,Xvars[[i]]])}
            k <- k+length(yvars)
          }
          if(i==1){X <- Xcurr}
          else{X <- cbind(X,Xcurr)}
          i <- i+1
        }}}
    else{ #if !staggerZeroes
      if(blockByPheno){
        while(i <= length(yvars)){
          Xcurr <- as.matrix(data[ ,Xvars[[i]] ])
          if(addOnes){
            Xcurr <- cbind(1,Xcurr)
            colnames(Xcurr) <- c("Intrcpt",Xvars[[1]])
          }
          else{colnames(Xcurr) <- Xvars[[1]]}
          if(i==1){X <- Xcurr}
          else{X <- rbind(X,Xcurr)}
          i <- i+1
        }}
      else{
        for(i in 1:nrow(data)){
          Xcurr <- matrix(NA,nrow=length(yvars),ncol=length(Xvars[[1]]))
          for(j in 1:length(yvars)){
            Xcurr[j,] <- as.vector(data[i,Xvars[[j]]])
          }
          if(addOnes){Xcurr <- cbind(1,Xcurr)}
          if(i==1){
            if(addOnes){colnames(Xcurr) <- c("Intrcpt",Xvars[[1]])}
            else{colnames(Xcurr) <- Xvars[[1]]}
            X <- Xcurr
          }
          else{X <- rbind(X,Xcurr)}
  }}}}
  if( length(unique(colnames(X))) < ncol(X) ){
    warning("resulting 'X' matrix has some redundant column names; is it full rank?")
  }
  
  y.colnames <- colnames(y)
  X.colnames <- colnames(X)
  
  #Identify which subjects have incomplete data:
  whichHaveNA <- which(as.logical(rowSums(is.na(cbind(y,X)))))
  if(length(whichHaveNA)){
    y <- as.matrix(y[-whichHaveNA,])
    colnames(y) <- y.colnames
    X <- as.matrix(X[-whichHaveNA,])
    colnames(X) <- X.colnames
  }
  
  return(list(yX=cbind(y,X),casesToDrop=whichHaveNA))
}



#TODO: Make sure the fixed-effects results look OK in summary() output for multigroup model
GREMLFixEffList <- function(model) {
  if(length(model@submodels) > 0) {
    ptable <- vector("list",length(model@submodels)+1)
    names(ptable) <- c(model$name, names(model@submodels))
    ptable[2:length(ptable)] <- lapply(model@submodels,GREMLFixEffList)
    ptable[[1]] <- GREMLFixEffListHelper(model)
    ptable <- ptable[sapply(ptable,function(x){length(x)>0})]
  } 
  else{
    ptable <- GREMLFixEffListHelper(model)
  }
  if(length(ptable)==0){ptable <- NULL}
  return(ptable)
}

GREMLFixEffListHelper <- function(model) {
  ptable <- NULL
  if( (length(model@output)==0) || (length(model$expectation$b)==0) ){ return(ptable) }
  if(length(model$expectation$yXcolnames) > 0){
    ptable <- data.frame(name = model$expectation$yXcolnames[-1], stringsAsFactors = F)
  }
  else{ptable <- data.frame(name = paste("x", 0:(length(model$expectation$b)-1), sep=""), 
                            stringsAsFactors = F)}
  ptable$coeff <- model$expectation$b
  ptable$se <- sqrt(diag(model$expectation$bcov))
  return(ptable)
}
