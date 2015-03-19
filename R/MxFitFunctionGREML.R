#
#   Copyright 2007-2015 The OpenMx Project
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

setClass(Class = "MxFitFunctionGREML", 
         slots=c(
           casesToDrop="integer",
           dropNAfromV = "logical",
           numObs = "integer"),
         contains = "MxBaseFitFunction")


setMethod("initialize", "MxFitFunctionGREML",
          function(.Object, name = 'fitfunction', casesToDrop=integer(0), dropNAfromV=logical(0)) {
            .Object@name <- name
            .Object@vector <- FALSE
            .Object@numObs <- 0L
            .Object@casesToDrop <- casesToDrop
            .Object@dropNAfromV <- dropNAfromV
            return(.Object)
          }
)


setMethod("qualifyNames", signature("MxFitFunctionGREML"), 
          function(.Object, modelname, namespace) {
            .Object@name <- imxIdentifier(modelname, .Object@name)
            return(.Object)
          })

setMethod("genericFitConvertEntities", "MxFitFunctionGREML",
          function(.Object, flatModel, namespace, labelsData) {
            
            name <- .Object@name
            modelname <- imxReverseIdentifier(flatModel, .Object@name)[[1]]
            expectName <- paste(modelname, "expectation", sep=".")
            
            expectation <- flatModel@expectations[[expectName]]
            dataname <- expectation@data		
            
            return(flatModel)
          })


setMethod("genericFitFunConvert", "MxFitFunctionGREML", 
          function(.Object, flatModel, model, labelsData, defVars, dependencies) {
            name <- .Object@name
            modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
            expectName <- paste(modelname, "expectation", sep=".")
            if (expectName %in% names(flatModel@expectations)) {
              expectIndex <- imxLocateIndex(flatModel, expectName, name)
            } else {
              expectIndex <- as.integer(NA)
            }
            .Object@expectation <- expectIndex
            return(.Object)
          })


setMethod("genericFitInitialMatrix", "MxFitFunctionGREML",
          function(.Object, flatModel) {return(matrix(as.double(NA), 1, 1))})


mxFitFunctionGREML <- function(casesToDrop=integer(0), dropNAfromV=TRUE){
  return(new("MxFitFunctionGREML",casesToDrop=casesToDrop,dropNAfromV=dropNAfromV))
}


mxGREMLStarter <- function(model, data, Xdata, ydata, Xname="X", yname="y", 
                           addOnes=TRUE, blockByPheno=TRUE, staggerZeroes=TRUE, dropNAfromV=TRUE){
  
  #Input checks:
  blockByPheno <- as.logical(blockByPheno)[1]
  staggerZeroes <- as.logical(staggerZeroes)[1]
  dropNAfromV <- as.logical(dropNAfromV)[1]
  addOnes <- as.logical(addOnes)[1]
  if( !(class(model) %in% c("character","MxModel")) ){
    stop("argument 'model' must be either a character string or an MxModel object")
  }
  if( !is.matrix(data) && !is.data.frame(data) ){
    stop("argument 'data' must be either a matrix or dataframe")
  }
  if(!length(colnames(data))){stop("data must have column names")}
  if( !is.list(Xdata) ){
    if(length(ydata)==1){Xdata <- list(Xdata)}
    else{stop("argument 'Xdata' must be provided as a list when argument 'ydata' is of length greater than 1")}
  }
  if( !all(sapply(Xdata,is.character)) ){
    stop("elements of argument 'Xdata' must be of type 'character' (the data column names of the covariates)")
  }
  if(!is.character(ydata)){
    stop("argument 'ydata' must be of type 'character' (the data column names of the phenotypes)")
  }
  if( !all(ydata %in% colnames(data)) ){
    badname <- (ydata[!(ydata %in% colnames(data))])[1]
    msg <- paste("'",badname,"' in argument 'ydata' is not among the data column names",sep="")
    stop(msg)
  }
  if( !all(unlist(Xdata) %in% colnames(data)) ){
    badname <- (unlist(Xdata)[!(unlist(Xdata) %in% colnames(data))])[1]
    msg <- paste("'",badname,"' in argument 'Xdata' is not among the data column names",sep="")
    stop(msg)
  }
  if(!is.character(Xname)){
    stop("argument 'Xname' is not of type 'character' (the name for the matrix of covariates)")
  }
  if(!is.character(yname)){
    stop("argument 'yname' is not of type 'character' (the name for the column vector of phenotypes)")
  }
  if( ("MxModel" %in% class(model)) && (Xname %in% c(names(model@matrices),names(model@algebras))) ){
    msg <- paste("already an MxMatrix or MxAlgebra named '",Xname,"' in model '",model$name,"'",sep="")
    stop(msg)
  }
  if( ("MxModel" %in% class(model)) && (yname %in% c(names(model@matrices),names(model@algebras))) ){
    msg <- paste("already an MxMatrix or MxAlgebra named '",yname,"' in model '",model$name,"'",sep="")
    stop(msg)
  }
  if( !staggerZeroes && !all(sapply(Xdata,length)==sapply(Xdata,length)[1]) ){
    stop("all phenotypes must have the same number of covariates when staggerZeroes=TRUE")
  }
  if(length(Xdata)!=length(ydata)){
    #In the polyphenotype case, the same covariates will often be used for all phenotypes:
    if(length(Xdata)<length(ydata)){Xdata <- rep(Xdata,length.out=length(ydata))}
    else{stop("conflicting number of phenotypes specified by arguments 'Xdata' and 'ydata'")}
  }
  
  #Handle phenotypes:
  y <- NULL
  i <- 1
  if(blockByPheno){ #Stack phenotypes
    while(i <= length(ydata)){
      y <- rbind(y,as.matrix(data[,ydata[i]]))
      i <- i+1
  }}
  else{ #Stack cases
    for(i in 1:nrow(data)){
      y <- rbind(y,as.matrix(data[i,ydata]))
  }}
  
  #Assemble matrix of covariates:
  X <- NULL
  i <- 1
  if(staggerZeroes){
    if(blockByPheno){
      while(i <= length(ydata)){
        ncolprev <- ncol(X)
        Xcurr <- as.matrix(data[ ,Xdata[[i]] ])
        if(addOnes){Xcurr <- cbind(1,Xcurr)}
        if(i==1){X <- Xcurr}
        else{
          X <- rbind(
            cbind( X, matrix(0,nrow(X),ncol(Xcurr)) ),
            cbind( matrix(0,nrow(Xcurr),ncol(X)), Xcurr )
          )
        }
        if(length(ydata)==1){
          if(addOnes){colnames(X) <- c("1",Xdata[[1]])}
          else{colnames(X) <- Xdata[[1]]}
        }
        else{
          if(i==1){
            if(addOnes){colnames(X) <- paste(ydata[1], c("1",Xdata[[1]]), sep="_")}
            else{colnames(X) <- paste(ydata[1], c(Xdata[[1]]), sep="_")}
          }
          else{
            if(addOnes){colnames(X)[(ncolprev+1):ncol(X)] <- paste(ydata[i], c("1",Xdata[[i]]), sep="_")}
            else{colnames(X)[(ncolprev+1):ncol(X)] <- paste(ydata[i], c(Xdata[[i]]), sep="_")}
        }}
        i <- i+1
    }}
    else{ #if !blockByPheno
      while(i <= length(ydata)){
        k <- i
        Xcurr <- matrix(0,nrow(data)*length(ydata),ncol=length(Xdata[[i]])+addOnes)
        if(length(ydata)==1){
          if(addOnes){colnames(Xcurr) <- c("1",Xdata[[1]])}
          else{colnames(Xcurr) <- Xdata[[1]]}
        }
        else{
          if(addOnes){colnames(Xcurr) <- paste(ydata[i], c("1",Xdata[[i]]), sep="_")}
          else{colnames(Xcurr) <- paste(ydata[i], c(Xdata[[i]]), sep="_")}
        }
        for(j in 1:nrow(data)){
          if(addOnes){Xcurr[k,] <- c(1,as.vector(data[j,Xdata[[i]]]))}
          else{Xcurr[k,] <- as.vector(data[j,Xdata[[i]]])}
          k <- k+length(ydata)
        }
        if(i==1){X <- Xcurr}
        else{X <- cbind(X,Xcurr)}
        i <- i+1
  }}}
  else{ #if !staggerZeroes
    if(blockByPheno){
      while(i <= length(ydata)){
        Xcurr <- as.matrix(data[ ,Xdata[[i]] ])
        if(addOnes){
          Xcurr <- cbind(1,Xcurr)
          colnames(Xcurr) <- c("1",Xdata[[1]])
        }
        else{colnames(Xcurr) <- Xdata[[1]]}
        if(i==1){X <- Xcurr}
        else{X <- rbind(X,Xcurr)}
        i <- i+1
      }}
    else{
      if(length(ydata)){
        for(i in 1:nrow(data)){
          Xcurr <- matrix(NA,nrow=length(ydata),ncol=length(Xdata[[1]]))
          for(j in 1:length(ydata)){
            Xcurr[j,] <- as.vector(data[i,Xdata[[j]]])
          }
          if(addOnes){Xcurr <- cbind(1,Xcurr)}
          if(i==1){
            if(addOnes){colnames(Xcurr) <- c("1",Xdata[[1]])}
            else{colnames(Xcurr) <- Xdata[[1]]}
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
  if( ("MxModel" %in% class(model)) && !is.null(model$fitfunction) ){
    msg <- paste("not adding MxFitFunctionGREML because model '",model$name,"' already contains a fitfunction",
                 sep="")
    warning(msg)
    gff <- NULL
  }
  else{gff <- mxFitFunctionGREML(casesToDrop=whichHaveNA, dropNAfromV=dropNAfromV)}
  
  #Create dummy data object (needed for technical reasons):
  if( ("MxModel" %in% class(model)) && !is.null(model$data) ){
    msg <- paste("not adding dummy MxData object because model '",model$name,"' already contains an MxData object",sep="")
    warning(msg)
    dummydata <- NULL
  }
  else{
    dummydata <- 
         mxData(observed = matrix(as.double(NA),1,1,dimnames = list("dummyData","dummyData")), type="raw")
    dummydata@numObs <- 0
  }
  
  #Assemble MxModel to be returned:
  model.out <- mxModel(
      model,
      dummydata,
      mxMatrix(type="Full",nrow=nrow(y),ncol=1,free=F,values=y,dimnames=list(NULL,yname),name=yname,
               condenseSlots=T),
      mxMatrix(type="Full",nrow=nrow(X),ncol=ncol(X),free=F,values=X,dimnames=list(NULL,colnames(X)),name=Xname,
               condenseSlots=T),
      gff
  )
  return(model.out)
}


#TODO: Make sure the fixed-effects results look OK in summary() output for multigroup model
GREMLFixEffList <- function(model) {
  if(length(model@submodels) > 0) {
    ptable <- vector("list",length(model@submodels)+1)
    names(ptable) <- c(model$name, names(model@submodels))
    ptable[2:length(ptable)] <- lapply(model@submodels,GREMLFixEffList)
    ptable[[1]] <- GREMLFixEffListHelper(model)
    ptable <- ptable[sapply(ptable,function(x){!is.null(x)})]
  } 
  else{
    ptable <- GREMLFixEffListHelper(model)
  }
  return(ptable)
}

GREMLFixEffListHelper <- function(model) {
  ptable <- NULL
  if( (length(model@output)==0) || (length(model$fitfunction$info$b)==0) ){ return(ptable) }
  if(length(colnames(model[[model$expectation$X]]))>0){
    ptable <- data.frame(name = colnames(model[[model$expectation$X]]), stringsAsFactors = F)
  }
  else{ptable <- data.frame(name = paste("x", 0:(ncol(model[[model$expectation$X]])-1), sep=""), 
                            stringsAsFactors = F)}
  ptable$coeff <- model$fitfunction$info$b
  ptable$se <- sqrt(diag(model$fitfunction$info$bcov))
  return(ptable)
}
