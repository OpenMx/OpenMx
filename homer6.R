#options(error = utils::recover)   # uncomment for more help with debugging
library(OpenMx)
library(mvtnorm)
library(plyr)

mxOption(NULL, "Number of Threads", 1)

set.seed(1)

gen.data <- function(n) {
    data.frame(A=rnorm(n))
}

fanout <- 3

school.data <- cbind(id=1, gen.data(1))
#school.data$C <- school.data$id * 1000
teacher.data <- cbind(schoolId=1, id=seq(1,fanout), gen.data(fanout))

manifests<-c("A")
student <- mxModel("student",
                 type="RAM",
                 manifestVars = manifests,
                 mxPath(from="A",to=c("A"), free=c(TRUE), value=c(.9), arrows=2,
                        label=c("VAR_A") ),
                   mxPath(from="one", to=c(manifests), value=.1, free=FALSE)
)

relabel <- function(m, prefix) {
  for (mat in c("A","S")) {
    lab <- m@matrices[[mat]]@labels
    lab[!is.na(lab)] <- paste0(prefix, lab[!is.na(lab)])
    m@matrices[[mat]]@labels <- lab
  }
  m
}

teacher <- relabel(mxModel(student, name="teacher"), "tea_")
school <- relabel(mxModel(student, name="school"), "sch_")
student <- relabel(student, "st_")

print(school.data)
print(teacher.data)

school <- mxModel(school, mxData(school.data, type="raw", primaryKey="id"))
teacher <- mxModel(teacher,
                   mxData(teacher.data, type="raw",
                          primaryKey="id",
                          foreignKeys=list(c('schoolId', 'school.id'))))

district <- mxModel("district",
                     type="RAM",
                     school, teacher,
                     
		    mxPath(from="school.A", to=c("teacher.A"), free=FALSE, value=1),
                     
                     mxExpectationRAM(H=paste0("H",1), HomerTransform=TRUE),
                     mxFitFunctionML()
)

district <- mxRun(district, useOptimizer=FALSE)

#district@submodels[[1]]@matrices$S@values
#district@submodels[[2]]@matrices$S@values
#district@submodels[[3]]@matrices$S@values

district@output$minimum

stop("here")

modelLL <- function(obs, mmean, A, S) {
	iA <- solve(diag(nrow(A)) - A)
	dmvnorm(t(obs - iA %*% mmean), sigma=iA %*% S %*% t(iA), log=TRUE)
}

m1LL <- function(m) {
    obs <- m$data@observed[,rownames(m$S@values)]
    modelMean <- matrix(m$M@values, byrow=TRUE, ncol=length(m$M@values), nrow=nrow(obs))
    modelLL(obs, modelMean,
            kronecker(diag(nrow(obs)), m$A@values),
            kronecker(diag(nrow(obs)), m$S@values))
}

m1LL(district$school)

debug(modelLL)

m <- district$teacher
