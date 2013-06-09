options(error = utils::recover)   # uncomment for more help with debugging
library(OpenMx)
library(mvtnorm)
library(plyr)

set.seed(1)

#more.noise <- 0
more.noise <- 1

gen.data <- function(n) {
  data.cov <- matrix(c(1, .2, .2, 1), byrow=TRUE, nrow=2)
  latent <- rmvnorm(n, mean=c(0,0), sigma=data.cov)
  colnames(latent) <- c("A","B")
  latent <- as.data.frame(latent)
  df <- data.frame(C=latent$A + latent$B,
                   D=latent$A - latent$B)
  if (more.noise) {
    df$C <- df$C + rnorm(1, sd=more.noise)
    df$D <- df$D + rnorm(1, sd=more.noise)
  }
  df
}

fanout <- 20

school.data <- cbind(id=1:fanout, gen.data(fanout))
#school.data$C <- school.data$id * 1000
teacher.data <- cbind(schoolId=1:fanout, id=seq(1,fanout^2), gen.data(fanout^2))
#teacher.data$C <- teacher.data$id * 100
student.data <- cbind(teacherId=seq(1,fanout^2), id=seq(1,fanout^3), gen.data(fanout^3))

stack.data <- function (key, upper, lower) {  
  lower <- ddply(lower, key, function(slice) {
    id <- unique(slice[[key]])
    parent <- upper[upper$id == id,]
    slice$C <- slice$C + parent$C
    slice$D <- slice$D + parent$C
    slice
  })
}
teacher.data <- stack.data("schoolId", school.data, teacher.data)
student.data <- stack.data("teacherId", teacher.data, student.data)

manifests<-c("C","D")
latents<-c("A","B")
student <- mxModel("student",
                 type="RAM",
                 manifestVars = manifests,
                 latentVars = latents,
                 mxPath(from="A",to=c("C","D"), free=c(FALSE,FALSE), value=c(1,1) , arrows=1,
                        label=c("A_TO_C","A_TO_D") ),
                 mxPath(from="B",to=c("C","D"), free=c(FALSE,FALSE), value=c(1,-1) ,
                        arrows=1, label=c("B_TO_C","B_TO_D") ),
                 mxPath(from="A",to=c("A","B"), free=c(TRUE,TRUE), value=c(1,0) , arrows=2,
                        label=c("VAR_A","COV_A_B") ),
                 mxPath(from="B",to=c("B"), free=c(TRUE), value=c(1) , arrows=2,
                        label=c("VAR_B") ),
                 mxPath(from="C",to=c("C"), free=as.logical(more.noise), value=more.noise, arrows=2,
                        label=c("VAR_C") ),
                 mxPath(from="D",to=c("D"), free=as.logical(more.noise), value=more.noise, arrows=2,
                        label=c("VAR_D") ),
                   mxPath(from="one", to=c(manifests, latents), value=0, free=FALSE)
);

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

school <- mxModel(school, mxData(school.data, type="raw", primaryKey="id"))
teacher <- mxModel(teacher,
                   mxData(teacher.data, type="raw",
                          primaryKey="id",
                          foreignKeys=list(c('schoolId', 'school.id'))))

student <- mxModel(student,
                   mxData(student.data, type="raw",
                          primaryKey="id",
                          foreignKeys=list(c('teacherId', 'teacher.id'))))

sonly <- mxModel(student, type="RAM", name="sonly")
sonly <- mxRun(sonly)
summary(sonly)

district <- mxModel("district",
                     type="RAM",
                     school, teacher, student,
                     
		    mxPath(from="school.C", to=c("teacher.A"), free=FALSE, value=1),
		    mxPath(from="teacher.C", to="student.A", free=FALSE, value=1),
                     
                     mxExpectationRAM(H=paste0("H",1:2), HomerTransform=TRUE),
                     mxFitFunctionML()
)

district <- mxRun(district)
district@submodels[[1]]@matrices$S@values
district@submodels[[2]]@matrices$S@values
district@submodels[[3]]@matrices$S@values

