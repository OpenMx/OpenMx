# This is the original test case that Timo & I wrote back in Spring 2013.

#options(error = utils::recover)   # uncomment for more help with debugging
library(OpenMx)
library(mvtnorm)

set.seed(1)

more.noise <- 0
#more.noise <- 1

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

fanout <- 5

school.data <- cbind(id=1:fanout, gen.data(fanout))
#school.data$C <- school.data$id * 1000
teacher.data <- cbind(schoolId=1:fanout, id=seq(1,fanout^2),
		      gen.data(fanout^2))
#teacher.data$C <- teacher.data$id * 100
student.data <- cbind(teacherId=seq(1,fanout^2),
		      id=seq(1,fanout^3), gen.data(fanout^3))

stack.data <- function(key, upper, lower) {
	for (pk in upper$id) {
		mask <- lower[[key]] == pk
		for (col in c('C','D')) {
			lower[mask, col] <-
				lower[mask, col] + upper[upper$id == pk, 'C']
		}
	}
	lower
}
teacher.data <- stack.data("schoolId", school.data, teacher.data)
student.data <- stack.data("teacherId", teacher.data, student.data)

manifests<-c("C","D")
latents<-c("A","B")
student <- mxModel(
    "student", type="RAM",
    manifestVars = manifests,
    latentVars = latents,
    mxPath(from="A",to=c("C","D"), free=c(FALSE,FALSE),
	   value=c(1,1), arrows=1,
	   label=c("A_TO_C","A_TO_D") ),
    mxPath(from="B",to=c("C","D"), free=c(FALSE,FALSE), value=c(1,-1) ,
	   arrows=1, label=c("B_TO_C","B_TO_D") ),
    mxPath(from="A",to=c("A","B"), free=c(TRUE,TRUE),
	   value=c(1,0), arrows=2,
	   label=c("VAR_A","COV_A_B") ),
    mxPath(from="B",to=c("B"), free=c(TRUE), value=c(1) , arrows=2,
	   label=c("VAR_B") ),
    mxPath(from="C",to=c("C"), free=as.logical(more.noise),
	   value=more.noise, arrows=2, label=c("VAR_C") ),
    mxPath(from="D",to=c("D"), free=as.logical(more.noise),
	   value=more.noise, arrows=2, label=c("VAR_D") ),
    mxPath(from="one", to=c(manifests, latents), value=0, free=FALSE)
);

relabel <- function(m, prefix) {
  for (mat in c("A","S")) {
    lab <- m[[mat]]$labels
    lab[!is.na(lab)] <- paste0(prefix, lab[!is.na(lab)])
    m[[mat]]$labels <- lab
  }
  m
}

teacher <- relabel(mxModel(student, name="teacher"), "tea_")
school <- relabel(mxModel(student, name="school"), "sch_")
student <- relabel(student, "st_")

school <- mxModel(
    school,
    mxData(school.data, type="raw", primaryKey="id"))

teacher <- mxModel(
    teacher, school,
    mxData(teacher.data, type="raw", primaryKey="id"),
    mxPath('school.C', 'A', free=FALSE, value=1, joinKey="schoolId"))

student <- mxModel(
    student, teacher,
    mxData(student.data, type="raw", primaryKey="id"),
    mxPath('teacher.C', 'A', free=FALSE, value=1, joinKey="teacherId"))

#student$expectation$verbose <- 1L

student$expectation$.rampartCycleLimit <- 0L
pt1 <- mxRun(mxModel(
    student,
    mxComputeSequence(list(
	mxComputeOnce('fitfunction', 'fit'),
	mxComputeNumericDeriv(checkGradient=FALSE,
			      iterations=2, hessian=FALSE),
	mxComputeReportDeriv(),
	mxComputeReportExpectation()))))

student$expectation$.rampartCycleLimit <- as.integer(NA)
pt2 <- mxRun(mxModel(
    student,
    mxComputeSequence(list(
	mxComputeOnce('fitfunction', 'fit'),
	mxComputeNumericDeriv(checkGradient=FALSE,
			      iterations=2, hessian=FALSE),
	mxComputeReportDeriv(),
	mxComputeReportExpectation()))))

omxCheckCloseEnough(pt2$expectation$debug$rampartUsage,
		    c((fanout-1)*fanout^2, (fanout-1)*fanout), 1)
omxCheckCloseEnough(pt2$expectation$debug$numGroups, 3)

if (0) {
	layout <- pt2$expectation$debug$layout
	head(layout[layout$group==3, ],n=20)
}

omxCheckCloseEnough(pt1$output$fit, pt2$output$fit, 1e-7)
omxCheckCloseEnough(pt1$output$gradient, pt2$output$gradient, 1e-6)

student <- mxRun(student)
if (!more.noise) {
	omxCheckCloseEnough(student$output$fit, 1055.161, 1e-2)
} else {
	omxCheckCloseEnough(student$output$fit, 1132.713, 1e-2)  # but code RED
}
#print(student$expectation$debug$rampartUsage)

if (.Platform$OS.type != 'windows' && parallel::detectCores() > 1) {
	omxCheckTrue(student$compute$steps[['GD']]$output$maxThreads > 1)
}

if (0) {
	ex <- student$expectation
	eo = ex$output
	ed = ex$debug
	ed$layout
}

got <- mxGenerateData(student)
omxCheckTrue(setequal(names(got), c("school", "teacher", "student")))
omxCheckTrue(setequal(colnames(got[['school']]),
                      colnames(student$school$data$observed)))
omxCheckTrue(all(got[['school']]$C != student$school$data$observed$C))

omxCheckError(mxGenerateData(student, 10, returnModel=TRUE),
	      paste("Specification of the number of rows",
		    "is not supported for relational models"))

got <- mxGenerateData(student, returnModel=TRUE)
omxCheckTrue(is(got, "MxModel"))
omxCheckTrue(all(got$school$data$observed$C != student$school$data$observed$C))
