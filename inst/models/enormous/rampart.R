library(OpenMx)
library(mvtnorm)

set.seed(3)

numIndicators <- 5

numSchools <- 7
numTeachers <- 3
numStudents <- 5

genStructure <- function(upper, fanout, keyname) {
	lowerData <- NULL
	for (sx in 1:nrow(upper)) {
		extraFanout <- sample.int(fanout, 1)
		lowerData <- rbind(lowerData, data.frame(
		    upper=upper[sx,1], skill=rnorm(fanout + extraFanout, mean=upper[sx, 'skill'])))
	}
	colnames(lowerData)[[1]] <- colnames(upper)[[1]]
	lowerData[[keyname]] <- 1:nrow(lowerData)
	lowerData <- lowerData[,c(3,1,2)]
	lowerData
}

dataEnv <- new.env()

assign("schoolData", data.frame(schoolID=1:numSchools, skill=rnorm(numSchools)), envir=dataEnv)
assign("teacherData", genStructure(dataEnv$schoolData, numTeachers, 'teacherID'), envir=dataEnv)
assign("studentData", genStructure(dataEnv$teacherData, numStudents, 'studentID'), envir=dataEnv)

createIndicators <- function(latentSkill, indicatorMean, indicatorVariance) {
    if (missing(indicatorMean)) {
        indicatorMean <- runif(numIndicators,min=-1,max=1)
    }
    if (missing(indicatorVariance)) {
        indicatorVariance <- rlnorm(numIndicators) / 8
    }
    ind <- matrix(NA, length(latentSkill), length(indicatorVariance))
    for (ix in 1:length(latentSkill)) {
        ind[ix,] <- sapply(indicatorVariance,
                           function(sd) rnorm(1, mean=latentSkill[ix], sd=sd))
    }
    ind <- t(t(ind) + indicatorMean)
    colnames(ind) <- paste0('i', 1:length(indicatorVariance))
    as.data.frame(ind)
}

for (tbl in paste0(c('school', 'teacher', 'student'), 'Data')) {
    dataEnv[[tbl]] <- cbind(dataEnv[[tbl]], createIndicators(dataEnv[[tbl]]$skill))
}

dataEnv$studentData$i1[runif(nrow(dataEnv$studentData)) > .8] <- NA
#teacherData$i4[runif(nrow(teacherData)) > .8] <- NA

mkSingleFactor <- function(latent=c()) {
	mxModel('template', type='RAM',
		manifestVars = paste0('i', 1:numIndicators),
		latentVars = c("skill",latent),
		mxPath('skill', arrows=2, labels="Var", values=rlnorm(1), lbound=.01),
		mxPath(paste0('i',1:numIndicators), arrows=2,
		       values=rlnorm(1), labels="Err", lbound=.01),
		mxPath("one", paste0('i',1:numIndicators), free=TRUE, values=rnorm(4)),
		mxPath('skill', paste0('i',1:numIndicators),
		       labels=paste0('L',1:numIndicators), lbound=0,
		       values=c(1, runif(numIndicators-1, .5,1.5)),
		       free=c(FALSE, rep(TRUE,numIndicators-1)))
		)
}

singleFactor <- mkSingleFactor(NULL)
			
relabel <- function(m, prefix) {
  for (mat in c("A","S")) {
    lab <- m[[mat]]$labels
    lab[!is.na(lab)] <- paste0(prefix, lab[!is.na(lab)])
    m[[mat]]$labels <- lab
  }
  mxModel(m, name=prefix)
}

schMod <- mxModel(relabel(mkSingleFactor(), "school"),
		  mxData(type="raw", observed=dataEnv$schoolData, primaryKey="schoolID", sort=FALSE))

tMod <- mxModel(relabel(singleFactor, "teacher"), schMod,
		  mxData(type="raw", observed=dataEnv$teacherData, primaryKey="teacherID", sort=FALSE),
		  mxPath('school.skill', 'skill', joinKey="schoolID", values=runif(1)))

sMod <- mxModel(relabel(singleFactor, "student"), tMod,
		  mxData(type="raw", observed=dataEnv$studentData, primaryKey="studentID", sort=FALSE),
		  mxPath('teacher.skill', 'skill', joinKey="teacherID", values=runif(1)))

interest <- c('wallTime', 'infoDefinite', 'conditionNumber', 'fit', 'timestamp')

if (1) {
    result <- expand.grid(rampart=c(TRUE,FALSE), rep=1:200, gradient=NA)
    for (e1 in names(coef(sMod))) result[[e1]] <- NA
    for (i1 in interest) result[[i1]] <- NA
} else {
    load("/tmp/rampart.rda")
}

plan <- mxComputeSequence(list(
    mxComputeGradientDescent(),
    mxComputeNumericDeriv(iterations=2L),
    mxComputeHessianQuality(),
    mxComputeReportDeriv()
))

for (rrow in 1:nrow(result)) {
    if (!is.na(result[rrow, 'wallTime'])) next
#    if (!result[rrow, 'rampart']) next

    if (result[rrow, 'rampart']==FALSE &&
        !result[result$rep == result[rrow, 'rep'] & result$rampart==TRUE, 'infoDefinite']) {
        print("skip")
        next
    }

    set.seed(result[rrow, 'rep'])
    trial <- mxGenerateData(sMod, returnModel=TRUE)

    if (result[rrow, 'rampart']) {
        trial$expectation$.rampart <- as.integer(NA)
    } else {
        trial$expectation$.rampart <- 0L
        trial$fitfunction$parallel <- TRUE
    }
    trialFit <- mxRun(mxModel(trial, plan))

    result[rrow, names(coef(trialFit))] <- coef(trialFit)
    result[rrow, interest] <- trialFit$output[interest]
    result[rrow, 'gradient'] <- max(abs(trialFit$output$gradient))

    save(result, file="/tmp/rampart.rda")
}

sum(!is.na(result[result$rampart==TRUE, 'conditionNumber']))
sum(!is.na(result[result$rampart==FALSE, 'conditionNumber']))

cnMask <- (result$conditionNumber < median(result$conditionNumber, na.rm=TRUE) + 5 * mad(result$conditionNumber, na.rm=TRUE))
bothOkay <- cnMask[result$rampart==TRUE] & cnMask[result$rampart==FALSE]
length(which(bothOkay))

good <- result[result$rep %in% which(bothOkay),]
good[,c("rep",'rampart', "conditionNumber", 'gradient')]
cor(good[good$rampart==TRUE,"conditionNumber"], good[good$rampart==FALSE,"conditionNumber"])
cor(good[good$rampart==TRUE,"fit"], good[good$rampart==FALSE,"fit"])
summary <- c(rMean=norm(colMeans(good[good$rampart==TRUE, names(coef(sMod))]) - coef(sMod), "2"),
             fMean=norm(colMeans(good[good$rampart==FALSE, names(coef(sMod))]) - coef(sMod), "2"),
             rVar=norm(apply(good[good$rampart==TRUE, names(coef(sMod))], 2, var), "2"),
             fVar=norm(apply(good[good$rampart==FALSE, names(coef(sMod))], 2, var), "2"))
print(summary)

library(ggplot2)
ggplot(good) + geom_histogram(aes(wallTime)) + facet_wrap(~rampart, scales="free_x")

