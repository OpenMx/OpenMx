library(OpenMx)

numIndicators <- 4

genData <- function(upper, fanout, variation, keyname) {
	lowerData <- NULL
	for (sx in 1:nrow(upper)) {
		extraFanout <- 0L
		if (variation) {
			extraFanout <- sample.int(fanout, 1)
		}
		lowerData <- rbind(lowerData, data.frame(
		    upper=upper[sx,1], skill=rnorm(fanout + extraFanout,
					   mean=upper[sx, 'skill'])))
	}
	colnames(lowerData)[[1]] <- colnames(upper)[[1]]
	lowerData[[keyname]] <- 1:nrow(lowerData)
	lowerData <- lowerData[,c(3,1,2)]
	lowerData
}

createIndicators <- function(latentSkill, indicatorVariance) {
	if (missing(indicatorVariance)) {
		indicatorVariance <- rep(1, numIndicators)
					#rlnorm(numIndicators) / 8
	}
	ind <- matrix(NA, length(latentSkill), length(indicatorVariance))
	for (ix in 1:length(latentSkill)) {
		ind[ix,] <-
		  sapply(indicatorVariance,
			 function(sd) rnorm(1, mean=latentSkill[ix], sd=sd))
	}
	# per indicator mean
#	ind <- t(t(ind) + runif(numIndicators,min=-1,max=1))
	colnames(ind) <- paste0('i', 1:length(indicatorVariance))
	as.data.frame(ind)
}

mkSingleFactor <- function(prefix, ...) {
	oneLevel <- mxModel(prefix, type='RAM', ...,
		manifestVars = paste0('i', 1:numIndicators),
		latentVars = c("skill"),
		mxPath(from='skill', arrows=2, labels="Var",
		       values=rlnorm(1), lbound=.01),
		mxPath(from=paste0('i',1:numIndicators), arrows=2,
		       values=rlnorm(1), labels="Err", lbound=.01),
		mxPath(from="one", to=paste0('i',1:numIndicators),
		       free=TRUE, values=rnorm(4)),
		mxPath(from='skill', to=paste0('i',1:numIndicators),
		       labels=paste0('L',1:numIndicators), lbound=.1,
		       values=c(1, runif(numIndicators-1, .5,1.5)),
		       free=c(FALSE, rep(TRUE,numIndicators-1)))
		)
	for (mat in c("A","S")) {
		lab <- oneLevel[[mat]]$labels
		lab[!is.na(lab)] <- paste0(prefix, lab[!is.na(lab)])
		oneLevel[[mat]]$labels <- lab
	}
	oneLevel$fitfunction <- NULL
	oneLevel
}

buildModel <- function(numDistricts, numSchools, numTeachers, numStudents, v, naprob) {
	districtData <- data.frame(districtID=1:numDistricts,
				   skill=rnorm(numDistricts))
	schoolData <- genData(districtData, numSchools, v[1], 'schoolID')
	teacherData <- genData(schoolData, numTeachers, v[2], 'teacherID')
	studentData <- genData(teacherData, numStudents, v[3], 'studentID')

	districtData <- cbind(districtData, createIndicators(districtData$skill))
	schoolData <- cbind(schoolData, createIndicators(schoolData$skill))
	teacherData <- cbind(teacherData, createIndicators(teacherData$skill))
	studentData <- cbind(studentData, createIndicators(studentData$skill))

	studentData$i4[runif(nrow(studentData)) > 1-naprob] <- NA
	#teacherData$i4[runif(nrow(teacherData)) > .8] <- NA

	dMod <- mkSingleFactor(
		"district",
		mxData(type="raw", observed=districtData,
		       primaryKey="districtID"))

	schMod <- mkSingleFactor(
		"school", dMod,
		mxData(type="raw", observed=schoolData,
		       primaryKey="schoolID"),
		mxPath(from='district.skill', to='skill',
		       joinKey="districtID", values=runif(1)))

	tMod <- mkSingleFactor(
		"teacher", schMod,
		mxData(type="raw", observed=teacherData,
		       primaryKey="teacherID"),
		mxPath(from='school.skill', to='skill',
		       joinKey="schoolID", values=runif(1)))

	sMod <- mkSingleFactor(
		"student", tMod,
		mxData(type="raw", observed=studentData,
		       primaryKey="studentID"),
		mxPath(from='teacher.skill', to='skill',
		       joinKey="teacherID", values=runif(1)))

	sMod$fitfunction <- mxFitFunctionML()

	sMod
}

checkSinglePoint <- function(origModel, cycleStart, perUnit=FALSE) {
	plan <- mxComputeSequence(list(
		mxComputeOnce('fitfunction', 'fit'),
		mxComputeReportExpectation()))

	mod <- mxModel(origModel, plan)
	fit <- list()
	maxRampart <- 3
	layoutLen <- 2L
	for (rampart in rev(cycleStart:maxRampart)) {
		if (!perUnit) {
			mod$expectation$.rampartCycleLimit <- rampart
			fit1 <- mxRun(mod, silent=TRUE)
			print(fit1$expectation$debug$rampartUsage)
			#print(c(rampart, limit, fit1$output$fit))
			fit[[1 + length(fit)]] <- fit1
		} else {
			limit <- 1L
			while (limit < layoutLen) {
				mod$expectation$.rampartCycleLimit <- rampart
				mod$expectation$.rampartUnitLimit <- limit
				fit1 <- mxRun(mod, silent=TRUE)
				layoutLen <- nrow(fit1$expectation$debug$layout)
				#print(fit1$expectation$debug$rampartUsage)
				print(c(rampart, limit, fit1$output$fit))
				fit[[1 + length(fit)]] <- fit1
				limit <- limit + 1L
				if (rampart == 0L) break
			}
		}
	}
	fitVec <- sapply(fit, function(m) m$output$fit)
	print(fitVec)
	fitDiff = diff(fitVec)
	omxCheckCloseEnough(fitDiff, rep(0,length(fitDiff)), 1e-8)
}

set.seed(1)

sMod <- buildModel(5,4,3,5, rep(TRUE,3), .2)

checkSinglePoint(sMod, 0)

fit1 <- mxRun(mxModel(sMod, mxComputeSequence(list(
  mxComputeGradientDescent(),
  mxComputeReportExpectation()))))

print(summary(fit1))

omxCheckCloseEnough(fit1$output$fit, 17144.43, .01)
#We've check it before. If fit is good then the gradient would be good.
#omxCheckCloseEnough(max(abs(fit1$output$gradient)), 0, .01)
ed <- fit1$expectation$debug
omxCheckCloseEnough(ed$rampartUsage, c(902, 16))
omxCheckCloseEnough(ed$numGroups, 14L)
omxCheckCloseEnough(
    sapply(sprintf("g%02d", 1:14),
	   function(x) nrow(ed[[x]]$layout) %/% ed[[x]]$clumpSize),
    c(97L, 805L, 2L, 2L, 2L, 4L, 3L, 1L, 2L, 1L, 1L, 1L,  1L, 1L))

plan <- mxComputeSequence(list(
    mxComputeOnce('fitfunction', 'fit'),
    mxComputeReportExpectation()
))
slow <- fit1
slow$expectation$.rampartCycleLimit <- 0L
slowEx <- mxRun(mxModel(slow, plan))
ed <- slowEx$expectation$debug
omxCheckTrue(length(ed$rampartUsage)==0)
# each (entire) district is an independent unit
omxCheckCloseEnough(
    sapply(unique(ed$layout$group),
	   function(x) nrow(ed$layout[ed$layout$group==x,]) %/% ed[[paste0('g',x)]]$clumpSize),
		    rep(1L,5))
omxCheckCloseEnough(fit1$output$fit - slowEx$output$fit, 0, 1e-8)

if (0) {
	options(width=120)
	plan <- mxComputeSequence(list(
	    mxComputeOnce('fitfunction', 'fit'),
	    mxComputeNumericDeriv(checkGradient=FALSE,
				  hessian=FALSE, iterations=2),
	    mxComputeReportDeriv(),
	    mxComputeReportExpectation()
	))

	sMod$expectation$.rampartCycleLimit <- 0L
	square <- mxRun(mxModel(sMod, plan))

	sMod$expectation$.rampartCycleLimit <- 2L
	rotated <- mxRun(mxModel(sMod, plan))
	
	ex <- square$expectation
	ex <- rotated$expectation
	eo <- ex$output
	ed <- ex$debug
	print(ed$layout)
	print(ed$rampartUsage)
	print(ed$numGroups)
	table(ed$layout$group)
	head(ed$layout[ed$layout$group == 1, ], n=20)
	#print(round(ed$A[1:20,1:20],2))
	#print(round(ed$rA[1:20,1:20],2))
	#print(ed$mean)

	#omxCheckCloseEnough(ed$rampartUsage, c(11064L, 317L, 198L, 2L), 1L)
	print(abs(rotated$output$fit - square$output$fit))
	print(max(abs(rotated$output$gradient - square$output$gradient)))
#	omxCheckCloseEnough(rotated$output$gradient,
#             square$output$gradient, 1e-4)
}

# -------------------------------------------------

# ulimit -m 8388608
# ulimit -v 8388608

set.seed(1)
trivial <- buildModel(4,2,2,2, c(FALSE,FALSE,FALSE), 0)
checkSinglePoint(trivial, 0, TRUE)

set.seed(1)
bigMod <- buildModel(35,4,4,4, rep(TRUE,3), 0)
checkSinglePoint(bigMod, 1)

if (0) {
	plan <- mxComputeSequence(list(
	    mxComputeOnce('fitfunction', 'fit'),
	    mxComputeReportExpectation()
	))
	bigMod <- mxModel(bigMod, plan)

	bigMod$expectation$.useSufficientSets <- FALSE
	bigMod$expectation$.rampartUnitLimit <- 38L  # correct
	fit1 <- mxRun(bigMod)
	bigMod$expectation$.useSufficientSets <- TRUE
	bigMod$expectation$.rampartUnitLimit <- 38L
	fit2 <- mxRun(bigMod)

	fit1$output$fit - fit2$output$fit

	ed1 <- fit1$expectation$debug
	ed2 <- fit2$expectation$debug
	l1 <- ed1$layout
	l2 <- ed2$layout

	ed1$g1$fit - ed2$g1$fit #OK
	ed1$g2$fit - ed2$g2$fit #bad
	ed1$g3$fit - ed2$g3$fit #OK
	ed1$g4$fit - ed2$g4$fit #OK
	ed1$g5$fit - ed2$g5$fit #OK

	str(ed1$g1)
	str(ed2$g1)
	ed1$g1$layout
	ed2$g1$layout

	length(ed2$g2$dataVec)
	ed2$g2$ss1
	rowMeans(matrix(ed2$g2$dataVec[8 + 1:16], nrow=8))
	ed2$g2$ss2
	rowMeans(matrix(ed2$g2$dataVec[8 + 17:32], nrow=8))
}
