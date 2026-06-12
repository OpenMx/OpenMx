library(OpenMx)

mgen <- mxModel('mg', type='RAM', manifestVars = c('a','b'),
              mxPath(c('a','b'), arrows=2, values=1, labels=paste0(c('a','b'),'Var')),
              mxPath('a', 'b', values=.5, labels="reg"))

set.seed(42)
data1 <- mxGenerateData(mgen, nrows = 400)
mgen$S$values['a','a'] <- .5
mgen$S$values['b','b'] <- .5
data2 <- mxGenerateData(mgen, nrows = 200)

m1 <- mxModel(mgen, mxData(data1, 'raw'), mxFitFunctionWLS(), name='g1')
m2 <- mxModel(mgen, mxData(data2, 'raw'), mxFitFunctionWLS(), name='g2')
mg <- mxModel('two', m1, m2, mxFitFunctionMultigroup(paste0('g', 1:2)))

# Turn off Hessian and Standard Errors
mg <- mxOption(mg, "Calculate Hessian", "No")
mg <- mxOption(mg, "Standard Errors", "No")

mg <- mxRun(mg)

message("--- Multigroup WLS (No SE/Hessian) output values ---")
cat("chi: ", mg$output$chi, "\n")
cat("chiDoF: ", mg$output$chiDoF, "\n")
cat("chiM: ", mg$output$chiM, "\n")
cat("chiMV: ", mg$output$chiMV, "\n")
cat("chiMadjust: ", mg$output$chiMadjust, "\n")
cat("chiMVadjust: ", mg$output$chiMVadjust, "\n")
cat("chiDoFstar: ", mg$output$chiDoFstar, "\n")
cat("fitUnits: ", mg$output$fitUnits, "\n")
