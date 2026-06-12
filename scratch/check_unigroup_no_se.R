library(OpenMx)

mgen <- mxModel('mg', type='RAM', manifestVars = c('a','b'),
              mxPath(c('a','b'), arrows=2, values=1, labels=paste0(c('a','b'),'Var')),
              mxPath('a', 'b', values=.5, labels="reg"))

set.seed(42)
data1 <- mxGenerateData(mgen, nrows = 400)
m1 <- mxModel(mgen, mxData(data1, 'raw'), mxFitFunctionWLS(), name='g1')

# Turn off Hessian and Standard Errors
m1 <- mxOption(m1, "Calculate Hessian", "No")
m1 <- mxOption(m1, "Standard Errors", "No")

m1 <- mxRun(m1)

message("--- Single-Group WLS (No SE/Hessian) output values ---")
cat("chi: ", m1$output$chi, "\n")
cat("chiDoF: ", m1$output$chiDoF, "\n")
cat("chiM: ", m1$output$chiM, "\n")
cat("chiMV: ", m1$output$chiMV, "\n")
cat("chiMadjust: ", m1$output$chiMadjust, "\n")
cat("chiMVadjust: ", m1$output$chiMVadjust, "\n")
cat("chiDoFstar: ", m1$output$chiDoFstar, "\n")
cat("fitUnits: ", m1$output$fitUnits, "\n")
