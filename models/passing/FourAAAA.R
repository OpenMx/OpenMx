A <- mxMatrix(nrow = 1, ncol = 1, labels = 'data.A', name = 'A')
model1 <- mxModel('model1', A)

B <- mxMatrix(nrow = 2, ncol = 1, free = TRUE, labels = c('A', 'data.A'), name = 'B')
data <- mxData(matrix(1, dimnames = list(c(), c('A'))), type = 'raw')
model2 <- mxModel('model2', B, data)

superModel <- mxModel('superModel')
data <- mxData(matrix(0, dimnames = list(c(), c('A'))), type = 'raw')
superModel <- mxModel(superModel, model1, model2, data)

namespace <- omxGenerateNamespace(superModel)
flatModel <- omxFlattenModel(superModel, namespace)
