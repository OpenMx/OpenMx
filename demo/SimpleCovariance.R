library(OpenMx)

# Define a model
model <- mxModel()
model <- mxModel(model, mxMatrix("Full", c(0,0.2,0,0), name="A", nrow=2, ncol=2))
model <- mxModel(model, mxMatrix("Full", c(0.8,0,0,0.8), name="S", nrow=2, ncol=2, free=TRUE))
model <- mxModel(model, mxMatrix("Full", c(1,0,0,1), name="F", nrow=2, ncol=2))
model <- mxModel(model, mxBounds(c("apple", "banana"), 0, NA))

model[["A"]]@specification[2,1] <- NA
model[["S"]]@specification[2,1] <- 0
model[["S"]]@specification[1,2] <- 0
model[["S"]]@specification[1,1] <- "apple"
model[["S"]]@specification[2,2] <- "banana"

# Define the objective function
objective <- mxRAMObjective("A", "S", "F", "objective")

# Define the observed covariance matrix
covMatrix <- matrix( c(0.77642931, 0.39590663, 0.39590663, 0.49115615), nrow = 2, ncol = 2, byrow = TRUE)

# Add the objective function and the data to the model
model <- mxModel(model, objective, covMatrix)

# Run the job
model <- mxJobRun(model)

print(model@output)
