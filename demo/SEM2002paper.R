
manifest <- c("x1", "x2", "y1", "y2", "y3")
latent <- c("L")

model <- mxModel(manifestVars = manifest, latentVars = latent)

model <- mxModel(model, 
	mxPath(from = manifest, arrows = 2, free = TRUE))

model <- mxModel(model, 
	mxPath(from = latent, arrows = 2, free = TRUE))
	
model <- mxModel(model,
	mxPath(from = "x1", to = "x2", arrows = 2, free = TRUE))

model <- mxModel(model,
	mxPath(from = c("x1", "x2"), to = "L", arrows = 1, free = TRUE))

model <- mxModel(model,
	mxPath(from = "L", to = c("y1","y2","y3"), 
		arrows = 1, free = TRUE))
