allinputs   <- Sys.getenv("R_SWIFT_ARGS")
model_obj <- noquote(strsplit(allinputs," ")[[1]][1])

library(OpenMx)
library(R.oo)
model = mxModel()
model = load(model_obj)
model = saveLoadReference
model = mxRun(model)
print(model_obj)
save.Object(model,file=paste("results/",model_obj,sep=""))
