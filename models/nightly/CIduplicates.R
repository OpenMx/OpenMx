library(OpenMx)
data(myFADataRaw, package="OpenMx")
manifests = paste0("x",1:3)
myFADataRaw = myFADataRaw[, manifests]
latents = c("G")
m1 <- mxModel("m1", type="RAM",
	manifestVars = manifests, latentVars   = latents,
	mxPath(from = latents, to = manifests),
	mxPath(from = manifests, arrows = 2, labels = paste0(manifests, "_resid")),
	mxPath(from = latents, arrows = 2, free = F, values = 1), # latents fixed at 1
	mxData(cov(myFADataRaw, use="complete"), type = "cov", numObs = nrow(myFADataRaw))
)
m1 = mxRun(m1)

tmp = mxRun(mxModel(m1, mxCI("x1_resid")), intervals = T)
# Add the same CI with a bracket name: no complaint
tmp = mxRun(mxModel(tmp, mxCI("S[1,1]")), intervals = T)
omxCheckCloseEnough(nrow(tmp$output$confidenceIntervals), 1)

omxCheckCloseEnough(tmp$output$confidenceIntervals[1, c('lbound', 'ubound')],
                    c(.280, .430), .01)

q()

# Add and run 1 CI im each of two models. then concatenate the two model's CIs
tmp1 = mxRun(mxModel(m1, mxCI("x1_resid")), intervals = T)
tmp2 = mxRun(mxModel(m1, mxCI("x1_resid")), intervals = T)
a = tmp1$output$confidenceIntervals
b = tmp2$output$confidenceIntervals
a_names = attr(a, "dimnames")[[1]]
b_names = attr(b, "dimnames")[[1]]
all_names = c(a_names, b_names)
all_CIs = rbind(a, b)
tmp1@output$confidenceIntervals = all_CIs
summary(tmp1)
# confidence intervals:
#      lbound  estimate    ubound note
# 1 0.2802408 0.3525381 0.4308769
# 2 0.2802408 0.3525381 0.4308769
# 3 0.2802408 0.3525381 0.4308769
# 4 0.2802408 0.3525381 0.4308769

# Warning message:
# In data.row.names(row.names, rowsi, i) :
#   some row.names duplicated: 3,4 --> row.names NOT used
