library(OpenMx)

pdat <- data.frame(length=c(1,2), ID=1L:2L, parentID=1L:2L)

m1 <- mxModel("slug", type="RAM",
              manifestVars = c("length"),
              mxData(type="raw", observed=pdat, primaryKey = 'ID'),
              mxPath("one", "length"),
              mxPath("length", arrows=2, values=1),
              mxPath("slug.length", "length", joinKey='parentID', labels="fromParent"))

omxCheckError(mxRun(m1), "slug.expectation cycle detected: 'slug.data' row 1 joins against itself")

