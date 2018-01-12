library(OpenMx)

omxCheckEquals(capture.output(mxPath("a","b")),
               c('mxPath',
                 "a -> b [value=0, free=TRUE]"))
               
omxCheckEquals(capture.output(mxPath("a","b",lbound=0)),
               c('mxPath',
                 "a -> b [value=0, free=TRUE, lbound=0]"))
                              
omxCheckEquals(capture.output(mxPath("model.a","b",ubound=0)),
               c("mxPath", "model.a -> b [value=0, free=TRUE, ubound=0]"))

omxCheckEquals(capture.output(mxPath("a","b", free=FALSE)),
               c("mxPath", "a -> b [value=0, free=FALSE]"))

omxCheckEquals(capture.output(mxPath("G", paste0("i",1:5))),
               c("mxPath", "G -> i1 [value=0, free=TRUE]", "G -> i2 [value=0, free=TRUE]",  "G -> i3 [value=0, free=TRUE]", "G -> i4 [value=0, free=TRUE]",  "G -> i5 [value=0, free=TRUE]"))

omxCheckEquals(capture.output(mxPath("G", arrows=2, values=1)),
               c("mxPath", "G <-> G [value=1, free=TRUE]"))

omxCheckEquals(capture.output(mxPath(paste0('i',1:3), connect='all.pairs')),
               c("mxPath", "i1 -> i1 [value=0, free=TRUE]", "i1 -> i2 [value=0, free=TRUE]",  "i1 -> i3 [value=0, free=TRUE]", "i2 -> i1 [value=0, free=TRUE]",  "i2 -> i2 [value=0, free=TRUE]", "i2 -> i3 [value=0, free=TRUE]",  "i3 -> i1 [value=0, free=TRUE]", "i3 -> i2 [value=0, free=TRUE]",  "i3 -> i3 [value=0, free=TRUE]"))

omxCheckEquals(capture.output(mxPath(paste0('i',1:3), connect='unique.pairs')),
               c(c("mxPath", "i1 -> i1 [value=0, free=TRUE]", "i1 -> i2 [value=0, free=TRUE]",  "i1 -> i3 [value=0, free=TRUE]", "i2 -> i2 [value=0, free=TRUE]",  "i2 -> i3 [value=0, free=TRUE]", "i3 -> i3 [value=0, free=TRUE]" )))

omxCheckEquals(capture.output(mxPath(paste0('i',1:3), connect='all.bivariate')),
               c("mxPath", "i1 -> i2 [value=0, free=TRUE]", "i1 -> i3 [value=0, free=TRUE]",  "i2 -> i1 [value=0, free=TRUE]", "i2 -> i3 [value=0, free=TRUE]",  "i3 -> i1 [value=0, free=TRUE]", "i3 -> i2 [value=0, free=TRUE]" ))

omxCheckEquals(capture.output(mxPath(paste0('i',1:3), connect='unique.bivariate')),
               c("mxPath", "i1 -> i2 [value=0, free=TRUE]", "i1 -> i3 [value=0, free=TRUE]",  "i2 -> i3 [value=0, free=TRUE]"))

