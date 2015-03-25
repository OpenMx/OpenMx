omxSymbolTable <- read.table("util/omxSymbolTable.tab", stringsAsFactors=FALSE, header=TRUE)
save(omxSymbolTable, file="R/sysdata.rda")
