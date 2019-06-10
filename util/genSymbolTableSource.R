header <- paste(
	"#include \"omxSymbolTable.h\"",
	"#include \"AlgebraFunctions.h\"",

	"const omxAlgebraTableEntry omxAlgebraSymbolTable[] = {",
	sep = "\n"
)

symbolTable <- 'util/omxSymbolTable.tab'

table <- read.table(symbolTable, header = TRUE,
	stringsAsFactors = FALSE)

strings <- apply(table, 1, function(x) {
	checkFun <- 'NULL'
	calcFun <- 'NULL'
	if (x[[5]] != 'NULL') {
		checkFun <- paste0(x[[5]],'Check')
		calcFun <- x[[5]]
	}
	paste('{', x[[1]], ',\t"', x[[2]], '",\t"', x[[3]],
		'",\t', x[[4]], ",\t",checkFun, ",\t", calcFun, '},', sep = '')
})

strings <- paste(strings, collapse = '\n')

output <- paste(header, strings, '};', sep = '\n')

cat(c(output, '\n'))
