header <- paste(
	"#include \"omxSymbolTable.h\"",
	"#include \"AlgebraFunctions.h\"",

	"const omxAlgebraTableEntry omxAlgebraSymbolTable[] = {",
	sep = "\n"
)

symbolTable <- 'data/omxSymbolTable.tab'
if (file.access(symbolTable) < 0) {
	symbolTable <- 'data/omxSymbolTable.tab.gz'
} 

table <- read.table(symbolTable, header = TRUE,
	stringsAsFactors = FALSE)

strings <- apply(table, 1, function(x) {
	paste('{', x[[1]], ',\t"', x[[2]], '",\t"', x[[3]],
		'",\t', x[[4]], ",\t",x[[5]], '},', sep = '')
})

strings <- paste(strings, collapse = '\n')

output <- paste(header, strings, '};', sep = '\n')

cat(output)
