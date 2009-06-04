header <- paste(
	"#include <R.h>",
	"#include <Rinternals.h>",
	"#include <Rdefines.h>",
	"#include <R_ext/Rdynload.h>",
	"#include <R_ext/BLAS.h>",
	"#include <R_ext/Lapack.h>",
	"#include \"omxSymbolTable.h\"",
#	"#include \"omxAlgebraFunctions.h\"",

	"const omxAlgebraTableEntry omxAlgebraSymbolTable[omxSymbolTableLength] = {", 
	sep = "\n"
)

table <- read.table('data/omxSymbolTable.tab', header = TRUE,
	stringsAsFactors = FALSE)

strings <- apply(table, 1, function(x) {
	paste('{', x[[1]], ',\t"', x[[2]], '",\t"', x[[3]],
		'",\t', x[[4]], ",\t(void*)",x[[5]], '},', sep = '')
})

strings <- paste(strings, collapse = '\n')

output <- paste(header, strings, '};', sep = '\n')

cat(output)
