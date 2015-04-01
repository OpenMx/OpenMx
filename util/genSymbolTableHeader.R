symbolTable <- 'util/omxSymbolTable.tab'
table <- read.table(symbolTable, header = TRUE)

declares <- apply(table, 1, function(x) {		# Generates a declaration line for each of the functions in the table.
	if(x[[5]] != 'NULL') {
		paste('void ', x[[5]], '(omxMatrix** args, int numArgs, omxMatrix* result);', sep="")
	} else {
		''
	}
})

declares <- paste(declares, collapse = '\n')

output <- paste(
"#ifndef _OMX_SYMBOL_TABLE_",
"#define _OMX_SYMBOL_TABLE_ TRUE",

"typedef struct omxAlgebraTableEntry omxAlgebraTableEntry;",
"#include \"omxMatrix.h\"",
"#include \"algebraOp.h\"",
"struct omxAlgebraTableEntry {",

"	unsigned int number;",
"	const char opName[32];",
"	const char rName[32];",
"	int numArgs;",
"	algebra_op_t funWrapper;",

"};",

declares,

"extern const omxAlgebraTableEntry omxAlgebraSymbolTable[];",
"#endif", sep = "\n")

cat(output)
cat("\n")
