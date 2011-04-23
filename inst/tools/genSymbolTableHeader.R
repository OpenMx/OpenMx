symbolTable <- 'data/omxSymbolTable.tab'
if (file.access(symbolTable) < 0) {
	symbolTable <- 'data/omxSymbolTable.tab.gz'
} 
table <- read.table(symbolTable, header = TRUE)

numEntries <- dim(table)[[1]]

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

"#include <R.h>",
"#include <Rinternals.h>",
"#include <Rdefines.h>",
"#include <R_ext/Rdynload.h>",
"#include <R_ext/BLAS.h>",
"#include <R_ext/Lapack.h>",
"typedef struct omxAlgebraTableEntry omxAlgebraTableEntry;",
"#include \"omxMatrix.h\"",

"struct omxAlgebraTableEntry {",

"	unsigned int number;",
"	const char opName[250];",
"	const char rName[250];",
"	short int numArgs;",
"	void* funWrapper;",

"};",

declares,

paste("#define omxSymbolTableLength", numEntries),

"const omxAlgebraTableEntry omxAlgebraSymbolTable[omxSymbolTableLength];",


"#endif", sep = "\n")

cat(output)
