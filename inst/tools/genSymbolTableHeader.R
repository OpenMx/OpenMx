table <- read.table('data/omxSymbolTable.tab', header = TRUE)

numEntries <- dim(table)[[1]]

output <- paste(
"#ifndef _OMX_SYMBOL_TABLE_",
"#define _OMX_SYMBOL_TABLE_ TRUE",

"#include <R.h>",
"#include <Rinternals.h>",
"#include <Rdefines.h>",
"#include <R_ext/Rdynload.h>",
"#include <R_ext/BLAS.h>",
"#include <R_ext/Lapack.h>",

"typedef struct {",

"	unsigned int number;",
"	const char opName[250];",
"	const char rName[250];",
"	short int numArgs;",
"	void* funWrapper;",

"} omxAlgebraTableEntry;",

paste("#define omxSymbolTableLength", numEntries),

"const omxAlgebraTableEntry omxAlgebraSymbolTable[omxSymbolTableLength];",


"#endif", sep = "\n")

cat(output)