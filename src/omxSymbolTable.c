/*
 *  Copyright 2007-2009 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include <R.h> 
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxSymbolTable.h"
	
const omxAlgebraTableEntry omxAlgebraSymbolTable[14] = {
	{0, 	"*SPECIAL*", 				"*NONE*", 	0, 		(void*)NULL},
	{1, 	"Inversion",				"solve", 	1, 		(void*)&omxMatrixInvert       },
	{2, 	"Transposition", 			"t",		1, 		(void*)&omxMatrixTranspose    },
	{3, 	"Element powering",			"%^%", 		2, 		(void*)&omxElementPower       },
	{4, 	"Multiplication", 			"%*%", 		2, 		(void*)&omxMatrixMult         },
	{5, 	"Dot product", 				"*", 		2, 		(void*)&omxMatrixDot          },
	{6, 	"Kronecker product",		"%x%",		2, 		(void*)&omxKroneckerProd      },
	{7, 	"Quadratic product",		"%&%",		2, 		(void*)&omxQuadraticProd      },
	{8, 	"Element division",			"/",		2, 		(void*)&omxElementDivide      },
	{9, 	"Addition",					"+",		2, 		(void*)&omxMatrixAdd          },
	{10,	"Subtraction (binary)",		"-",		2, 		(void*)&omxMatrixSubtract     },
	{11,	"Subtraction (unary)",		"-",		1, 		(void*)&omxUnaryMinus         },
	{12,	"Horizontal adhesion",		"cbind",	-1,		(void*)&omxMatrixHorizCat     },
	{13,	"Vertical adhesion",		"rbind",	-1,		(void*)&omxMatrixVertCat      },
};
