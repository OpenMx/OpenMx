/*
 *  Copyright 2007-2009 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxDefines.h"
#include "omxAlgebraFunctions.h"
#include "omxSymbolTable.h"
#include "omxData.h"

#ifndef _OMX_OBJECTIVE_METADATA_
#define _OMX_OBJECTIVE_METADATA_ TRUE
#define numObjectiveMetadatas 1

// void* omxInitRAMMetaData(SEXP rObj, omxFIMLObjective* ofo, omxState* currentState);

/* Objective Metadata Table */
typedef struct omxObjectiveMetadataContainer {
  omxMatrix *cov, *means;
  void* subObjective;			// Inner Objective Object
	void (*covarianceMeansFunction)(void* subObjective, omxMatrix* cov, omxMatrix* means);
								// Inner Objective Function
	void (*destroySubObjective)(void* subObjective, struct omxObjectiveMetadataContainer* omc);
} omxObjectiveMetadataContainer;

typedef struct omxObjectiveMetadata {
	char name[50];
	 void (*initFunction)(SEXP input, omxObjectiveMetadataContainer* omc, omxState* currentState);
} omxObjectiveMetadata;


const omxObjectiveMetadata omxObjectiveMetadataTable[numObjectiveMetadatas];

#endif /* _OMX_OBJECTIVE_METADATA_ */
