/*
 *  Copyright 2007-2013 The OpenMx Project
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

#include <stdarg.h>
#include <errno.h>

#include "omxState.h"
#include "Compute.h"
#include "omxOpenmpWrap.h"

struct omxGlobal *Global = NULL;

/* Initialize and Destroy */
	void omxInitState(omxState* state) {
		state->numConstraints = 0;
		state->childList = NULL;
		state->conList = NULL;

		state->majorIteration = 0;
		state->minorIteration = 0;
		state->startTime = 0;
		state->endTime = 0;
		state->numCheckpoints = 0;
		state->checkpointList = NULL;
		state->chkptText1 = NULL;
		state->chkptText2 = NULL;

		state->computeCount = 0;
		state->currentRow = -1;

		strncpy(state->statusMsg, "", 1);
	}

	void omxSetMajorIteration(omxState *state, int value) {
		state->majorIteration = value;
		if (!state->childList) return;
		for(int i = 0; i < Global->numChildren; i++) {
			omxSetMajorIteration(state->childList[i], value);
		}
	}

	void omxSetMinorIteration(omxState *state, int value) {
		state->minorIteration = value;
		if (!state->childList) return;
		for(int i = 0; i < Global->numChildren; i++) {
			omxSetMinorIteration(state->childList[i], value);
		}
	}
	
	void omxDuplicateState(omxState* tgt, omxState* src) {
		tgt->dataList			= src->dataList;
		
		for(size_t mx = 0; mx < src->matrixList.size(); mx++) {
			// TODO: Smarter inference for which matrices to duplicate
			tgt->matrixList.push_back(omxDuplicateMatrix(src->matrixList[mx], tgt));
		}

		tgt->numConstraints     = src->numConstraints;
		tgt->conList			= (omxConstraint*) R_alloc(tgt->numConstraints, sizeof(omxConstraint));
		for(int j = 0; j < tgt->numConstraints; j++) {
			tgt->conList[j].size   = src->conList[j].size;
			tgt->conList[j].opCode = src->conList[j].opCode;
			tgt->conList[j].lbound = src->conList[j].lbound;
			tgt->conList[j].ubound = src->conList[j].ubound;
			tgt->conList[j].result = omxDuplicateMatrix(src->conList[j].result, tgt);
		}

		for(size_t j = 0; j < src->algebraList.size(); j++) {
			// TODO: Smarter inference for which algebras to duplicate
			tgt->algebraList.push_back(omxDuplicateMatrix(src->algebraList[j], tgt));
		}

		for(size_t j = 0; j < src->expectationList.size(); j++) {
			// TODO: Smarter inference for which expectations to duplicate
			tgt->expectationList.push_back(omxDuplicateExpectation(src->expectationList[j], tgt));
		}

		for(size_t j = 0; j < tgt->algebraList.size(); j++) {
			omxDuplicateAlgebra(tgt->algebraList[j], src->algebraList[j], tgt);
		}

		for(size_t j = 0; j < src->expectationList.size(); j++) {
			// TODO: Smarter inference for which expectations to duplicate
			omxCompleteExpectation(tgt->expectationList[j]);
		}

		tgt->childList 			= NULL;

		tgt->majorIteration 	= 0;
		tgt->minorIteration 	= 0;
		tgt->startTime 			= src->startTime;
		tgt->endTime			= 0;
		
		// TODO: adjust checkpointing based on parallelization method
		tgt->numCheckpoints     = 0;
		tgt->checkpointList 	= NULL;
		tgt->chkptText1 		= NULL;
		tgt->chkptText2 		= NULL;
                                  
		tgt->computeCount 		= src->computeCount;
		tgt->currentRow 		= src->currentRow;

		strncpy(tgt->statusMsg, "", 1);
	}

	omxMatrix* omxLookupDuplicateElement(omxState* os, omxMatrix* element) {
		if(os == NULL || element == NULL) return NULL;

		if (element->hasMatrixNumber) {
			int matrixNumber = element->matrixNumber;
			if (matrixNumber >= 0) {
				return(os->algebraList[matrixNumber]);
			} else {
				return(os->matrixList[-matrixNumber - 1]);
			}
		}

		omxConstraint* parentConList = globalState->conList;

		for(int i = 0; i < os->numConstraints; i++) {
			if(parentConList[i].result == element) {
				if(os->conList[i].result != NULL) {   // Not sure of proper failure behavior here.
            	return(os->conList[i].result);
				} else {
                    omxRaiseError(os, -2, "Initialization Copy Error: Constraint required but not yet processed.");
            }
			}
		}

		return NULL;
	}
	
void omxFreeChildStates(omxState *state)
{
	if (!state->childList || Global->numChildren == 0) return;

	for(int k = 0; k < Global->numChildren; k++) {
		// Data are not modified and not copied. The same memory
		// is shared across all instances of state. We only need
		// to free the data once, so let the parent do it.
		state->childList[k]->dataList.clear();

		omxFreeState(state->childList[k]);
	}
	Free(state->childList);
	state->childList = NULL;
	Global->numChildren = 0;
}

	void omxFreeState(omxState *state) {
		omxFreeChildStates(state);

		for(size_t ax = 0; ax < state->algebraList.size(); ax++) {
			if(OMX_DEBUG) { mxLog("Freeing Algebra %d at 0x%x.", ax, state->algebraList[ax]); }
			omxFreeAllMatrixData(state->algebraList[ax]);
		}

		if(OMX_DEBUG) { mxLog("Freeing %d Matrices.", state->matrixList.size());}
		for(size_t mk = 0; mk < state->matrixList.size(); mk++) {
			if(OMX_DEBUG) { mxLog("Freeing Matrix %d at 0x%x.", mk, state->matrixList[mk]); }
			omxFreeAllMatrixData(state->matrixList[mk]);
		}
		
		if(OMX_DEBUG) { mxLog("Freeing %d Model Expectations.", state->expectationList.size());}
		for(size_t ex = 0; ex < state->expectationList.size(); ex++) {
			if(OMX_DEBUG) { mxLog("Freeing Expectation %d at 0x%x.", ex, state->expectationList[ex]); }
			omxFreeExpectationArgs(state->expectationList[ex]);
		}

		if(OMX_DEBUG) { mxLog("Freeing %d Constraints.", state->numConstraints);}
		for(int k = 0; k < state->numConstraints; k++) {
			if(OMX_DEBUG) { mxLog("Freeing Constraint %d at 0x%x.", k, state->conList[k]); }
			omxFreeAllMatrixData(state->conList[k].result);
		}

		if(OMX_DEBUG) { mxLog("Freeing %d Data Sets.", state->dataList.size());}
		for(size_t dx = 0; dx < state->dataList.size(); dx++) {
			if(OMX_DEBUG) { mxLog("Freeing Data Set %d at 0x%x.", dx, state->dataList[dx]); }
			omxFreeData(state->dataList[dx]);
		}

		if(OMX_DEBUG) { mxLog("Freeing %d Checkpoints.", state->numCheckpoints);}
		for(int k = 0; k < state->numCheckpoints; k++) {
			if(OMX_DEBUG) { mxLog("Freeing Data Set %d at 0x%x.", k, state->checkpointList[k]); }
			omxCheckpoint oC = state->checkpointList[k];
			switch(oC.type) {
				case OMX_FILE_CHECKPOINT:
					fclose(oC.file);
					break;
				case OMX_CONNECTION_CHECKPOINT:	// NYI :::DEBUG:::
					// Do nothing: this should be handled by R upon return.
					break;
			}
			if(state->chkptText1 != NULL) {
				Free(state->chkptText1);
			}
			if(state->chkptText2 != NULL) {
				Free(state->chkptText2);
			}
			// Checkpoint list itself is freed by R.
		}

		delete state;

		if(OMX_DEBUG) { mxLog("State Freed.");}
	}

	void omxResetStatus(omxState *state) {
		int numChildren = Global->numChildren;
		state->statusMsg[0] = '\0';
		if (!state->childList) return;
		for(int i = 0; i < numChildren; i++) {
			omxResetStatus(state->childList[i]);
		}
	}

std::string string_snprintf(const std::string fmt, ...)
{
    int size = 100;
    std::string str;
    va_list ap;
    while (1) {
        str.resize(size);
        va_start(ap, fmt);
        int n = vsnprintf((char *)str.c_str(), size, fmt.c_str(), ap);
        va_end(ap);
        if (n > -1 && n < size) {
            str.resize(n);
            return str;
        }
        if (n > -1)
            size = n + 1;
        else
            size *= 2;
    }
    return str;
}

void mxLogBig(const std::string str)   // thread-safe
{
	ssize_t len = ssize_t(str.size());
	ssize_t wrote = 0;
#pragma omp critical(stderr)
	while (1) {
		ssize_t got = write(2, str.data() + wrote, len - wrote);
		if (got == EINTR) continue;
		if (got <= 0) error("mxLogBig failed with errno=%d", got);
		wrote += got;
		if (wrote == len) break;
	}
}

void mxLog(const char* msg, ...)   // thread-safe
{
	const int maxLen = 240;
	char buf1[maxLen];
	char buf2[maxLen];

	va_list ap;
	va_start(ap, msg);
	vsnprintf(buf1, maxLen, msg, ap);
	va_end(ap);

	int len = snprintf(buf2, maxLen, "[%d] %s\n", omx_absolute_thread_num(), buf1);

	ssize_t wrote = 0;
#pragma omp critical(stderr)
	while (1) {
		ssize_t got = write(2, buf2 + wrote, len - wrote);
		if (got == EINTR) continue;
		if (got <= 0) error("mxLog failed with errno=%d", got);
		wrote += got;
		if (wrote == len) break;
	}
}

void _omxRaiseError()
{
	// keep for debugger breakpoints
}

void omxRaiseErrorf(omxState *state, const char* errorMsg, ...)
{
	_omxRaiseError();
	va_list ap;
	va_start(ap, errorMsg);
	int fit = vsnprintf(state->statusMsg, MAX_STRING_LEN, errorMsg, ap);
	va_end(ap);
	if(OMX_DEBUG) {
		if (!(fit > -1 && fit < MAX_STRING_LEN)) {
			mxLog("Error exceeded maximum length: %s", errorMsg);
		} else {
			mxLog("Error raised: %s", state->statusMsg);
		}
	}
}

	void omxRaiseError(omxState *state, int errorCode, const char* errorMsg) { // DEPRECATED
		_omxRaiseError();
		if(OMX_DEBUG && errorCode) { mxLog("Error %d raised: %s", errorCode, errorMsg);}
		if(OMX_DEBUG && !errorCode) { mxLog("Error status cleared."); }
		strncpy(state->statusMsg, errorMsg, 249);
		state->statusMsg[249] = '\0';
	}

	void omxStateNextRow(omxState *state) {
		state->currentRow++;
	};

	void omxStateNextEvaluation(omxState *state) {
		state->currentRow = -1;
		state->computeCount++;
	};

static void omxWriteCheckpointHeader(omxState *os, omxCheckpoint* oC) {
	// rewrite with std::string TODO
	std::vector< omxFreeVar* > &vars = Global->freeGroup[0]->vars;
	size_t numParam = vars.size();

		os->chkptText1 = (char*) Calloc((24 + 15 * numParam), char);
		os->chkptText2 = (char*) Calloc(1.0 + 15.0 * numParam*
			(numParam + 1.0) / 2.0, char);
		if (oC->type == OMX_FILE_CHECKPOINT) {
			fprintf(oC->file, "iterations\ttimestamp\tobjective\t");
			for(size_t j = 0; j < numParam; j++) {
				if(strcmp(vars[j]->name, CHAR(NA_STRING)) == 0) {
					fprintf(oC->file, "%s", vars[j]->name);
				} else {
					fprintf(oC->file, "\"%s\"", vars[j]->name);
				}
				if (j != numParam - 1) fprintf(oC->file, "\t");
			}
			fprintf(oC->file, "\n");
			fflush(oC->file);
		}
	}
 
void omxWriteCheckpointMessage(char *msg) {
	std::vector< omxFreeVar* > &vars = Global->freeGroup[0]->vars;
	size_t numParam = vars.size();

		omxState *os = globalState;
		for(int i = 0; i < os->numCheckpoints; i++) {
			omxCheckpoint* oC = &(os->checkpointList[i]);
			if(os->chkptText1 == NULL) {    // First one: set up output
				omxWriteCheckpointHeader(os, oC);
			}
			if (oC->type == OMX_FILE_CHECKPOINT) {
				fprintf(oC->file, "%d \"%s\" NA ", os->majorIteration, msg);
				for(size_t j = 0; j < numParam; j++) {
					fprintf(oC->file, "NA ");
				}
				fprintf(oC->file, "\n");
			}
		}
	}

void omxSaveCheckpoint(double* x, double f, int force) {
	// rewrite with std::string TODO
	std::vector< omxFreeVar* > &vars = Global->freeGroup[0]->vars;
	size_t numParam = vars.size();

		omxState *os = globalState;
		time_t now = time(NULL);
		int soFar = now - os->startTime;		// Translated into minutes
		int n;
		for(int i = 0; i < os->numCheckpoints; i++) {
			n = 0;
			omxCheckpoint* oC = &(os->checkpointList[i]);
			// Check based on time            
			if((oC->time > 0 && (soFar - oC->lastCheckpoint) >= oC->time) || force) {
				oC->lastCheckpoint = soFar;
				n = 1;
			}
			// Or iterations
			if((oC->numIterations > 0 && (os->majorIteration - oC->lastCheckpoint) >= oC->numIterations) || force) {
				oC->lastCheckpoint = os->majorIteration;
				n = 1;
			}

			if(n) {		//In either case, save a checkpoint.
				if(os->chkptText1 == NULL) {	// First one: set up output
					omxWriteCheckpointHeader(os, oC);
				}
				char tempstring[25];
				sprintf(tempstring, "%d", os->majorIteration);

				if(strncmp(os->chkptText1, tempstring, strlen(tempstring))) {	// Returns zero if they're the same.
					struct tm * nowTime = localtime(&now);						// So this only happens if the text is out of date.
					strftime(tempstring, 25, "%b %d %Y %I:%M:%S %p", nowTime);
					sprintf(os->chkptText1, "%d \"%s\" %9.5f", os->majorIteration, tempstring, f);
					for(size_t j = 0; j < numParam; j++) {
						sprintf(tempstring, " %9.5f", x[j]);
						strncat(os->chkptText1, tempstring, 14);
					}
				}

				if(oC->type == OMX_FILE_CHECKPOINT) {
					fprintf(oC->file, "%s", os->chkptText1);
					if(oC->saveHessian)
						fprintf(oC->file, "%s", os->chkptText2);
					fprintf(oC->file, "\n");
					fflush(oC->file);
				} else if(oC->type == OMX_CONNECTION_CHECKPOINT) {
					warning("NYI: R_connections are not yet implemented.");
					oC->numIterations = 0;
					oC->time = 0;
				}
			}
		}
	}

void omxExamineFitOutput(omxState *state, omxMatrix *fitMatrix, int *mode)
{
	if (!R_FINITE(fitMatrix->data[0])) {
		omxRaiseErrorf(state, "Fit function returned %g at iteration %d.%d",
			       fitMatrix->data[0], state->majorIteration, state->minorIteration);
		*mode = -1;
	}
}
