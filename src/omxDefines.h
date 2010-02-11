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
 *
 */

/***********************************************************
 * 
 *  omxDefines.h
 *
 *  Created: Timothy R. Brick 	Date: 2009-09-23
 *
 *	Contains #define information for debugging purposes.
 *
 **********************************************************/
#ifndef _OMXDEFINES_H_
#define _OMXDEFINES_H_

#define OMX_DEFAULT_MAX_PTS 100*cov->rows*cov->rows
#define MAX_STRING_LEN 250
#define EPSILON 1e-16
#define TRUE 1
#define FALSE 0

#ifdef DEBUGMX
#define OMX_DEBUG 1
#define OMX_VERBOSE 1
#else
#ifdef VERBOSEMX
#define OMX_DEBUG 0
#define OMX_VERBOSE 1
#else
#define OMX_DEBUG 0
#define OMX_VERBOSE 0
#endif /* VERBOSEMX */
#endif /* DEBUGMX */

#ifdef DEBUGMX_ROWS
#define OMX_DEBUG_ROWS (row > 3 || row % 10 == 0)
#else
#define OMX_DEBUG_ROWS 0
#endif /* DEBUGMX_ROWS */

#ifdef DEBUGNPSOL
#define OMX_DEBUG_OPTIMIZER 1
#else
#define OMX_DEBUG_OPTIMIZER 0
#endif /* DEBUGNPSOL */

#ifdef DEBUGMX_MATRIX
#define OMX_DEBUG_MATRIX 1
#else
#define OMX_DEBUG_MATRIX 0
#endif /* DEBUGMX_MATRIX */

#ifdef DEBUGMX_ALGEBRA
#define OMX_DEBUG_ALGEBRA 1
#else
#define OMX_DEBUG_ALGEBRA 0
#endif /* DEBUGMX_ALGEBRA */

#endif /* _OMXDEFINES_H_ */