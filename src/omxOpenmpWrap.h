/*
 *  Copyright 2007-2011 The OpenMx Project
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

#ifndef _OMX_OPENMP_WRAP_H
#define _OMX_OPENMP_WRAP_H

#include "Rconfig.h"
#include "omxDefines.h"
#include "omxSadmvnWrapper.h"
#include "omxRObjective.h"

#ifdef _OPENMP

#include <omp.h>

static OMXINLINE int omx_omp_get_thread_num(void) {
   return(omp_get_thread_num());
}

static OMXINLINE void omx_omp_init() {
   omp_init_lock(&sadmvn_lock);
   omp_init_lock(&robjective_lock);
}

static OMXINLINE void omx_omp_set_lock(omp_lock_t* lock) {
   omp_set_lock(lock);
}

static OMXINLINE void omx_omp_unset_lock(omp_lock_t* lock) {
   omp_unset_lock(lock);
}

#else

static OMXINLINE int omx_omp_get_thread_num(void) {
   return(0);
}

static OMXINLINE void omx_omp_init() {}

static OMXINLINE void omx_omp_set_lock(void** __attribute__((unused)) lock) {}

static OMXINLINE void omx_omp_unset_lock(void** __attribute__((unused)) lock) {}


#endif // #ifdef _OPENMP




#endif // #ifndef _OMX_OPENMP_WRAP_H
