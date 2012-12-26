/*
 *  Copyright 2007-2012 The OpenMx Project
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
#include "omxRFitFunction.h"

#ifdef _OPENMP

#include <omp.h>

#if _OPENMP <= 200505

static OMXINLINE int omx_absolute_thread_num(void) {
   return(omp_get_thread_num());
}

#else

static OMXINLINE int omx_absolute_thread_num(void) {
   int retval = 0;
   int level = omp_get_level();
   int scale = 1;
   for(int i = level; i > 0; i--) {
       retval += scale * omp_get_ancestor_thread_num(i);
       scale *= omp_get_team_size(i);
   }
   return(retval);
}

#endif

static OMXINLINE void omx_omp_init() {
   omp_init_lock(&rfitfunction_lock);
#if _OPENMP <= 200505
   omp_set_nested(0);
#endif
}

static OMXINLINE void omx_omp_set_lock(omp_lock_t* lock) {
   omp_set_lock(lock);
}

static OMXINLINE void omx_omp_unset_lock(omp_lock_t* lock) {
   omp_unset_lock(lock);
}

#else

static OMXINLINE int omx_absolute_thread_num(void) {
   return(0);
}

static OMXINLINE void omx_omp_init() {}

static OMXINLINE void omx_omp_set_lock(void** __attribute__((unused)) lock) {}

static OMXINLINE void omx_omp_unset_lock(void** __attribute__((unused)) lock) {}


#endif // #ifdef _OPENMP




#endif // #ifndef _OMX_OPENMP_WRAP_H
