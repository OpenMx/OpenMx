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

#include "omxState.h"

/*
 * I would prefer this global variable didn't exist at all,
 * except that the NPSOL API appears to be non thread-safe. See
 * the function "F77_SUB(objectiveFunction)" to convince
 * yourself of this fact.
 *
 * It may be possible to create multiple copies of "objectiveFunction"
 * in memory, but even then we are unsure if the NPSOL internals
 * are thread-safe.
 *
 * So for now, we can't perform any shared memory
 * parallel computations that invoke the optimizer. -mspiegel
 */

omxState* globalState;			// Current State of optimization
