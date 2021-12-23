/*
 *  Copyright 2021 by the individuals mentioned in the source code history
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
#ifndef u_OMX_CUDA_SPECIFIC_H_
#define u_OMX_CUDA_SPECIFIC_H_

void gpuCholeskyInvertAndDiag(double* h_input, double* h_result, double* h_diag, int N, int* h_devinfo);
void initIdentityMatrix(double* array, int N);
void getDiagonalFromDevice(double* d_array, double* h_diagvec, int N);

#endif // #define u_OMX_CUDA_SPECIFIC_H
