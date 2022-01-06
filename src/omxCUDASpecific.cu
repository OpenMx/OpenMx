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
#include <iostream>
#include <stdio.h>

#include "cuda_runtime.h"
#include <cusolverDn.h>

#include "omxCUDASpecific.h"

void gpuCholeskyInvertAndDiag(double* h_input, double* h_result, double* h_diag, int N, int* h_devinfo)
{
  // Instantiate solver GPU context
  cusolverDnHandle_t solver_handle;
  cusolverDnCreate(&solver_handle);

  // Allocate memory and copy input array to GPU
  double *d_input;
  cudaMalloc(&d_input, N*N*sizeof(double));
  cudaMemcpy(d_input, h_input, N*N*sizeof(double), cudaMemcpyHostToDevice);

  // Determine block allocations for optimized Cholesky factorization
  int worksize = 0;
  cusolverDnDpotrf_bufferSize(solver_handle, CUBLAS_FILL_MODE_UPPER, N, d_input, N, &worksize);
  double *work;
  cudaMalloc(&work, worksize*sizeof(double));
  int* d_devinfo;
  cudaMalloc(&d_devinfo, sizeof(int));

  // Run Cholesky Factorization
  cusolverDnDpotrf(solver_handle, CUBLAS_FILL_MODE_UPPER, N, d_input, N, work, worksize, d_devinfo);
  cudaMemcpy(h_devinfo, d_devinfo, sizeof(int), cudaMemcpyDeviceToHost);
  //std::cout << *h_devinfo <<std::endl;
  if (*h_devinfo != 0){
    std::cout << "GPU Cholesky Factorization failed!" << std::endl;
    return;
  }
  getDiagonalFromDevice(d_input, h_diag, N);

  // Setup Identity Matrix & Solve For Inversion
  // Inversion result modifies d_identity in place
  double *h_identity = (double*)malloc(N*N*sizeof(double));
  initIdentityMatrix(h_identity, N);
  double *d_identity;
  cudaMalloc(&d_identity, N*N*sizeof(double));
  cudaMemcpy(d_identity, h_identity, N*N*sizeof(double), cudaMemcpyHostToDevice);
  cusolverDnDpotrs(solver_handle, CUBLAS_FILL_MODE_UPPER, N, N, d_input, N, d_identity, N, d_devinfo);
  cudaMemcpy(h_devinfo, d_devinfo, sizeof(int), cudaMemcpyDeviceToHost);
  if (*h_devinfo != 0){
    std::cout << "GPU Solve for Inversion failed!" << std::endl;
    return;
  }
  cudaMemcpy(h_result, d_identity, N*N*sizeof(double), cudaMemcpyDeviceToHost);
  return;
}

void initIdentityMatrix(double* array, int N){
    for (int i=0; i < N*N; i++)
    {
    	if ((i-i/N)%N == 0) array[i] = 1.0;
    	else array[i] = 0.0;
    }
    return;
}

void getDiagonalFromDevice(double* d_array, double* h_diagvec, int N){
	for(int i = 0; i < N; i++){
		cudaMemcpy(&h_diagvec[i], &d_array[i+i*N], sizeof(double), cudaMemcpyDeviceToHost);
	}
  return;
}
