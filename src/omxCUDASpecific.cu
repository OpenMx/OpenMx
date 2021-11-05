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
#include "omxCUDASpecific.h"
#include <stdio.h>

// Squares each cell value in an array of integers
__global__ void square_array(int *array, int arrayCount)
{
 int idx = threadIdx.x + blockIdx.x * blockDim.x;
 if (idx < arrayCount)
 {
   array[idx] *= array[idx];
 }
}

void callHelloWorld(){
  printf("Hello, World from .cu!\n");

  int n = 1024;

  int *array;
  size_t bytes = n*sizeof(int);
  array = (int*)malloc(bytes);

  for (int i =0; i <n; i++){
    array[i] = i;
  }

  int *d_a;

  cudaMalloc(&d_a, n*sizeof(int));

  int blockSize;
  int minGridSize;
  int gridSize;

  cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize,
                                      square_array);
  // Round up according to array size
  gridSize = (n + blockSize - 1) / blockSize;

  cudaMemcpy(d_a, array, n*sizeof(int), cudaMemcpyHostToDevice);

  square_array<<< gridSize, blockSize >>>(d_a, n);

  cudaDeviceSynchronize();

  cudaMemcpy(array, d_a, n*sizeof(int), cudaMemcpyDeviceToHost);

  for (int i =0; i <n; i+=100){
    printf("Square of %d is %d\n", i, array[i]);
  }
}
