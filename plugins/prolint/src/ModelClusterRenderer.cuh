#pragma once
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <thrust/scan.h>
#include <thrust/unique.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>



typedef unsigned int uint;


/*
 * greater than comparison for float 
 */
struct greater_float {
    __device__ bool operator()(const float& lhs, const float& rhs) const;
};


extern "C" cudaError SortTrianglesDevice(
    uint triaCnt, uint3* sortedIndices, float* d_TriangleToCamDistances);

extern "C" cudaError getTriangleToCamDistancesCUDA(int block, int threads, uint triaCnt, float* vertices, uint* indices,  float3* camPos, float* d_TriangleToCamDistances);

__device__ __host__ float distanceTwoVectors(float3 v1, float3 v2);