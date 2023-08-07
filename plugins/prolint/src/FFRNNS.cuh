#pragma once

/*
 * FFRNNS.cuh
 *
 * Copyright (C) 2020 by University of Tuebingen (BDVA).
 * All rights reserved.
 */

/**
 * This is a CUDA implementation of 'Fast Fixed Radius Nearest Neighbour Search'
 */
#ifndef FFRNNS_CUH_INCLUDED
#define FFRNNS_CUH_INCLUDED

#define max_num_neigbors 2500

#include <cuda_runtime.h>

// helper functions
/**
 * pointDistance of FFRNNS
 *
 * @info: calculates the distance/length between to points
 */
float pointDistance(float4 v1, float4 v2);

namespace megamol {
namespace prolint {

/**
 * prolint::countingSortCUDA of FFRNNS
 *
 * @info: insert points into bins and count points per bin (atomicAdd)
 */
extern "C" void gridAssignPointsCUDA(int blocks, int threads, int* d_grid, int2* d_insert_grid, float4* d_pointPos,
    float3 move_origin_pdb, int3 griddims, int pointCnt, float3 grid_step_size);


/**
 * prolint::prefixSumCUDA of FFRNNS
 *
 * @param	start: start index to use
 * @param	end: end index to use
 * @param	dst: destination for storing the prefix sum
 */
extern "C" void prefixSumCUDA(int* start, int* end, int* dst);


/**
 * prolint::countingSortCUDA of FFRNNS
 *
 * @info: performs a counting sort on GPU
 */
extern "C" void countingSortCUDA(
    int blocks, int threads, int2* d_insertGrid, int2* d_sorted_insertGrid, int pointCnt, int* d_pre_sum);


/**
 * prolint::reIndexCUDA of FFRNNS
 *
 * @info: reIndexing on GPU to get cell_start- and cell_end index of "data points" assigned to the grid
 */
extern "C" void reIndexCUDA(int blocks, int threads, int2* d_sorted_insertGrid, int2* d_cellStartEnd, int pointCnt);


/**
 * prolint::NeighbourSearchCUDA of FFRNNS
 *
 * @info: performs the neighbour search (currently 9x9 "grid cells" arround the "data point" of interest)
 */
extern "C" void NeighbourSearchCUDA(int blocks, int threads, float4* d_searchPointPos, float4* d_pointPos,
    int searchPointCnt, float3 moveOriginPoint, int3 gridDimensions, int2* d_cellStartEnd, int2* d_sorted_insertGrid,
    unsigned int* d_neighbours, unsigned int* d_neighbourCounts, unsigned int* d_neighbourCountsIndexSum,
    float searchRadius, float3 gridCellSize, bool saveNeighbours = false);

} // namespace prolint
} // namespace megamol

#endif // FFRNNS_CUH_INCLUDED
