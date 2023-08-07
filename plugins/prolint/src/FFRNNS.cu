/*
 * FFRNNS.cu
 *
 * Copyright (C) 2020 by University of Tuebingen (BDVA).
 * All rights reserved.
 */

#include "FFRNNS.cuh"
#include <cuda_runtime_api.h>
#include <device_atomic_functions.h>
#include <device_functions.h>
#include <device_launch_parameters.h>
#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/unique.h>


// helper functions
/**
 * pointDistance of FFRNNS
 *
 * @info: calculates the distance/length between to points
 */
__host__ __device__ float pointDistance(float4 v1, float4 v2) {
    float dis = sqrtf(powf(v2.x - v1.x, 2) + powf(v2.y - v1.y, 2) + powf(v2.z - v1.z, 2));
    return dis;
}


/**
 * prolint::countingSortCUDA of FFRNNS
 *
 * @info: insert points into bins and count points per bin (atomicAdd)
 */
__global__ void gridAssignPoints(int* d_grid, int2* d_insert_grid, float4* d_pointPos, float3 move_origin_pdb,
    int3 griddims, int pointCnt, float3 grid_step_size) {

    unsigned int x, y, z, grid_cell = 0;
    int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    int pointID = blockId * blockDim.x + threadIdx.x;

    if (pointID >= pointCnt) return;

    float3 grid_spacing = { 1.0f / grid_step_size.x, 1.0f / grid_step_size.y, 1.0f / grid_step_size.z};
    float4 point = d_pointPos[pointID];

    x = (unsigned int)floorf((point.x - move_origin_pdb.x) * grid_spacing.x);
    y = (unsigned int)floorf((point.y - move_origin_pdb.y) * grid_spacing.y);
    z = (unsigned int)floorf((point.z - move_origin_pdb.z) * grid_spacing.z);

    // get cell_number of grid
    grid_cell = x * (griddims.y * griddims.z) + y * griddims.z + z;

    // add +1 for cell_counter
    int internIndex = atomicAdd(&d_grid[grid_cell], 1);


    d_insert_grid[pointID] = {(int)grid_cell, internIndex};

}


/**
 * prolint::countingSortCUDA of FFRNNS
 *
 * @info: performs a counting sort on GPU
 */
__global__ void countingSort(int2* d_insertGrid, int2* d_sorted_insertGrid, int pointCnt, int* d_pre_sum) {

    int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    int pointID = blockId * blockDim.x + threadIdx.x;

    if (pointID >= pointCnt) return;
    /* d_insert_grid:	 
	*		index = pointID; 
	*		x	  = gridcell
	*		y	  = intern index of point in this cell
	*
	* newIDx = prefixSum[gridCellID of point] - intern index in gridCell - 1
	*/
    int newIDx = d_pre_sum[d_insertGrid[pointID].x] - d_insertGrid[pointID].y - 1;
    d_sorted_insertGrid[newIDx].x = pointID;
    d_sorted_insertGrid[newIDx].y = d_insertGrid[pointID].x;
}


/**
 * prolint::reIndexCUDA of FFRNNS
 *
 * @info: reIndexing on GPU to get cell_start- and cell_end index of "data points" assigned to the grid
 */
__global__ void reIndex(int2* d_sorted_insertGrid, int2* d_cellStartEnd, int pointCnt) {

    int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    int j = blockId * blockDim.x + threadIdx.x;

    if (j >= pointCnt) {
        return;
    }

    if (j > 0) { // last cell index != current cell index  =>(start position of cell)
        if (d_sorted_insertGrid[j - 1].y != d_sorted_insertGrid[j].y) {
            d_cellStartEnd[(d_sorted_insertGrid[j].y)].x = j; // start position
            d_cellStartEnd[(d_sorted_insertGrid[j - 1].y)].y = j; // end position before
        }
    } else {
		// case j = 0 (first entry)
        d_cellStartEnd[(d_sorted_insertGrid[0].y)].x = 0;
		// case j = pointCnt - 1 (last entry)
        d_cellStartEnd[d_sorted_insertGrid[pointCnt - 1].y].y = pointCnt;
    }
}


/**
 * prolint::NeighbourSearchCUDA of FFRNNS
 *
 * @info: performs the neighbour search (currently 9x9 "grid cells" arround the "data point" of interest)
 */
__global__ void nearestNeighbourSearch(float4* d_searchPointPos, float4* d_pointPos, int searchPointCnt,
    float3 moveOriginPoint, int3 gridDimensions, int2* d_cellStartEnd, int2* d_sorted_insertGrid,
    unsigned int* d_neighbours, unsigned int* d_neighbourCounts, unsigned int* d_neighbourCountsIndexSum,
    float searchRadius, float3 gridCellSize, bool saveNeighbours) {

    int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    int i = blockId * blockDim.x + threadIdx.x;

	// don't run out of bounds
    if (i >= searchPointCnt) return;
   
    float3 grid_spacing = {1.0f / gridCellSize.x, 1.0f / gridCellSize.y, 1.0f / gridCellSize.z};
    int atomX, atomY, atomZ, start, end, cell, neighbourCount = 0;

    atomX = ((unsigned int)floorf((d_searchPointPos[i].x - moveOriginPoint.x) * grid_spacing.x)) - 1;
    atomY = ((unsigned int)floorf((d_searchPointPos[i].y - moveOriginPoint.y) * grid_spacing.y)) - 1;
    atomZ = ((unsigned int)floorf((d_searchPointPos[i].z - moveOriginPoint.z) * grid_spacing.z)) - 1;

    for (int x = atomX; x < atomX + 3; x++) {
        for (int y = atomY; y < atomY + 3; y++) {
            for (int z = atomZ; z < atomZ + 3; z++) {
                if (x >= 0 && y >= 0 && z >= 0 && x < gridDimensions.x && y < gridDimensions.y &&
                    z < gridDimensions.z) { // boundaries of grid
                    cell = x * (gridDimensions.y * gridDimensions.z) + y * gridDimensions.z + z;
                    start = d_cellStartEnd[cell].x;
                    end = d_cellStartEnd[cell].y;
                    for (int k = start; k < end; k++) {
                        int atom_pos_index = d_sorted_insertGrid[k].x;
						// TODO: should it still be allowed: self containing in neighbour list?
                        //if (i != atom_pos_index) {
                            // check for radius
                            float detectionRadius =
                                d_searchPointPos[i].w + d_pointPos[atom_pos_index].w + searchRadius - 0.001f;
                        float atomDistance = pointDistance(d_searchPointPos[i], d_pointPos[atom_pos_index]);
                            if (atomDistance < detectionRadius) {
                                if (saveNeighbours) {   
									d_neighbours[d_neighbourCountsIndexSum[i] - neighbourCount] =
                                        atom_pos_index;
								}
								neighbourCount++;
                            }
                       // }
                    }
                }
            }
        }
    }
    if (saveNeighbours) {
		d_neighbourCounts[i] = neighbourCount;
    } else {
        d_neighbourCountsIndexSum[i] = neighbourCount;
    }
    
}



namespace megamol {
namespace prolint {

/*************************************
****** EXTERN C DEFINITIONS **********
**************************************/

extern "C" void gridAssignPointsCUDA(int blocks, int threads, int* d_grid, int2* d_insert_grid, float4* d_pointPos,
    float3 move_origin_pdb, int3 griddims, int pointCnt, float3 grid_step_size) {

    // insert points into bins and count points per bin (atomicAdd)
    gridAssignPoints<<<blocks, threads>>>(
        d_grid, d_insert_grid, d_pointPos, move_origin_pdb, griddims, pointCnt, grid_step_size);
}


extern "C" void prefixSumCUDA(int* start, int* end, int* dst) {
    /**
     * prolint::prefixSumCUDA of FFRNNS
     *
     * @param	start: start index to use
     * @param	end: end index to use
     * @param	dst: destination for storing the prefix sum
     */
    thrust::inclusive_scan(
        thrust::device_ptr<int>(start), thrust::device_ptr<int>(end), thrust::device_ptr<int>(dst)); // in-place scan
}


extern "C" void countingSortCUDA(
    int blocks, int threads, int2* d_insertGrid, int2* d_sorted_insertGrid, int pointCnt, int* d_pre_sum) {

    countingSort<<<blocks, threads>>>(d_insertGrid, d_sorted_insertGrid, pointCnt, d_pre_sum);
}


extern "C" void reIndexCUDA(int blocks, int threads, int2* d_sorted_insertGrid, int2* d_cellStartEnd, int pointCnt) {

    reIndex<<<blocks, threads>>>(d_sorted_insertGrid, d_cellStartEnd, pointCnt);
}

extern "C" void NeighbourSearchCUDA(int blocks, int threads, float4* d_searchPointPos, float4* d_pointPos,
    int searchPointCnt, float3 moveOriginPoint, int3 gridDimensions, int2* d_cellStartEnd, int2* d_sorted_insertGrid,
    unsigned int* d_neighbours, unsigned int* d_neighbourCounts, unsigned int* d_neighbourCountsIndexSum,
    float searchRadius, float3 gridCellSize, bool saveNeighbours) {

	nearestNeighbourSearch<<<blocks, threads>>>(d_searchPointPos, d_pointPos, searchPointCnt, moveOriginPoint,
        gridDimensions, d_cellStartEnd, d_sorted_insertGrid, d_neighbours, d_neighbourCounts, d_neighbourCountsIndexSum,
        searchRadius, gridCellSize, saveNeighbours);

}

	
} // namespace prolint
} // namespace megamol
