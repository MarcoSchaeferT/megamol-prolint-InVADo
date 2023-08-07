#pragma once
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>
#include <thrust/device_ptr.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/unique.h>


#define mem_reserve 1.1f
#define max_num_neigbors 150
//#define TEST

//#define reduceAtoms    


// insert atoms into bins and count atoms per bin (atomicAdd)
extern "C" void gridInsertCuda(int blocks, int threads, int* d_grid, int2* d_insert_grid, float4* d_atomPos,
    float3 move_origin_pdb, int3 griddims, int atom_count, float grid_step_size);

extern "C" void gridInsertCuda_intersection(int blocks, int threads, int* d_grid, int2* d_insert_grid,
    float3* d_intersections, float3 move_origin_pdb, int3 griddims, int atom_count, float grid_step_size);

// counting sort on GPU
extern "C" void countingSortCuda(
    int blocks, int threads, int2* d_insert_grid, int2* d_sorted, int atom_count, int* d_pre_sum);

// reindexing on GPU (cell_start, cell_end)
extern "C" void reindexCuda(int blocks, int threads, int2* d_sorted, int2* d_cell_start_end, int atom_count);

// searching neighbours for each atom on GPU (cell_start, cell_end)
extern "C" void NeighbourSearchCuda(int blocks, int threads, float4* d_atomPos, int atom_count, float3 move_origin_pdb,
    int3 griddims, int2* d_cell_start_end, int2* d_sorted, unsigned int* d_neighbours,
    unsigned int* d_reduced_neighbours, unsigned int* d_neighbourCounts, unsigned int* d_reduced_neighbourCounts,
    float probeRad, float grid_step_size);

extern "C" void NeighbourSearch4f_intersectionsCuda(int blocks, int threads, float4* d_atomPos, int filtered_dim,
    float3 move_origin_pdb, int3 griddims, int2* d_cell_start_end, int2* d_sorted, float4* d_neighbours,
    int* d_neighbourCounts, float probeRad, float grid_step_size);

extern "C" void prepare_Number_of_CombinationsCuda(int blocks, int threads, int atom_count,
    unsigned int* d_reduced_neighbours, unsigned int* d_reduced_neighbourCounts, int* d_neigbour_combination_numbers);

extern "C" void check_3Spheres_OverlapCuda(int blocks, int threads, int* atomicIdx, int3* d_SpheresComb3,
    const float4* d_atomPos, int d_calcInter_On_off, unsigned int* d_neighbours, unsigned int* d_neighbourCounts,
    int number_of_3atomCombos, float probeRad, int memSize, float4* d_intersections, int4* d_filtered_SpheresComb3, int* d_atomPosIdx);

extern "C" void calc_3spheresCombosCuda(int blocks, int threads, int atom_count, unsigned int* d_reduced_neighbours,
    unsigned int* d_reduced_neighbourCounts, int* number_of_3atomCombos, int3* d_SpheresComb3);

extern "C" void ScanCuda(int* start, int* end, int* dst);

extern "C" void writeTorusAxesCUDA(
    int blocks, int threads, unsigned int numSphericalTriangles, int4* sphericalTriangles, int2* d_torusAxes);

extern "C" void cutSphereTrianglesCUDA(int blocks, int threads, unsigned int numAtoms, float4* atomPos,
    unsigned int sphereVertCnt, float3* sphereVert, unsigned int sphereTriaCnt, uint3* sphereTria,
    unsigned int* neighbours, unsigned int* neighbourCounts,
    uint3* outerSphereTria, uint3* cutSphereTria);

extern "C" cudaError sortAndUniqueTorusAxesCUDA(unsigned int& numTorusAxes, int2* d_torusAxes);
extern "C" cudaError SortAndUniqueAtomPosIdxCUDA(unsigned int& sizeAtomPosIdx, int* d_torusAxes);

__device__ __host__ float distance2Vectors(float3 v1, float3 v2);


/*
 * greater than comparison for int2 (compares x values)
 */
struct greaterInt2X {
    __host__ __device__ bool operator()(const int2& lhs, const int2& rhs) const {
        return (lhs.x > rhs.x) || ((lhs.x == rhs.x) && (lhs.y > rhs.y));
    }
};

/*
 * less than comparison for int2 (compares x values)
 */
struct lessInt2X {
    __host__ __device__ bool operator()(const int2& lhs, const int2& rhs) const {
        return (lhs.x < rhs.x) || ((lhs.x == rhs.x) && (lhs.y < rhs.y));
    }
};

/*
 * equal_to comparison adapted from thrust for int2
 */
struct equalInt2 {
    __host__ __device__ bool operator()(const int2& lhs, const int2& rhs) const {
        return (lhs.x == rhs.x && lhs.y == rhs.y);
    }
};

/*
 * less than comparison for int2 (compares x values)
 */
struct lessIntX {
    __host__ __device__ bool operator()(const unsigned int& lhs, const unsigned int& rhs) const {
        return (lhs < rhs);
    }
};

/*
 * equal_to comparison adapted from thrust for int2
 */
struct equalInt {
    __host__ __device__ bool operator()(const unsigned int& lhs, const unsigned int& rhs) const {
        return (lhs == rhs);
    }
};
