#pragma once

/*
 * FFRNNS.h
 *
 * Copyright (C) 2020 by University of Tuebingen (BDVA).
 * All rights reserved.
 */

/**
 * This is a CUDA implementation of 'Fast Fixed Radius Nearest Neighbour Search'
 */
#ifndef FFRNNS_H_INCLUDED
#define FFRNNS_H_INCLUDED

#include <glm/glm.hpp>
#include <vector>
#include "vislib/Array.h"
#include "vislib/ArrayAllocator.h"
#include "vislib/math/Vector.h"
#include "vislib/math/Cuboid.h"

#include "FFRNNS.cuh"
#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>


typedef unsigned int uint;
typedef vislib::Array<vislib::math::Vector<float, 4>> visFloat4;

namespace megamol {
namespace prolint {
/***********************************************************
 *** FFRNNS - Fast Fixed Radius Nearest Neighbour Search ***
 ***********************************************************/

class FFRNNS {
public:
	/** Ctor. */
	FFRNNS(void);

	/** Dtor. */
	virtual ~FFRNNS(void);

protected:
	/**
	 * Implementation of 'create'.
	 */
	virtual void create(void);

	/**
	 * Implementation of 'release'.
	 */
	virtual void release(void);

private:
	/**
	 * variables
	 */
	
	// the data points(x,y,z,w->radius) 
    visFloat4 points;

	// the radius around a particle to search for neighbours
    float searchRadius;
	
	// the raw data size (points x 4)
    uint dataSize;

	// the number of data points (used to build search structure)
    uint pointCnt;

	// the number of search data points (used to search neighbours for)
	uint searchPointCnt;

	// the bounding box of the given data
	vislib::math::Cuboid<float> bbox;

	// the size of a grid cell in x,y,z direction
    float3 gridCellSize;

	// the extent of the whole grid (max x,y,z)
	int3 gridDimensions;

	// the minimum values for moving the data to let the min point (lowest x,y,z) start at (0,0,0)
    float3 moveOriginPoint;

	// the neihbour list for each data point
    std::vector<std::vector<uint>> neighbours;

	// bool param which indicates whether grid paramters has changed
	bool gridParamsChanged;

	// bool param which indicates whether the point data has changed
    bool pointDataChanged;

	// number of blocks to use for CUDA kernel call
	int blocks;

	// number of threads per block to use for CUDA kernel call
    int threads;

	// grid  insert:
	int* d_grid;

	/* the size of the gird depending on the chosen gridCellSize
	 * (numberOfCell in X dircetion * numberOfCell Y * numberOfCell Z)
	 */
    int gridSize;

	//**** for GRID INSERT *****
	/* the assignment of data points to grid cells 
	 * d_insertGrid[]:	data point ID
	 * d_insertGrid.x:	grid cell ID
	 * d_insertGrid.y:	intern index of a "data point" in a "grid cell"
	 */
    int2* d_insertGrid;

	// a device copy of pointPos (host array) (x,y,z,w->radius) 
    float4* d_pointPos;
    float4* d_searchPointPos;

	//**** for COUNTING SORT ****
	/* the assignment of data points to grid cells sorted on "grid cell ID"
     * d_sorted_insertGrid[]:	no meaning
     * d_sorted_insertGrid.x:	grid cell ID
     * d_sorted_insertGrid.y:	intern index of a "data point" in a "grid cell"
     */
    int2* d_sorted_insertGrid;

	//**** for REINDEXING ****
    /* the assignment of grid cells to a range of data points in d_sorted_insertGrid
     * d_sorted_insertGrid[]:	grid cell ID
     * d_sorted_insertGrid.x:	start position index
     * d_sorted_insertGrid.y:	end position index
     */
    int2* d_cellStartEnd;

	//**** for NEIGHBOUR SEARCH ****
	// plain neighbour list on (device/host)
    uint* d_neighbours;
    uint* h_neighbours;

    // number of neighbours for each data point ID (device/host)
    uint* d_neighbourCounts;
    uint* h_neighbourCounts;

    // indexSum of the number of neighbours for each data point ID (device/host)
    uint* d_neighbourCountsIndexSum;
    uint* h_neighbourCountsIndexSum;

	// the total number of neighbours
    int totalNeighbourCnt;

	// bool param which indicates whether neighbours are to be calculated
	bool doCalcNeigbours;


/**
 * setter functions
 */

public:
    /**
	 * Sets the point data (x,y,z,w) used to build search structure for neighbour search
	 *
	 * @param	data (size points x 4)
	 * @param	searchRadius (distance for searching neighbours around the point in addtion to its radius)
     * @src		prolint::FFRNNS::setData
	 */
	bool setData(vislib::Array<float> data, float searchRadius);

	/** 
     * Sets the point data (x,y,z,w) used to build search structure for neighbour search
     *
     * @param	data (Array of vector<float,4>)
     * @param	searchRadius (distance for searching neighbours around the point in addtion to its radius)
     * @src		prolint::FFRNNS::setData
     */
    bool setData(visFloat4 data, float searchRadius);

	/**
     * Sets the point data (x,y,z,w) used to build search structure for neighbour search
     *
     * @param	data (Array of vector<float,4>)
     * @param	searchRadius (distance for searching neighbours around the point in addtion to its radius)
     * @param	searchRadius (distance for searching neighbours around the point in addtion to its radius)
     * @src		prolint::FFRNNS::setData
     */
    bool setData(visFloat4 data, float searchRadius, float pointRadius);

	/**
     * Sets the point data (x,y,z) and radius (w) used to build search structure for neighbour search
     *
     * @param	data (size points x 3)
     * @param	searchRadius (distance for searching neighbours around the point in addtion to its radius)
     * @param	pointRadius (fixed radius set for all points)
	 * @src		prolint::FFRNNS::setData
     */
    bool setData(vislib::Array<float> data, float searchRadius, float pointRadius);


	/**
     * optional: sets grid cell size of neighbour search grid (cubic voxels)
     *			 if not set it is automatically determined [ceil(2 x maxRad + 2x searchRadius)]
	 *
     * @param	gridCellSize (used for x,y,z extent of grid cells)
     * @src		prolint::FFRNNS::setGridCellSize
     */
	void setGridCellSize(int gridCellSize);

	/**
     * optional: sets grid cell size of neighbour search grid (also non cubic voxels)
     *			 if not set it is automatically determined [ceil(2 x maxRad +2 x searchRadius)]
     *
     * @param	gridCellSize (used for x,y,z extent of grid cells)
     * @src		prolint::FFRNNS::setGridCellSize
     */
    template <typename...> struct always_false { static constexpr bool value = false; };
    template <typename... Ts>
    void setGridCellSize(int gridCellSize_x, int gridCellSize_y, int gridCellSize_z, Ts&&... );

/**
 * getter functions
 */
public:
// TODO get pointer to data
// TODO get grid Cell Size
// TODO get grid dimensions
// TODO get bbox

	/**
	 * uses the points set with setData() to search neighbours 
	 *
	 * @this-param	this->points
     * @return		this->neighbours (a list)
     * @src			prolint::FFRNNS::getNeigbours
     */
    std::vector<std::vector<uint>>& getNeigbours();

    /**
	 * uses the given searchPoints to search neighbours (from the points set with setData()) 
	 *
     * @param	searchPoints (size points x 4) [(x,y,z,w) * n]
     * @return		this->neighbours (a list)
     * @src		prolint::FFRNNS::getNeigbours (a list)
     */
    std::vector<std::vector<uint>>& getNeigbours(vislib::Array<float> searchPoints);

	/**
	 * uses the given searchPoints to search neighbours (from the points set with setData()) 
	 *
     * @param	searchPoints (Array of vector<float,4>) [(x,y,z,w) * n]
	 * @return	this->neighbours
     * @src		prolint::FFRNNS::getNeigbours (a list)
     */
    std::vector<std::vector<uint>>& getNeigbours(visFloat4 searchPoints);
    

/** 
 * intern functions
 */
//*** private:
private:
	/**
     * calculates the bounding box of the all point data (x,y,z)
     *
     * @this-param	this->points
     * @this-param	this->bbox
     * @src			prolint::FFRNNS::calcBoundingBox
     */
	void calcBoundingBox();

	/**
     * calculates paramters of the grid(dimensions, grid cell size, min values of bbox)
     *
     * @this-param	this->points
     * @this-param	this->bbox
     * @src			prolint::FFRNNS::calcGridParams
     */
	bool calcGridParams();

	/**
     * performs the neighbour search on GPU
     *
     * @param	d_searchPointPos
     * @src		prolint::FFRNNS::calcNeigbours
     */
	bool calcNeigbours(float4 * d_searchPointPos);

	/**
     * CUDA error report (for standard CUDA functions)
     *
     * @param	code: the CUDA error code
     * @param	file: the CUDA file where the error occured
	 * @param   line: the line of the CUDA file where the error occured
     * @src		prolint::FFRNNS::check_CUDA_err
	 *
	 * inspiration:	https://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
     */
	void check_CUDA_err(cudaError_t code, const char* file, int line);

	/**
     * CUDA error report (for custom CUDA kernels/functions)
	 *
	 * @info	just returning the last CUDA ERROR that occured
     * @src		prolint::FFRNNS::check_CUDA_err
     */
    void check_CUDA_err(void);

	/**
     * deletes allocated device and host arrays which only needed during the neighbour seach
	 *
     * @this-param	this->d_grid
     * @this-param	this->d_insertGrid
     * @this-param	this->d_pointPos
     * @this-param	this->d_sorted_insertGrid
	 * @this-param	this->d_cellStartEnd
     * @this-param	this->d_neighbours
	 * @this-param	this->d_neighbourCounts
     * @this-param	this->d_neighbourCountsIndexSum
     * @this-param	this->h_neighbourCounts
     * @this-param	this->h_neighbours
     * @this-param	this->h_neighbourCountsIndexSum
     * @src			prolint::FFRNNS::cleanUp
     */
	void cleanUp(void);

	

};

} // namespace prolint
} // namepsace megamol


#endif // FFRNNS_H_INCLUDED