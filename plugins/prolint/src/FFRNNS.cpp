/*
 * FFRNNS.cpp
 *
 * Copyright (C) 2020 by University of Tuebingen (BDVA).
 * All rights reserved.
 */

#include "FFRNNS.h"
#include "FFRNNS.cuh"
#include <algorithm>    // std::sort only used for testing results

#define NUM_THREADS 256


namespace megamol {
namespace prolint {

/**
 * prolint::FFRNNS::FFRNNS (CTOR)
 */
FFRNNS::FFRNNS(void) { this->create(); }


/**
 * prolint::FFRNNS::~FFRNNS (DTOR)
 */
FFRNNS::~FFRNNS(void) { this->release(); }


/**
 * prolint::FFRNNS::create
 */
void FFRNNS::create(void) {
    this->points.Resize(0);
    this->searchRadius = 0.f;
    this->dataSize = 0;
    this->pointCnt = 0;
    this->bbox.Set(0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
    this->gridCellSize = {0.f, 0.f, 0.f};
    this->gridDimensions = {0, 0, 0};
    this->moveOriginPoint = {0.f, 0.f, 0.f};
    this->gridParamsChanged = true;
    this->blocks = 0;
    this->threads = NUM_THREADS;
    this->pointDataChanged = true;
    this->doCalcNeigbours = true;
};


/**
 * prolint::FFRNNS::release
 */
void FFRNNS::release(void) {
	// delete copy of point data
    this->points.Resize(0);
	// free memory of neighbor vectors
    if (this->neighbours.size() > 0) {
		for (int i = 0; i < this->searchPointCnt; i++) {
			this->neighbours[i].resize(0);
		}
    }
    this->neighbours.resize(0);

};


/**
 * Sets the point data (x,y,z,w) used to build search structure for neighbour search
 *
 * @param	data (size points x 4)
 * @param	searchRadius (distance for searching neighbours around the point in addtion to its radius)
 * @src		prolint::FFRNNS::setData
 */
bool FFRNNS::setData(vislib::Array<float> data, float searchRadius) {
    this->pointDataChanged = true;
    // set everthing to zero/default
    create();
    this->dataSize = data.Count();

    if (this->dataSize % 4 != 0) {
        printf("ERROR in FFRNNS::setData: given data is not a multiple of 4 (x,y,z,w)!\n");
        return false;
    }

    this->pointCnt = this->dataSize / 4;
    this->points.SetCount(this->pointCnt);

#pragma omp parallel for
    for (int i = 0; i < this->pointCnt; i++) {
        this->points[i].SetX(data[i * 4 + 0]);
        this->points[i].SetY(data[i * 4 + 1]);
        this->points[i].SetZ(data[i * 4 + 2]);
        this->points[i].SetW(data[i * 4 + 3]);
    }
    
    this->searchRadius = searchRadius;
    calcBoundingBox();
    return true;
}




/**
 * Sets the point data (x,y,z,w) used to build search structure for neighbour search
 *
 * @param	data (Array of vector<float,4>)
 * @param	searchRadius (distance for searching neighbours around the point in addtion to its radius)
 * @src		prolint::FFRNNS::setData
 */
bool FFRNNS::setData(visFloat4 data, float searchRadius) {
    this->pointDataChanged = true;
    // set everthing to zero/default
    create();
    this->dataSize = data.Count() * 4;
    this->pointCnt = data.Count();
    this->points.SetCount(this->pointCnt);

#pragma omp parallel for
    for (int i = 0; i < this->pointCnt; i++) {
        points[i].SetX(data[i].GetX());
        points[i].SetY(data[i].GetY());
        points[i].SetZ(data[i].GetZ());
        points[i].SetW(data[i].GetW());
    }
    this->searchRadius = searchRadius;
    calcBoundingBox();
    return true;
}


/**
 * Sets the point data (x,y,z,w) used to build search structure for neighbour search
 *
 * @param	data (Array of vector<float,4>)
 * @param	searchRadius (distance for searching neighbours around the point in addtion to its radius)
 * @param	searchRadius (distance for searching neighbours around the point in addtion to its radius)
 * @src		prolint::FFRNNS::setData
 */
bool FFRNNS::setData(visFloat4 data, float searchRadius, float pointRadius) {
    this->pointDataChanged = true;
    // set everthing to zero/default
    create();
    this->dataSize = data.Count() * 4;
    this->pointCnt = data.Count();
    this->points.SetCount(this->pointCnt);

#pragma omp parallel for
    for (int i = 0; i < this->pointCnt; i++) {
        points[i].SetX(data[i].GetX());
        points[i].SetY(data[i].GetY());
        points[i].SetZ(data[i].GetZ());
        points[i].SetW(pointRadius);
    }
    this->searchRadius = searchRadius;
    calcBoundingBox();
    return true;
}


/**
 * Sets the point data (x,y,z) and radius (w) used to build search structure for neighbour search
 *
 * @param	data (size points x 3)
 * @param	pointRadius (fixed radius set for all points)
 * @param	searchRadius (distance for searching neighbours around the point in addtion to its radius)
 * @src		prolint::FFRNNS::setData
 */
bool FFRNNS::setData(vislib::Array<float> data, float searchRadius, float pointRadius) {
    this->pointDataChanged = true;
    // set everthing to zero/default
    create();
    this->dataSize = data.Count() / 3 * 4;

    if (data.Count() % 3 != 0) {
        printf("ERROR in FFRNNS::setData: given data is not a multiple of 3 (x,y,z)!\n");
        return false;
    }

    this->pointCnt = this->dataSize / 4;
    this->points.SetCount(this->pointCnt);

#pragma omp parallel for
    for (int i = 0; i < this->pointCnt; i++) {
        points[i].SetX(data[i * 4 + 0]);
        points[i].SetY(data[i * 4 + 1]);
        points[i].SetZ(data[i * 4 + 2]);
        points[i].SetW(pointRadius);
    }
    this->searchRadius = searchRadius;
    calcBoundingBox();
    return true;
    getNeigbours();
}


/**
 * uses the points set with setData() to search neighbours
 *
 * @this-param	this->points
 * @return		this->neighbours (a pointer to list)
 * @src			prolint::FFRNNS::getNeigbours
 */
std::vector<std::vector<uint>>& FFRNNS::getNeigbours() { 
	cudaMalloc((void**)&this->d_pointPos, this->pointCnt * sizeof(visFloat4));
    
	// cuda memcopy not needed: d_pointPos is filled with data in calcNeigbours()
	//cudaMemcpy(this->d_pointPos, &this->points, this->pointCnt * sizeof(visFloat4), cudaMemcpyHostToDevice);

	this->d_searchPointPos = this->d_pointPos;
    this->searchPointCnt = this->pointCnt;
	calcNeigbours(this->d_searchPointPos);

    return this->neighbours;
	
}

/**
 * uses the given searchPoints to search neighbours (from the points set with setData())
 *
 * @param	searchPoints (size points x 4) [(x,y,z,w) * n]
 * @return		this->neighbours (a list)
 * @src		prolint::FFRNNS::getNeigbours (a list)
 */
std::vector<std::vector<uint>>& FFRNNS::getNeigbours(vislib::Array<float> searchPoints) {
    cudaMalloc((void**)&this->d_pointPos, this->pointCnt * sizeof(visFloat4));
    if (searchPoints.Count() % 4 != 0) {
        printf("ERROR in FFRNNS::getNeigbours: given data is not a multiple of 4 (x,y,z,w)!\n");
        std::vector<std::vector<uint>> a;
        return a;
    }
    this->searchPointCnt = searchPoints.Count() / 4;
    cudaMalloc((void**)&this->d_searchPointPos, this->searchPointCnt * sizeof(float4));
    cudaMemcpy(this->d_searchPointPos, &searchPoints[0], this->searchPointCnt * sizeof(float4), cudaMemcpyHostToDevice);
   
	calcNeigbours(this->d_searchPointPos);
   
	return this->neighbours;
}

/**
 * uses the given searchPoints to search neighbours (from the points set with setData())
 *
 * @param	searchPoints (Array of vector<float,4>) [(x,y,z,w) * n]
 * @return	this->neighbours
 * @src		prolint::FFRNNS::getNeigbours (a list)
 */
std::vector<std::vector<uint>>& FFRNNS::getNeigbours(visFloat4 searchPoints) {
    cudaMalloc((void**)&this->d_pointPos, this->pointCnt * sizeof(visFloat4));

    this->searchPointCnt = searchPoints.Count();
	cudaMalloc((void**)&this->d_searchPointPos, this->searchPointCnt * sizeof(float4));
    cudaMemcpy(this->d_searchPointPos, &searchPoints[0][0], this->searchPointCnt * sizeof(float4), cudaMemcpyHostToDevice);

	calcNeigbours(this->d_searchPointPos);

    return this->neighbours;
}

/**
 * performs the neighbour search on GPU
 *
 * @param	d_searchPointPos
 * @src		prolint::FFRNNS::calcNeigbours
 */
bool FFRNNS::calcNeigbours(float4* d_searchPointPos) {

    if (this->pointCnt == 0) {
        printf("ERROR: no data in FFRNNS!\n");
        return false;
    }

    // get grid paramters
    calcGridParams();
    if (this->gridParamsChanged || this->pointDataChanged) {
        this->doCalcNeigbours = true;
    }

    if (this->doCalcNeigbours) {

        
        // determine number of blocks and threads for CUDA kernesl
        this->threads = threads < this->pointCnt ? threads : this->pointCnt;
        this->blocks = (this->pointCnt % threads != 0) ? (this->pointCnt / threads + 1) : (this->pointCnt / threads);

		

        /*****************************************************************
        ************ GRID INSERT - ATOMS INTO BINS AND COUNT *************
        ******************************************************************/
		cudaMalloc((void**)&this->d_grid, this->gridSize * sizeof(int));
		cudaMalloc((void**)&this->d_insertGrid, this->pointCnt * sizeof(int2));
		
       
		cudaMemcpy(this->d_pointPos, &this->points[0][0], this->pointCnt * sizeof(float4), cudaMemcpyHostToDevice);
        // needed to be set to 0, because using atomicAdd()
        cudaMemset(this->d_grid, 0, this->gridSize * sizeof(int));
        cudaMemset(this->d_insertGrid, 0, this->pointCnt * sizeof(int2));
        // CUDA: insert atoms into bins and count atoms per bin (atomicAdd)
        gridAssignPointsCUDA(this->blocks, this->threads, this->d_grid, this->d_insertGrid, this->d_pointPos,
            this->moveOriginPoint, this->gridDimensions, this->pointCnt, this->gridCellSize);


        /***************************************
        ************ COUNTING SORT *************
        ****************************************/
        // calculate Prefix Sum
        prefixSumCUDA(this->d_grid, this->d_grid + this->gridSize, this->d_grid); // in-place scan
        cudaMalloc((void**)&this->d_sorted_insertGrid, (this->pointCnt) * sizeof(int2));
        // counting sort on GPU
        countingSortCUDA(this->blocks, this->threads, this->d_insertGrid, this->d_sorted_insertGrid, this->pointCnt, d_grid);

        /************************************
        ************ REINDEXING *************
        *************************************/
        cudaMalloc((void**)&this->d_cellStartEnd, (this->gridSize) * sizeof(int2));
        // memset needed for empty cells (set start and end to zero)
        cudaMemset(this->d_cellStartEnd, 0, this->gridSize * sizeof(int2));
        // CUDA reindexing on GPU (cell_start, cell_end)
        reIndexCUDA(this->blocks, this->threads, this->d_sorted_insertGrid, this->d_cellStartEnd, this->pointCnt);


        /*****************************************
        ************ NEIGHBOUR SEARCH ************
        ******************************************/
        this->threads = threads < this->searchPointCnt ? threads : this->searchPointCnt;
        this->blocks = (this->searchPointCnt % threads != 0) ? (this->searchPointCnt / threads + 1)
                                                             : (this->searchPointCnt / threads);
       // printf("this->searchPointCnt: %d\n", this->searchPointCnt);
        
        this->h_neighbourCounts = (unsigned int*)malloc(this->searchPointCnt * sizeof(uint));
        cudaMalloc((void**)&this->d_neighbourCounts, (this->searchPointCnt) * sizeof(unsigned int));
        cudaMalloc((void**)&this->d_neighbourCountsIndexSum, (this->searchPointCnt) * sizeof(unsigned int));
        this->h_neighbourCountsIndexSum = (unsigned int*)malloc(this->searchPointCnt * sizeof(uint));
		
		// for testing
        //cudaMemset(this->d_neighbours, 0, (this->pointCnt * max_num_neigbors) * sizeof(uint));
        //cudaMemset(this->d_neighbourCounts, 0, (this->searchPointCnt) * sizeof(uint));
        
		/********************
        ***** FIRST RUN *****
        *********************/
        // CUDA: searching neighbours for each point on GPU (only counting neigbours)
        NeighbourSearchCUDA(this->blocks, this->threads, this->d_searchPointPos, this->d_pointPos, this->searchPointCnt,
            this->moveOriginPoint, this->gridDimensions, this->d_cellStartEnd, this->d_sorted_insertGrid,
            this->d_neighbours, this->d_neighbourCounts, this->d_neighbourCountsIndexSum, this->searchRadius,
            this->gridCellSize, false);
       

        /***** NEIGHBOUR ALLOC *****/
        //inclusive scan to get also the total number of neighbours
        prefixSumCUDA((int*)this->d_neighbourCountsIndexSum,
            (int*)this->d_neighbourCountsIndexSum + this->searchPointCnt,
            (int*)this->d_neighbourCountsIndexSum); // in-place scan

        cudaMemcpy(this->h_neighbourCountsIndexSum, this->d_neighbourCountsIndexSum,
            this->searchPointCnt * sizeof(unsigned int), cudaMemcpyDeviceToHost);

        this->totalNeighbourCnt = this->h_neighbourCountsIndexSum[this->searchPointCnt - 1] + 1;

		// allocate memory for neighbours
        cudaMalloc((void**)&this->d_neighbours, (this->totalNeighbourCnt * sizeof(unsigned int)));
        this->h_neighbours = (unsigned int*)calloc(this->totalNeighbourCnt, sizeof(uint));

		/*********************
        ***** SECOND RUN *****
        **********************/
        // Search and store the neigbours
        NeighbourSearchCUDA(this->blocks, this->threads, this->d_searchPointPos, this->d_pointPos, this->searchPointCnt,
            this->moveOriginPoint, this->gridDimensions, this->d_cellStartEnd, this->d_sorted_insertGrid,
            this->d_neighbours, this->d_neighbourCounts, this->d_neighbourCountsIndexSum, this->searchRadius,
            this->gridCellSize, true);
  
		cudaMemcpy(this->h_neighbours, this->d_neighbours, this->totalNeighbourCnt * sizeof(unsigned int),
                           cudaMemcpyDeviceToHost);
        cudaMemcpy(this->h_neighbourCounts, this->d_neighbourCounts, this->searchPointCnt * sizeof(unsigned int),
                           cudaMemcpyDeviceToHost);
      
		// create vector of vectors data structure of plain neighbours array "h_neighbours"
        this->neighbours.resize(this->searchPointCnt);
#pragma omp parallel for
        for (int i = 0; i < this->searchPointCnt; i++) {
            this->neighbours[i].resize(this->h_neighbourCounts[i]);
            for (int j = 0; j < this->h_neighbourCounts[i]; j++) {
                this->neighbours[i][j] = this->h_neighbours[this->h_neighbourCountsIndexSum[i] - j];
            }
		}

		// ONLY FOR CHECK RESULTS WITH NAIVE NEIGHBOR SEARCH IN COMPARISON
        bool verbose = false;
        if (false) {
            if (this->neighbours.size() == 7) {
                verbose = true;
			}
            vislib::Array<float> h_searchPoints;
            h_searchPoints.SetCount(this->searchPointCnt*4);

            cudaMemcpy(&h_searchPoints[0], this->d_searchPointPos, this->searchPointCnt * sizeof(float4),
                cudaMemcpyDeviceToHost);
			int mismatch = 0;
			float dist, detectionRadius;
			int cnt = 0;
            for (int i = 0; i < this->searchPointCnt; i++) {
				std::sort(this->neighbours[i].begin(), this->neighbours[i].end());
                if (verbose) {
                    printf("ID: %d\n", i);
                }
				cnt = 0;
                float4 a = {h_searchPoints[i * 4 + 0], h_searchPoints[i * 4 + 1], h_searchPoints[i * 4 + 2],
                    h_searchPoints[i * 4 + 3]};
                for (int j = 0; j < this->pointCnt; j++) {
					float4 b = {
						this->points[j].GetX(), this->points[j].GetY(), this->points[j].GetZ(), this->points[j].GetW()};
					dist = sqrtf(powf(b.x - a.x, 2) + powf(b.y - a.y, 2) + powf(b.z - a.z, 2));
					detectionRadius = a.w + b.w + 2.0f * this->searchRadius - 0.001f;
					if (dist < detectionRadius) {
                        if (verbose) {
                            printf("dist < detectionRadius: %f < %f\n", dist, detectionRadius);
                        }
                        if (this->neighbours[i].size() < cnt+1) {
                            if (verbose) {
                               printf("%d = %d \t i:%d j:%d\n", j, 0, i, cnt);
                           }
                            mismatch++;
                        } else {
                            if (verbose) {
                                printf("%d = %d \t i:%d j:%d\n", j, this->neighbours[i][cnt], i, cnt);
                            }
							if (j != this->neighbours[i][cnt]) {
								mismatch++;
                                printf("%d = %d \t i:%d j:%d\n", j, this->neighbours[i][cnt], i, cnt);
								if (verbose) {printf("-------------Mistmatches: %d\n", mismatch);}
							}
                        }
                        cnt++;
					}
				}
			   if (verbose) { printf("count: %d = %d\n", cnt, this->h_neighbourCounts[i]);}
			}
            
			printf("Mistmatches: %d\n", mismatch);
            if (mismatch != 0) {
                Sleep(5000);
            }
        }// END OF:  ONLY FOR CHECK RESULTS

        this->doCalcNeigbours = false;
        cleanUp();
    }
    return true ;
 
}

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
void FFRNNS::cleanUp() {
    // for GRID INSERT
    cudaFree(this->d_grid);
    cudaFree(this->d_insertGrid);
    cudaFree(this->d_pointPos);
    // for COUNTING SORT
    cudaFree(this->d_sorted_insertGrid);
    // for REINDEXING
    cudaFree(this->d_cellStartEnd);
    // for NEIGHBOUR SEARCH
    cudaFree(this->d_neighbours);
    cudaFree(this->d_neighbourCounts);
    cudaFree(this->d_neighbourCountsIndexSum);
    free(this->h_neighbourCounts);
    free(this->h_neighbours);
    free(this->h_neighbourCountsIndexSum);
}


/**
 * calculates the bounding box of the all point data (x,y,z)
 *
 * @this-param	this->points
 * @this-param	this->bbox
 * @src			prolint::FFRNNS::calcBoundingBox
 */
void FFRNNS::calcBoundingBox() {
    for (int i = 0; i < this->points.Count(); i++) {
        this->bbox.GrowToPoint(this->points[i].GetX(), this->points[i].GetY(), this->points[i].GetZ());
    }
    this->bbox.EnforcePositiveSize();
}


/**
 * calculates paramters of the grid(dimensions, grid cell size, min values of bbox)
 *
 * @this-param	this->points
 * @this-param	this->bbox
 * @src			prolint::FFRNNS::calcGridParams
 */
bool FFRNNS::calcGridParams() {

    // check if GridCellSize already calculated/set
    if (this->gridCellSize.x == 0 || this->gridParamsChanged) {
        // determine the max radius of all points
        float maxRadius = 0.f;
        for (int i = 0; i < this->points.Count(); i++) {
            maxRadius = maxRadius < this->points[i].GetW() ? this->points[i].GetW() : maxRadius;
        }
        int calc_gridCellSize = (int)ceilf(2.0f * maxRadius + 2.0f *this->searchRadius);

		// handle special case if points have no radius and the searchRadius is set to zero
        if( calc_gridCellSize <= 0){
            calc_gridCellSize  =  1;
            this->searchRadius = 0.001f;
		}
        setGridCellSize(calc_gridCellSize*1.10f);
    }

    // check if gridDimensions already calculated
    if ((this->gridDimensions.x == 0 && this->gridDimensions.y == 0 && this->gridDimensions.z == 0) ||
        this->gridParamsChanged) {

        float grid_x_max = this->bbox.GetRightTopBack().GetX();
        float grid_x_min = this->bbox.GetLeftBottomFront().GetX();

        float grid_y_max = this->bbox.GetRightTopBack().GetY();
        float grid_y_min = this->bbox.GetLeftBottomFront().GetY();

        float grid_z_min = this->bbox.GetRightTopBack().GetZ();
        float grid_z_max = this->bbox.GetLeftBottomFront().GetZ();

        float grid_x_range = grid_x_max - grid_x_min;
        float grid_y_range = grid_y_max - grid_y_min;
        float grid_z_range = grid_z_max - grid_z_min;

        // to be sure that grid dimensions always >= max of (x,y,z) of all data points
        // (enough gird cells in each direction (x,y,z))
        this->gridDimensions.x = (int)ceilf(grid_x_range / this->gridCellSize.x);
        this->gridDimensions.y = (int)ceilf(grid_y_range / this->gridCellSize.y);
        this->gridDimensions.z = (int)ceilf(grid_z_range / this->gridCellSize.z);
        this->gridSize = this->gridDimensions.x * this->gridDimensions.y * this->gridDimensions.z;

        this->moveOriginPoint.x = grid_x_min;
        this->moveOriginPoint.y = grid_y_min;
        this->moveOriginPoint.z = grid_z_min;
    }

    return true;
}

/**
 * CUDA error report (for standard CUDA functions)
 *
 * @param	code: the CUDA error code
 * @param	file: the CUDA file where the error occured
 * @param   line: the line of the CUDA file where the error occured
 * @src		prolint::FFRNNS::check_CUDA_err
 */
void FFRNNS::check_CUDA_err(cudaError_t code, const char* file, int line) {
    if (code != cudaSuccess) {
        std::string file_and_line = "(";
        file_and_line.append(file);
        file_and_line.append(";");
        file_and_line.append((const char*)(line));

        throw thrust::system_error(code, thrust::cuda_category(), file_and_line);
    }
}


/**
 * CUDA error report (for custom CUDA kernels/functions)
 *
 * @info	just returning the last CUDA ERROR that occured
 * @src		prolint::FFRNNS::check_CUDA_err
 */
void FFRNNS::check_CUDA_err(void) {
    cudaError err = cudaGetLastError();
    std::string errString = cudaGetErrorString(err);
    if (err != cudaSuccess) {
        printf("CUDA ERROR: %s\n", errString.c_str());
    }
}

/**
 * optional: sets grid cell size of neighbour search grid (cubic voxels)
 *			 if not set it is automatically determined [ceil(2 x maxRad + 2x searchRadius)]
 *
 * @param	gridCellSize (used for x,y,z extent of grid cells)
 * @src		prolint::FFRNNS::setGridCellSize
 */
void FFRNNS::setGridCellSize(int gridCellSize) { 
	// TODO: call directly setGridCellSize(x,y,z) when it is implemented 
	// instead of directly assign this->gridCellSize
	//setGridCellSize(gridCellSize, gridCellSize, gridCellSize);
    this->gridCellSize.x = gridCellSize;
    this->gridCellSize.y = gridCellSize;
    this->gridCellSize.z = gridCellSize;
    this->gridParamsChanged = true;
}

/**
 * optional: sets grid cell size of neighbour search grid (also non cubic voxels)
 *			 if not set it is automatically determined [ceil(2 x maxRad +2 x searchRadius)]
 *
 * @param	gridCellSize (used for x,y,z extent of grid cells)
 * @src		prolint::FFRNNS::setGridCellSize
 */
template <typename... Ts>
void FFRNNS::setGridCellSize(int gridCellSize_x, int gridCellSize_y, int gridCellSize_z, Ts&&...) {
    // TODO: needs to be really implemented
    static_assert(always_false<Ts...>::value, "setGridCellSize() for non cubic function needs to be implemented!");
    this->gridCellSize.x = gridCellSize_x;
    this->gridCellSize.y = gridCellSize_y;
    this->gridCellSize.z = gridCellSize_z;
    this->gridParamsChanged = true;
}




} // namespace prolint
} // namespace megamol
