// GTX Model | Memory|	Cores	| SMX	|	Blocks	| CN Threads |	 CNH Threads Linux	| CNH Threads Win
//	1080	|	8	|	2560	| 20	|  80/65536 |	32/1024	 |			25			|     20
#include "SESTrilaterationCUDA.cuh"


// helper...
__host__ __device__ float distance2Vectors(float3 v1, float3 v2) {
    float dis = sqrtf(powf(v2.x - v1.x, 2) + powf(v2.y - v1.y, 2) + powf(v2.z - v1.z, 2));
    return dis;
}

__host__ __device__ float distance2Vectors4f(float4 v1, float4 v2) {
    float dis = sqrtf(powf(v2.x - v1.x, 2) + powf(v2.y - v1.y, 2) + powf(v2.z - v1.z, 2));
    return dis;
}

// INSERT ATOMS into bins and count atoms per bin (atomicAdd)
__global__ void grid_insert(int* d_grid, int2* d_insert_grid, float4* d_atomPos, float3 move_origin_pdb, int3 griddims,
    int atom_count, float grid_step_size) {

    unsigned int x, y, z, grid_cell = 0;
    int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    int threadId = blockId * blockDim.x + threadIdx.x;

    if (threadId >= atom_count) return;

    float grid_spacing = 1.0 / grid_step_size;
    float4 atom = d_atomPos[threadId];

    x = (unsigned int)floorf((atom.x - move_origin_pdb.x) * grid_spacing);
    y = (unsigned int)floorf((atom.y - move_origin_pdb.y) * grid_spacing);
    z = (unsigned int)floorf((atom.z - move_origin_pdb.z) * grid_spacing);

    // get cell_number of grid
    grid_cell = x * (griddims.y * griddims.z) + y * griddims.z + z;

    // add +1 for cell_counter
    int internIndex = atomicAdd(&d_grid[grid_cell], 1);


    d_insert_grid[threadId] = {(int)grid_cell, internIndex};
}

__global__ void grid_insert_intersection(int* d_grid, int2* d_insert_grid, float3* d_intersections,
    float3 move_origin_pdb, int3 griddims, int atom_count, float grid_step_size) {

    unsigned int x, y, z, grid_cell = 0;
    int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    int threadId = blockId * blockDim.x + threadIdx.x;

    if (threadId >= atom_count) return;

    float grid_spacing = 1.0 / grid_step_size;
    float3 atom = d_intersections[threadId];

    x = (unsigned int)floorf((atom.x - move_origin_pdb.x) * grid_spacing);
    y = (unsigned int)floorf((atom.y - move_origin_pdb.y) * grid_spacing);
    z = (unsigned int)floorf((atom.z - move_origin_pdb.z) * grid_spacing);

    // get cell_number of grid
    grid_cell = x * (griddims.y * griddims.z) + y * griddims.z + z;

    // add +1 for cell_counter
    int internIndex = atomicAdd(&d_grid[grid_cell], 1);


    d_insert_grid[threadId] = {(int)grid_cell, internIndex};
}


// COUNTING SORT on GPU
__global__ void counting_sort(int2* d_insert_grid, int2* d_sorted, int atom_count, int* d_pre_sum) {

    int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    int atom = blockId * blockDim.x + threadIdx.x;

    if (atom >= atom_count) return;

    d_sorted[d_pre_sum[d_insert_grid[atom].x] - d_insert_grid[atom].y - 1].x = atom;
    d_sorted[d_pre_sum[d_insert_grid[atom].x] - d_insert_grid[atom].y - 1].y = d_insert_grid[atom].x;
}


// REINDEXING on GPU (cell_start, cell_end)
__global__ void reindex(int2* d_sorted, int2* d_cell_start_end, int atom_count) {

    int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    int j = blockId * blockDim.x + threadIdx.x;

    if (j >= atom_count) {
        return;
    }

    if (j > 0) { // last cell index != current cell index  =>(start position of cell)
        if (d_sorted[j - 1].y != d_sorted[j].y) {
            d_cell_start_end[(d_sorted[j].y)].x = j;     // start position
            d_cell_start_end[(d_sorted[j - 1].y)].y = j; // end position before
        }
    } else {
        d_cell_start_end[(d_sorted[0].y)].x = 0;
        d_cell_start_end[d_sorted[atom_count - 1].y].y = atom_count;
    }
}


// searching neighbours for each atom on GPU (cell_start, cell_end)
__global__ void NeighbourSearch(float4* d_atomPos, int atom_count, float3 move_origin_pdb, int3 griddims,
    int2* d_cell_start_end, int2* d_sorted, unsigned int* d_neighbours, unsigned int* d_reduced_neighbours,
    unsigned int* d_neighbourCounts, unsigned int* d_reduced_neighbourCounts, float probeRad, float grid_step_size) {

    int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    int i = blockId * blockDim.x + threadIdx.x;

    if (i >= atom_count) return;
    float grid_spacing = 1.0 / grid_step_size;
    int atomX, atomY, atomZ, start, end, cell, neighbourCount = 0, reduced_neighbourCount = 0;

    atomX = ((d_atomPos[i].x - move_origin_pdb.x) * grid_spacing) - 1;
    atomY = ((d_atomPos[i].y - move_origin_pdb.y) * grid_spacing) - 1;
    atomZ = ((d_atomPos[i].z - move_origin_pdb.z) * grid_spacing) - 1;

    for (int x = atomX; x < atomX + 3; x++) {
        for (int y = atomY; y < atomY + 3; y++) {
            for (int z = atomZ; z < atomZ + 3; z++) {
                if (x >= 0 && y >= 0 && z >= 0 && x < griddims.x && y < griddims.y &&
                    z < griddims.z) { // boundaries of grid
                    cell = x * (griddims.y * griddims.z) + y * griddims.z + z;
                    start = d_cell_start_end[cell].x;
                    end = d_cell_start_end[cell].y;
                    for (int k = start; k < end; k++) {
                        int atom_pos_index = d_sorted[k].x;
                        if (i != atom_pos_index) {
                            // check for radius
                            float detectionRadius =
                                d_atomPos[i].w + d_atomPos[atom_pos_index].w + 2.0f * probeRad - 0.001f;
                            float atomDistance = distance2Vectors4f(d_atomPos[atom_pos_index], d_atomPos[i]);
                            if (atomDistance < detectionRadius) {
                                d_neighbours[i * max_num_neigbors + neighbourCount] = atom_pos_index;
                                neighbourCount++;
                                // store only neighbors that have a smaller index than the current atom
                                // --> this avoids duplicate combinations
                                if (i < atom_pos_index) {
                                    d_reduced_neighbours[i * max_num_neigbors + reduced_neighbourCount] =
                                        atom_pos_index;
                                    reduced_neighbourCount++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    d_reduced_neighbourCounts[i] = reduced_neighbourCount;
    d_neighbourCounts[i] = neighbourCount;
}

// searching neighbours for each atom on GPU (cell_start, cell_end)
__global__ void NeighbourSearch(float4* d_atomPos, int atom_count, float3 move_origin_pdb, int3 griddims,
    int2* d_cell_start_end, int2* d_sorted, unsigned int* d_neighbours, unsigned int* d_neighbourCounts, float probeRad,
    float grid_step_size) {

    int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    int i = blockId * blockDim.x + threadIdx.x;

    if (i >= atom_count) return;
    float grid_spacing = 1.0 / grid_step_size;
    int atomX, atomY, atomZ, start, end, cell, neighbourCount = 0, reduced_neighbourCount = 0;

    atomX = ((d_atomPos[i].x - move_origin_pdb.x) * grid_spacing) - 1;
    atomY = ((d_atomPos[i].y - move_origin_pdb.y) * grid_spacing) - 1;
    atomZ = ((d_atomPos[i].z - move_origin_pdb.z) * grid_spacing) - 1;

    for (int x = atomX; x < atomX + 3; x++) {
        for (int y = atomY; y < atomY + 3; y++) {
            for (int z = atomZ; z < atomZ + 3; z++) {
                if (x >= 0 && y >= 0 && z >= 0 && x < griddims.x && y < griddims.y &&
                    z < griddims.z) { // boundaries of grid
                    cell = x * (griddims.y * griddims.z) + y * griddims.z + z;
                    start = d_cell_start_end[cell].x;
                    end = d_cell_start_end[cell].y;
                    for (int k = start; k < end; k++) {
                        int atom_pos_index = d_sorted[k].x;
                        if (i != atom_pos_index) {
                            // check for radius
                            float detectionRadius =
                                d_atomPos[i].w + d_atomPos[atom_pos_index].w + 2.0f * probeRad - 0.001f;
                            float atomDistance = distance2Vectors4f(d_atomPos[atom_pos_index], d_atomPos[i]);
                            if (atomDistance < detectionRadius) {
                                d_neighbours[i * max_num_neigbors + neighbourCount] = atom_pos_index;
                                neighbourCount++;
                            }
                        }
                    }
                }
            }
        }
    }
    d_neighbourCounts[i] = neighbourCount;
}


// searching neighbours for each atom on GPU (cell_start, cell_end)
__global__ void NeighbourSearch_intersections(float4* d_intersections, int filtered_dim, float3 move_origin_pdb,
    int3 griddims, int2* d_cell_start_end, int2* d_sorted, float4* d_neighbours, int* d_neighbourCounts, float probeRad,
    float grid_step_size) {

    int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    int i = blockId * blockDim.x + threadIdx.x;

    if (i >= filtered_dim) return;
    float grid_spacing = 1.0 / grid_step_size;
    int atomX, atomY, atomZ, start, end, cell, neighbourCount = 0, reduced_neighbourCount = 0;

    atomX = ((d_intersections[i].x - move_origin_pdb.x) * grid_spacing) - 1;
    atomY = ((d_intersections[i].y - move_origin_pdb.y) * grid_spacing) - 1;
    atomZ = ((d_intersections[i].z - move_origin_pdb.z) * grid_spacing) - 1;

    for (int x = atomX; x < atomX + 3; x++) {
        for (int y = atomY; y < atomY + 3; y++) {
            for (int z = atomZ; z < atomZ + 3; z++) {
                if (x >= 0 && y >= 0 && z >= 0 && x < griddims.x && y < griddims.y &&
                    z < griddims.z) { // boundaries of grid
                    cell = x * (griddims.y * griddims.z) + y * griddims.z + z;
                    start = d_cell_start_end[cell].x;
                    end = d_cell_start_end[cell].y;
                    for (int k = start; k < end; k++) {
                        int atom_pos_index = d_sorted[k].x;
                        if (i != atom_pos_index) {
                            // check for radius
                            float detectionRadius = 2.0f * probeRad - 0.001f;
                            float atomDistance =
                                distance2Vectors4f(d_intersections[atom_pos_index], d_intersections[i]);
                            if (atomDistance < detectionRadius) {
                                // printf("atomDistance < detectionRadius:  %f < %f \n", atomDistance, detectionRadius);
                                d_neighbours[i * 50 + neighbourCount].x = d_intersections[atom_pos_index].x;
                                d_neighbours[i * 50 + neighbourCount].y = d_intersections[atom_pos_index].y;
                                d_neighbours[i * 50 + neighbourCount].z = d_intersections[atom_pos_index].z;
                                d_neighbours[i * 50 + neighbourCount].w = probeRad;
                                neighbourCount++;
                            }
                        }
                    }
                }
            }
        }
    }
    d_neighbourCounts[i] = neighbourCount;
}

__global__ void prepare_Number_of_Combinations(int atom_count, unsigned int* d_reduced_neighbours,
    unsigned int* d_reduced_neighbourCounts, int* d_neighbour_combination_numbers) {

    int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    int i = blockId * blockDim.x + threadIdx.x;

    if (i >= atom_count) return;
    int n = d_reduced_neighbourCounts[i];
    d_neighbour_combination_numbers[i] = n * (n - 1) / 2;
}


// wikipedia.org/wiki/Trilateration
__device__ void check_Solution_exists1(int* atomicIdx, int thread, float4 m0, float4 m1, float4 m2,
    int3* d_SpheresComb3, const float4* d_atomPos, int d_calcInter_On_Off, unsigned int* d_neighbours,
    unsigned int* d_neighbourCounts, float4* d_intersections, int4* d_filtered_SpheresComb3, float probeRad,
    int memSize, int* d_atomPosIdx) {
    float3 ex, ey, ez;
    float i, j, d;
    float x, y, z;
    float3 s1, s2;

    /*************************/
    /******* PRE WORK ********/
    /*************************/
    d = distance2Vectors({m0.x, m0.y, m0.z}, {m1.x, m1.y, m1.z});
    ex = {(m1.x - m0.x) / d, (m1.y - m0.y) / d, (m1.z - m0.z) / d};
    i = ex.x * (m2.x - m0.x) + ex.y * (m2.y - m0.y) + ex.z * (m2.z - m0.z);

    // ey
    float3 m2_minus_m0 = {
        m2.x - m0.x,
        m2.y - m0.y,
        m2.z - m0.z,
    };
    float3 m2_minus_m0_minus_i_mult_ex = {m2_minus_m0.x - ex.x * i, m2_minus_m0.y - ex.y * i, m2_minus_m0.z - ex.z * i};
    float d_ey = distance2Vectors(m2_minus_m0_minus_i_mult_ex, {0, 0, 0});
    ey = {m2_minus_m0_minus_i_mult_ex.x / d_ey, m2_minus_m0_minus_i_mult_ex.y / d_ey,
        m2_minus_m0_minus_i_mult_ex.z / d_ey};

    j = ey.x * m2_minus_m0.x + ey.y * m2_minus_m0.y + ey.z * m2_minus_m0.z;

    /***********************************/
    /******* CACL INTERSECTIONs ********/
    /***********************************/

    x = (powf(m0.w, 2) - powf(m1.w, 2) + powf(d, 2)) / (2 * d);

    y = ((powf(m0.w, 2) - powf(m2.w, 2) + powf(i, 2) + powf(j, 2)) / (2 * j)) - ((i / j) * x);

    float z_input = powf(m0.w, 2) - powf(x, 2) - powf(y, 2);

    int s1Valid = 1;
    int s2Valid = 1;
    if (z_input > 0) {
        z = sqrt(z_input);


        ez = {
            ex.y * ey.z - ex.z * ey.y,
            ex.z * ey.x - ex.x * ey.z,
            ex.x * ey.y - ex.y * ey.x,
        };

        s1 = {m0.x + x * ex.x + y * ey.x + z * ez.x, m0.y + x * ex.y + y * ey.y + z * ez.y,
            m0.z + x * ex.z + y * ey.z + z * ez.z};

        s2 = {m0.x + x * ex.x + y * ey.x - z * ez.x, m0.y + x * ex.y + y * ey.y - z * ez.y,
            m0.z + x * ex.z + y * ey.z - z * ez.z};
    } else {
        // no real solution found
        return;
    }

    // check neighbors for intersection
    unsigned int nIdx;
    float dist, nRad;
    float3 nPos;

    for (unsigned int cnt = 0; cnt < d_neighbourCounts[d_SpheresComb3[thread].x]; cnt++) {
        nIdx = d_neighbours[d_SpheresComb3[thread].x * max_num_neigbors + cnt];
        nPos = {d_atomPos[nIdx].x, d_atomPos[nIdx].y, d_atomPos[nIdx].z};
        nRad = d_atomPos[nIdx].w;
        dist = distance2Vectors(s1, nPos);
        if (dist < (probeRad + nRad) && nIdx != (d_SpheresComb3[thread].y) && nIdx != (d_SpheresComb3[thread].z)) {
            s1Valid = -1;
        }
        dist = distance2Vectors(s2, nPos);
        if (dist < (probeRad + nRad) && nIdx != (d_SpheresComb3[thread].y) && nIdx != (d_SpheresComb3[thread].z)) {
            s2Valid = -1;
        }
    }

    int writeIdx;


    if (s1Valid > 0) {
        writeIdx = atomicAdd(atomicIdx, 1);
        if (writeIdx < memSize) {
            d_filtered_SpheresComb3[writeIdx] = {
                d_SpheresComb3[thread].x, d_SpheresComb3[thread].y, d_SpheresComb3[thread].z, 1};
            d_intersections[writeIdx] = {s1.x, s1.y, s1.z, probeRad};
#ifdef reduceAtoms
            d_atomPosIdx[d_SpheresComb3[thread].x] = d_SpheresComb3[thread].x;
            d_atomPosIdx[d_SpheresComb3[thread].y] = d_SpheresComb3[thread].y;
            d_atomPosIdx[d_SpheresComb3[thread].z] = d_SpheresComb3[thread].z;
#endif
        }
    }
    if (s2Valid > 0) {
        writeIdx = atomicAdd(atomicIdx, 1);
        if (writeIdx < memSize) {
            d_filtered_SpheresComb3[writeIdx] = {
                d_SpheresComb3[thread].x, d_SpheresComb3[thread].y, d_SpheresComb3[thread].z, 1};
            d_intersections[writeIdx] = {s2.x, s2.y, s2.z, probeRad};
#ifdef reduceAtoms
            d_atomPosIdx[d_SpheresComb3[thread].x] = d_SpheresComb3[thread].x;
            d_atomPosIdx[d_SpheresComb3[thread].y] = d_SpheresComb3[thread].y;
            d_atomPosIdx[d_SpheresComb3[thread].z] = d_SpheresComb3[thread].z;
#endif
        }
    }
}


__global__ void check_3Spheres_Overlap(int* atomicIdx, int3* d_SpheresComb3, const float4* d_atomPos,
    int d_calcInter_On_off, unsigned int* d_neighbours, unsigned int* d_neighbourCounts, int number_of_3atomCombos,
    float probeRad, int memSize, float4* d_intersections, int4* d_filtered_SpheresComb3, int* d_atomPosIdx) {


    int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    int i = blockId * blockDim.x + threadIdx.x;

    if (i >= number_of_3atomCombos) return;

    float3 v0 = {d_atomPos[d_SpheresComb3[i].x].x, d_atomPos[d_SpheresComb3[i].x].y, d_atomPos[d_SpheresComb3[i].x].z};
    float3 v1 = {d_atomPos[d_SpheresComb3[i].y].x, d_atomPos[d_SpheresComb3[i].y].y, d_atomPos[d_SpheresComb3[i].y].z};
    float3 v2 = {d_atomPos[d_SpheresComb3[i].z].x, d_atomPos[d_SpheresComb3[i].z].y, d_atomPos[d_SpheresComb3[i].z].z};

    float dis12 = distance2Vectors(v1, v2);

    // check whether sphere s2 and s3 are neighbors
    if (dis12 < d_atomPos[d_SpheresComb3[i].y].w + d_atomPos[d_SpheresComb3[i].z].w + 2 * probeRad) {

        float4 m0 = {v0.x, v0.y, v0.z, d_atomPos[d_SpheresComb3[i].x].w + probeRad};
        float4 m1 = {v1.x, v1.y, v1.z, d_atomPos[d_SpheresComb3[i].y].w + probeRad};
        float4 m2 = {v2.x, v2.y, v2.z, d_atomPos[d_SpheresComb3[i].z].w + probeRad};

        check_Solution_exists1(atomicIdx, i, m0, m1, m2, d_SpheresComb3, d_atomPos, d_calcInter_On_off, d_neighbours,
            d_neighbourCounts, d_intersections, d_filtered_SpheresComb3, probeRad, memSize, d_atomPosIdx);
    }
}


__global__ void calc_3spheresCombos(int atom_count, unsigned int* d_reduced_neighbours,
    unsigned int* d_reduced_neighbourCounts, int* d_neighbour_combination_numbers, int3* d_SpheresComb3) {

    int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    int i = blockId * blockDim.x + threadIdx.x;

    if (i >= atom_count) return;
    int cnt = 1;

    for (int k = 0; k < d_reduced_neighbourCounts[i]; k++) {
        for (int j = k + 1; j < d_reduced_neighbourCounts[i]; j++) {
            d_SpheresComb3[d_neighbour_combination_numbers[i] - cnt] = {
                i, static_cast<int>(d_reduced_neighbours[i * max_num_neigbors + k]), static_cast<int>(d_reduced_neighbours[i * max_num_neigbors + j])};
            cnt++;
        }
    }
}

__global__ void writeTorusAxes(unsigned int numSphericalTriangles, int4* sphericalTriangles, int2* torusAxes) {

    unsigned int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    unsigned int k = blockId * blockDim.x + threadIdx.x;
    // don't overshoot the possible solutions
    if (k >= numSphericalTriangles) return;
    // loop over all intersections and write them to "torusAxis" array (sorted indices);
    torusAxes[3 * k + 0].x =
        sphericalTriangles[k].x < sphericalTriangles[k].y ? sphericalTriangles[k].x : sphericalTriangles[k].y;
    torusAxes[3 * k + 0].y =
        sphericalTriangles[k].x < sphericalTriangles[k].y ? sphericalTriangles[k].y : sphericalTriangles[k].x;
    torusAxes[3 * k + 1].x =
        sphericalTriangles[k].x < sphericalTriangles[k].z ? sphericalTriangles[k].x : sphericalTriangles[k].z;
    torusAxes[3 * k + 1].y =
        sphericalTriangles[k].x < sphericalTriangles[k].z ? sphericalTriangles[k].z : sphericalTriangles[k].x;
    torusAxes[3 * k + 2].x =
        sphericalTriangles[k].y < sphericalTriangles[k].z ? sphericalTriangles[k].y : sphericalTriangles[k].z;
    torusAxes[3 * k + 2].y =
        sphericalTriangles[k].y < sphericalTriangles[k].z ? sphericalTriangles[k].z : sphericalTriangles[k].y;
}

__global__ void cutSphereTriangles(unsigned int numAtoms, float4* atomPos, unsigned int sphereVertCnt,
    float3* sphereVert, unsigned int sphereTriaCnt, uint3* sphereTria, unsigned int* neighbours,
    unsigned int* neighbourCounts, uint3* outerSphereTria, uint3* cutSphereTria) {

    int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    int atomIdx = blockId * blockDim.x + threadIdx.x;

    if (atomIdx >= numAtoms) return;

    float3 pos = make_float3(atomPos[atomIdx].x, atomPos[atomIdx].y, atomPos[atomIdx].z);
    float rad = atomPos[atomIdx].w;

	// loop over all triangles or the current atom sphere and move them to the absolute position (i.e., add current atom position)
	for (unsigned int triaIdx = 0; triaIdx < sphereTriaCnt; triaIdx++) {
        uint3 currentTriaIdx = sphereTria[triaIdx];
        float3 vert0 = sphereVert[currentTriaIdx.x];
        float3 vert1 = sphereVert[currentTriaIdx.y];
        float3 vert2 = sphereVert[currentTriaIdx.z];
		// scale unit sphere to atom radius
        vert0.x *= rad;
        vert0.y *= rad;
        vert0.z *= rad;
        vert1.x *= rad;
        vert1.y *= rad;
        vert1.z *= rad;
        vert2.x *= rad;
        vert2.y *= rad;
        vert2.z *= rad;
		// move unit sphere vertex to atom position
        vert0.x += pos.x;
        vert0.y += pos.y;
        vert0.z += pos.z;
        vert1.x += pos.x;
        vert1.y += pos.y;
        vert1.z += pos.z;
        vert2.x += pos.x;
        vert2.y += pos.y;
        vert2.z += pos.z;
		// write triangle
        unsigned int atomTriaIdx = atomIdx * sphereTriaCnt + triaIdx;
        outerSphereTria[atomTriaIdx] = currentTriaIdx;
    }

}

extern "C" void gridInsertCuda(int blocks, int threads, int* d_grid, int2* d_insert_grid, float4* d_atomPos,
    float3 move_origin_pdb, int3 griddims, int atom_count, float grid_step_size) {

    // insert atoms into bins and count atoms per bin (atomicAdd)
    grid_insert<<<blocks, threads>>>(
        d_grid, d_insert_grid, d_atomPos, move_origin_pdb, griddims, atom_count, grid_step_size);
}

extern "C" void gridInsertCuda_intersection(int blocks, int threads, int* d_grid, int2* d_insert_grid,
    float3* d_intersections, float3 move_origin_pdb, int3 griddims, int atom_count, float grid_step_size) {

    grid_insert_intersection<<<blocks, threads>>>(
        d_grid, d_insert_grid, d_intersections, move_origin_pdb, griddims, atom_count, grid_step_size);
}

extern "C" void countingSortCuda(
    int blocks, int threads, int2* d_insert_grid, int2* d_sorted, int atom_count, int* d_pre_sum) {

    counting_sort<<<blocks, threads>>>(d_insert_grid, d_sorted, atom_count, d_pre_sum);
}

extern "C" void reindexCuda(int blocks, int threads, int2* d_sorted, int2* d_cell_start_end, int atom_count) {

    reindex<<<blocks, threads>>>(d_sorted, d_cell_start_end, atom_count);
}

extern "C" void NeighbourSearchCuda(int blocks, int threads, float4* d_atomPos, int atom_count, float3 move_origin_pdb,
    int3 griddims, int2* d_cell_start_end, int2* d_sorted, unsigned int* d_neighbours,
    unsigned int* d_reduced_neighbours, unsigned int* d_neighbourCounts, unsigned int* d_reduced_neighbourCounts,
    float probeRad, float grid_step_size) {

    if (d_reduced_neighbours && d_reduced_neighbourCounts) {
        NeighbourSearch<<<blocks, threads>>>(d_atomPos, atom_count, move_origin_pdb, griddims, d_cell_start_end,
            d_sorted, d_neighbours, d_reduced_neighbours, d_neighbourCounts, d_reduced_neighbourCounts, probeRad,
            grid_step_size);
    } else {
        NeighbourSearch<<<blocks, threads>>>(d_atomPos, atom_count, move_origin_pdb, griddims, d_cell_start_end,
            d_sorted, d_neighbours, d_neighbourCounts, probeRad, grid_step_size);
    }
}

extern "C" void NeighbourSearch4f_intersectionsCuda(int blocks, int threads, float4* d_atomPos, int filtered_dim,
    float3 move_origin_pdb, int3 griddims, int2* d_cell_start_end, int2* d_sorted, float4* d_neighbours,
    int* d_neighbourCounts, float probeRad, float grid_step_size) {

    NeighbourSearch_intersections<<<blocks, threads>>>(d_atomPos, filtered_dim, move_origin_pdb, griddims,
        d_cell_start_end, d_sorted, d_neighbours, d_neighbourCounts, probeRad, grid_step_size);
}

extern "C" void prepare_Number_of_CombinationsCuda(int blocks, int threads, int atom_count,
    unsigned int* d_reduced_neighbours, unsigned int* d_reduced_neighbourCounts, int* d_neighbour_combination_numbers) {

    prepare_Number_of_Combinations<<<blocks, threads>>>(
        atom_count, d_reduced_neighbours, d_reduced_neighbourCounts, d_neighbour_combination_numbers);
}


extern "C" void calc_3spheresCombosCuda(int blocks, int threads, int atom_count, unsigned int* d_reduced_neighbours,
    unsigned int* d_reduced_neighbourCounts, int* d_neighbour_combination_numbers, int3* d_SpheresComb3) {

    calc_3spheresCombos<<<blocks, threads>>>(
        atom_count, d_reduced_neighbours, d_reduced_neighbourCounts, d_neighbour_combination_numbers, d_SpheresComb3);
}

extern "C" void ScanCuda(int* start, int* end, int* dst) {
    thrust::inclusive_scan(
        thrust::device_ptr<int>(start), thrust::device_ptr<int>(end), thrust::device_ptr<int>(dst)); // in-place scan
}

extern "C" void check_3Spheres_OverlapCuda(int blocks, int threads, int* atomicIdx, int3* d_SpheresComb3,
    const float4* d_atomPos, int d_calcInter_On_off, unsigned int* d_neighbours, unsigned int* d_neighbourCounts,
    int number_of_3atomCombos, float probeRad, int memSize, float4* d_intersections, int4* d_filtered_SpheresComb3,
    int* d_atomPosIdx) {

    check_3Spheres_Overlap<<<blocks, threads>>>(atomicIdx, d_SpheresComb3, d_atomPos, d_calcInter_On_off, d_neighbours,
        d_neighbourCounts, number_of_3atomCombos, probeRad, memSize, d_intersections, d_filtered_SpheresComb3,
        d_atomPosIdx);
}

extern "C" void writeTorusAxesCUDA(
    int blocks, int threads, unsigned int numSphericalTriangles, int4* sphericalTriangles, int2* d_torusAxes) {
    writeTorusAxes<<<blocks, threads>>>(numSphericalTriangles, sphericalTriangles, d_torusAxes);
}


extern "C" void cutSphereTrianglesCUDA(int blocks, int threads, unsigned int numAtoms, float4* atomPos,
    unsigned int sphereVertCnt, float3* sphereVert, unsigned int sphereTriaCnt, uint3* sphereTria,
    unsigned int* neighbours, unsigned int* neighbourCounts,
    uint3* outerSphereTria, uint3* cutSphereTria) {
    cutSphereTriangles<<<blocks, threads>>>(numAtoms, atomPos, sphereVertCnt, sphereVert, sphereTriaCnt, sphereTria, neighbours, neighbourCounts, outerSphereTria, cutSphereTria);
}

extern "C" cudaError sortAndUniqueTorusAxesCUDA(unsigned int& numTorusAxes, int2* d_torusAxes) {
    thrust::sort(
        thrust::device_ptr<int2>(d_torusAxes), thrust::device_ptr<int2>(d_torusAxes + numTorusAxes), lessInt2X());
    const int numberOfUniqueValues = thrust::unique(thrust::device_ptr<int2>(d_torusAxes),
                                         thrust::device_ptr<int2>(d_torusAxes + numTorusAxes), equalInt2()) -
                                     thrust::device_ptr<int2>(d_torusAxes);
    numTorusAxes = static_cast<unsigned int>(numberOfUniqueValues);

    return cudaGetLastError();
}

extern "C" cudaError SortAndUniqueAtomPosIdxCUDA(unsigned int& size_AtomPosIdx, int* d_torusAxes) {

    thrust::sort(
        thrust::device_ptr<int>(d_torusAxes), thrust::device_ptr<int>(d_torusAxes + size_AtomPosIdx), lessIntX());
    const int numberOfUniqueValues = thrust::unique(thrust::device_ptr<int>(d_torusAxes),
                                         thrust::device_ptr<int>(d_torusAxes + size_AtomPosIdx), equalInt()) -
                                     thrust::device_ptr<int>(d_torusAxes);
    size_AtomPosIdx = static_cast<unsigned int>(numberOfUniqueValues);

    return cudaGetLastError();
}
