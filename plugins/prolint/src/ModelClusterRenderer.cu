#include "ModelClusterRenderer.cuh"

    /*
 * greater than comparison float
 */
__device__ bool greater_float::operator()(const float& lhs, const float& rhs) const {
   
    return (lhs > rhs);
}

__global__ void getTriangleToCamDistances(
    uint triaCnt, float* vertices, uint* indices, float3* camPos, float* d_TriangleToCamDistances) {
    unsigned int blockId = blockIdx.x + blockIdx.y * gridDim.x;
    unsigned int k = blockId * blockDim.x + threadIdx.x;

	if (k >= triaCnt) return;
    float3 tri_0;
    tri_0.x = vertices[indices[k * 3 + 0] * 3 + 0];
    tri_0.y = vertices[indices[k * 3 + 0] * 3 + 1];
    tri_0.z = vertices[indices[k * 3 + 0] * 3 + 2];

    float3 tri_1;
    tri_1.x = vertices[indices[k * 3 + 1] * 3 + 0];
    tri_1.y = vertices[indices[k * 3 + 1] * 3 + 1];
    tri_1.z = vertices[indices[k * 3 + 1] * 3 + 2];

    float3 tri_2;
    tri_2.x = vertices[indices[k * 3 + 2] * 3 + 0];
    tri_2.y = vertices[indices[k * 3 + 2] * 3 + 1];
    tri_2.z = vertices[indices[k * 3 + 2] * 3 + 2];

    // use midpoint
	float3 trianglePos;
 
    trianglePos.x = (tri_0.x + tri_1.x + tri_2.x) / 3.0f;
    trianglePos.y = (tri_0.y + tri_1.y + tri_2.y) / 3.0f;
    trianglePos.z = (tri_0.z + tri_1.z + tri_2.z) / 3.0f;
    float dis = distanceTwoVectors(trianglePos, camPos[0]);

    d_TriangleToCamDistances[k] = dis;
    //d_TriangleToCamDistances[k] = 2.0f;
}

extern "C" cudaError SortTrianglesDevice(
    uint triaCnt, uint3* sortedIndices, float* d_TriangleToCamDistances) {

    thrust::sort_by_key(thrust::device_ptr<float>(d_TriangleToCamDistances),
        thrust::device_ptr<float>(d_TriangleToCamDistances + triaCnt), thrust::device_ptr<uint3>(sortedIndices),
        greater_float());
    
    cudaDeviceSynchronize();
    return cudaGetLastError();
}


extern "C" cudaError getTriangleToCamDistancesCUDA( int blocks, int threads, uint triaCnt, float* vertices, uint* indices,  float3* camPos, float* d_TriangleToCamDistances) {

	getTriangleToCamDistances<<<blocks, threads>>>(triaCnt, vertices, indices, camPos, d_TriangleToCamDistances);
    return cudaGetLastError();

}

 // helper...
__host__ __device__ float distanceTwoVectors(float3 v1, float3 v2) {
    float dis = sqrtf(powf(v2.x - v1.x, 2) + powf(v2.y - v1.y, 2) + powf(v2.z - v1.z, 2));
    return dis;
}
