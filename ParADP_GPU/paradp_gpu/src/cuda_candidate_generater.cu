#include<stdio.h>
#include<cuda_candidate_generater.cuh>

/* device functions */
__device__ uint log2_uint(uint value);

__device__ ullong log2_ull(ullong value);

/* MISC. */
/* Table for quick log uint_64*/
__constant__ uint log_tab_ull[64] = {
    63,  0, 58,  1, 59, 47, 53,  2,
    60, 39, 48, 27, 54, 33, 42,  3,
    61, 51, 37, 40, 49, 18, 28, 20,
    55, 30, 34, 11, 43, 14, 22,  4,
    62, 57, 46, 52, 38, 26, 32, 41,
    50, 36, 17, 19, 29, 10, 13, 21,
    56, 45, 25, 31, 35, 16,  9, 12,
    44, 24, 15,  8, 23,  7,  6,  5};


/* Table for quick log uint_32*/
__constant__ uint log_tab_uint[32] = {
    0,  9,  1, 10, 13, 21,  2, 29,
   11, 14, 16, 18, 22, 25,  3, 30,
    8, 12, 20, 28, 15, 17, 24,  7,
   19, 27, 23,  6, 26,  5,  4, 31};


/* generate MBR candidates */
__global__ void candidate_generater(ullong num_mbrs_1, double4 *rects_1, ullong *rects_1_weights, ullong last_mbr_pos_1, ullong num_mbrs_2, double4 *rects_2, 
                                    ullong *rects_2_weights, ullong last_mbr_pos_2, double4 *cells_1, ullong num_cells_1, double4 *cells_2, ullong num_cells_2, 
                                    ullong num_mbrs_per_cell, ullong *num_candidates, Candidates *dev_candidates, ullong *dev_global_weight)
{
    // if (blockIdx.x * blockDim.x + threadIdx.x == 0){
    //     int deviceId;
    //     cudaGetDevice(&deviceId);

    //     printf("FROM KERNEL candidate_grouper: Device ID = %d\n", deviceId);
    // }
    
    /* For each cell in layer 1, find intersections in every cell in layer 2*/
    ullong i = blockIdx.x;

    for (ullong j = threadIdx.x; j < num_cells_2; j += blockDim.x)  
    {	
        if( !((cells_1[i].x > cells_2[j].z) || (cells_2[j].x > cells_1[i].z) || 
            (cells_1[i].w < cells_2[j].y) || (cells_2[j].w < cells_1[i].y)) 
            )
        {	
            int max_pos_1 = num_mbrs_per_cell;
            if (i == (num_cells_1-1))
                max_pos_1 = last_mbr_pos_1;
            
            int max_pos_2 = num_mbrs_per_cell;
            if (j == (num_cells_2-1))
                max_pos_2 = last_mbr_pos_2;
            
            for (int m = num_mbrs_per_cell * i; m < num_mbrs_per_cell * i + max_pos_1; m++)  
            {
                for (int n = num_mbrs_per_cell * j; n < num_mbrs_per_cell * j + max_pos_2; n++)  
                {	
                    if( !((rects_1[m].x > rects_2[n].z) || (rects_2[n].x > rects_1[m].z) 
                        ||(rects_1[m].y > rects_2[n].w) || (rects_2[n].y > rects_1[m].w))
                        )
                    {
                        /* MBRs intersection */
                        double intMinX = rects_1[m].x > rects_2[n].x ? rects_1[m].x : rects_2[n].x;
                        double intMinY = rects_1[m].y > rects_2[n].y ? rects_1[m].y : rects_2[n].y;
                        double intMaxX = rects_1[m].z < rects_2[n].z ? rects_1[m].z : rects_2[n].z;
                        double intMaxY = rects_1[m].w < rects_2[n].w ? rects_1[m].w : rects_2[n].w;

                        double temp_x = (intMinX + intMaxX) / 2.0;
                        double temp_y = (intMinY + intMaxY) / 2.0;

                        ullong temp_weight = rects_1_weights[m] + rects_2_weights[n];
                        temp_weight = temp_weight * log2_ull(temp_weight);

                        /* Record candidates */
                        ullong atomic_pos = atomicAdd(&num_candidates[0], 1);
                        
                        if (atomic_pos < MAX_NUM_CANDIDATES)
                        {
                            dev_candidates[0].m_x[atomic_pos] = temp_x;
                            dev_candidates[0].m_y[atomic_pos] = temp_y;
                            dev_candidates[0].weights[atomic_pos] = temp_weight;
                        }
                            
                        atomicAdd(&dev_global_weight[0], temp_weight);
                    }
                }
            }
        }
    }
}

/* Fast log2 for 32-bit unsigned int */
__device__ uint log2_uint(uint value)
{
    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;
    return log_tab_uint[(uint)(value*0x07C4ACDD) >> 27];
}

/* Fast log2 for 64-bit unsigned int */
__device__ ullong log2_ull(ullong value)
{
    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;
    value |= value >> 32;
    return log_tab_ull[((ullong)((value - (value >> 1))*0x07EDD5E59A4E28C2)) >> 58];
}

