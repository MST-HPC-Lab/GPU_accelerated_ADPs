#include<cuda_candidate_grouper.cuh>

/* Group a sub-set of MBRs and compute the MBR that covers this sub-set */
__global__ void candidate_grouper (ullong num_mbrs, double4 *rects, double4 *cells, ullong num_cells, ullong last_mbr_pos, 
            ullong num_mbrs_per_cell)
{
    for (ullong i = blockIdx.x * blockDim.x + threadIdx.x; i < num_cells; i += blockDim.x * gridDim.x) 
    {
        //rect1 xmin1; rect2 ymin1; rect3 xmax1; rect4 ymax1;
        double min_x = rects[num_mbrs_per_cell*i].x;
        double min_y = rects[num_mbrs_per_cell*i].y; 
        double max_x = rects[num_mbrs_per_cell*i].z;
        double max_y = rects[num_mbrs_per_cell*i].w;
                    
        int max_size = num_mbrs_per_cell;

        /* if last chunk */
        if (i == (num_cells-1))
            max_size = last_mbr_pos;
        
        /* input is sorted by minx, so no need to compute minx*/
        for(int k = 0; k < max_size; ++k)
        {
            //min_x = fmin(min_x, rects[num_mbrs_per_cell*i+k].x);
            min_y = fmin(min_y, rects[num_mbrs_per_cell*i+k].y);
            max_x = fmax(max_x, rects[num_mbrs_per_cell*i+k].z);
            max_y = fmax(max_y, rects[num_mbrs_per_cell*i+k].w);
        }

        cells[i].x = min_x;
        cells[i].y = min_y;
        cells[i].z = max_x;
        cells[i].w = max_y;
    }
}
