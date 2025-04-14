#ifndef __CUDA_CANDIDATE_GENERATOR_CUH_INCLUDED__
#define __CUDA_CANDIDATE_GENERATOR_CUH_INCLUDED__
#include <global_var.h>

/* Internal libs */

/* kernal functions */
__global__ void candidate_generater(ullong num_mbrs_1, double4 *rects_1, ullong *rects_1_weights, ullong last_mbr_pos_1, 
    ullong num_mbrs_2, double4 *rects_2, ullong *rects_2_weights, ullong last_mbr_pos_2, double4 *cells_1, 
    ullong num_cells_1, double4 *cells_2, ullong num_cells_2, ullong num_mbrs_per_cell, ullong *num_candidates, 
    Candidates *dev_candidates, ullong *dev_global_weight);

#endif //ndef __CUDA_CANDIDATE_GENERATOR_CUH_INCLUDED__
