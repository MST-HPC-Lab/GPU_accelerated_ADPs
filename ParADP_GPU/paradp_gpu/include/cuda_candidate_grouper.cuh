#ifndef __CUDA_CANDIDATE_GROUPER_CUH_INCLUDED__
#define __CUDA_CANDIDATE_GROUPER_CUH_INCLUDED__
#include <global_var.h>

/* Internal libs */

/* Global functions */
__global__ void candidate_grouper(ullong num_mbrs, double4 *rects, double4 *cells, ullong num_cells, ullong last_mbr_pos, ullong num_mbrs_per_cell);

/* device functions */

#endif //ndef __CUDA_CANDIDATE_GROUPER_CUH_INCLUDED__
