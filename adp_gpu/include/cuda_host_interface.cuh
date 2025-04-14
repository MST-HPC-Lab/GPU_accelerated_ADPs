#ifndef __CUDA_HOST_INTERFACE_CUH_INCLUDED__
#define __CUDA_HOST_INTERFACE_CUH_INCLUDED__
#include <global_var.h>

/* External libs */
#include <sys/time.h>

/* Internal libs */
/* CUDA libs should be included inside cuda_host_interface.cu file */

/* Wrapper Functions */
void GPU_Utils_wakeup_wrapper(ullong num_wakeup);

ullong CUDA_Utils_gpus_wrapper(MBR_Array *mbr_arr_1, MBR_Array *mbr_arr_2, MBR global_mbr, ullong target_num_cells, 
                                MBR_Array *partition_results, int gpu_id, int num_gpus);

ullong CUDA_Utils_candidate_generate_wrapper(MBR_Array *mbr_arr_1, MBR_Array *mbr_arr_2, ullong num_mbrs_per_cell, 
                                Candidates *candidates, ullong *global_weight);

ullong* CUDA_Utils_quadtree_partition_wrapper(Candidates *candidates, ullong num_candidates, MBR global_mbr, ullong global_weight, 
                                ullong target_partitions, MBR_Array *partitioned_cells);

ullong CUDA_Utils_adp_wrapper(MBR_Array *mbr_arr_1, MBR_Array *mbr_arr_2, MBR global_mbr, 
                                ullong target_num_cells, MBR_Array *partition_results);

int gpu_count();

#endif //ndef __CUDA_HOST_INTERFACE_CUH_INCLUDED__
