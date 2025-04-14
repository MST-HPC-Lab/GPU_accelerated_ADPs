#ifndef __CUDA_UTILS_CUH_INCLUDED__
#define __CUDA_UTILS_CUH_INCLUDED__
#include <global_var.h>

/* External libs */
#include <sys/time.h>

/* Internal libs */

/* Sorter Structs*/
#ifndef CUDA_DOUBLE_4_X_SORT_STRUCT
#define CUDA_DOUBLE_4_X_SORT_STRUCT
struct sort_d4_x
{
    __host__ __device__
    bool operator()(const double4 &a, const double4 &b) const
    {
        return (a.x < b.x);
    }
};
#endif //CUDA_DOUBLE_4_X_SORT_STRUCT

/* Host functions*/
__host__ double get_current_time();
/* Global functions */
__global__ void wakeup_gpu(ulong num_wakeup);


#endif //ndef __CUDA_UTILS_CUH_INCLUDED__
