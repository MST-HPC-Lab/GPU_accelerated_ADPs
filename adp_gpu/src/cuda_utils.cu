#include<cuda_utils.cuh>

/* Wake up gpu before other work */
__global__ void wakeup_gpu(ulong num_wakeup)
{	
	for (ulong i = blockIdx.x * blockDim.x + threadIdx.x; i < num_wakeup; i += blockDim.x * gridDim.x) 
    {
		int dummy = WAKEUP_VAR;
		dummy = dummy * WAKEUP_VAR + WAKEUP_DUMMY;
    }
}

/* A naive timer*/
__host__ double get_current_time()
{
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);

    return ( (double) tv.tv_sec + (double) tv.tv_usec * 1.e-6 );
}
