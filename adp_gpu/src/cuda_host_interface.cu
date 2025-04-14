#include<cuda_host_interface.cuh>
#include<cuda_utils.cuh>
#include<cuda_candidate_grouper.cuh>
#include<cuda_candidate_generater.cuh>
#include<cuda_quadtree.cuh>

#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/device_malloc.h>
#include <thrust/device_ptr.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/gather.h>

#include <vector>


/* Global wake up methods */
void GPU_Utils_wakeup_wrapper(ullong num_wakeup)
{
	int num_blocks = (num_wakeup + DEFAULT_THREADS_PER_BLOCK - 1) / DEFAULT_THREADS_PER_BLOCK;				
	
	wakeup_gpu <<<num_blocks, DEFAULT_THREADS_PER_BLOCK>>> (num_wakeup);
	
	cudaDeviceSynchronize();
	return;
}

/* adp in gpu approach */
int gpu_count()
{
    int num_devices;
    cudaGetDeviceCount(&num_devices);

    return num_devices;
}

ullong CUDA_Utils_gpus_wrapper(MBR_Array *mbr_arr_1, MBR_Array *mbr_arr_2, MBR global_mbr, ullong target_num_cells, MBR_Array *partition_results, int gpu_id, int num_gpus)
{
    //size_t free, total;
    //printf("\n");
    //cudamemgetinfo(&free,&total);   
    //printf("%d kb free of total %d kb\n",free/1024,total/1024);

    int numGPUs = num_gpus;
    cudaSetDevice(gpu_id % numGPUs);//gpu_id);

    MBR_Array local_mbr_arr_1, local_mbr_arr_2;

    local_mbr_arr_1.size = mbr_arr_1->size / numGPUs;
    local_mbr_arr_2.size = mbr_arr_2->size;

    local_mbr_arr_1.mbrs = (MBR*)malloc(sizeof(MBR) * local_mbr_arr_1.size);
    memcpy(local_mbr_arr_1.mbrs, mbr_arr_1->mbrs + gpu_id * (local_mbr_arr_1.size), sizeof(MBR) * local_mbr_arr_1.size);
    local_mbr_arr_1.weights = (ullong*)malloc(sizeof(ullong) * local_mbr_arr_1.size);
    memcpy(local_mbr_arr_1.weights, mbr_arr_1->weights + gpu_id * (local_mbr_arr_1.size), sizeof(ullong) * local_mbr_arr_1.size);
    
    local_mbr_arr_2.mbrs = (MBR*)malloc(sizeof(MBR) * local_mbr_arr_2.size);
    memcpy(local_mbr_arr_2.mbrs, mbr_arr_2->mbrs, sizeof(MBR) * local_mbr_arr_2.size);
    local_mbr_arr_2.weights = (ullong*)malloc(sizeof(ullong) * local_mbr_arr_2.size);
    memcpy(local_mbr_arr_2.weights, mbr_arr_2->weights, sizeof(ullong) * local_mbr_arr_2.size);

    cudaDeviceSynchronize();
    
    //printf("%p %p %p %p \n", local_mbr_arr_1.mbrs, local_mbr_arr_2.mbrs, local_mbr_arr_1.weights, local_mbr_arr_2.weights);
    auto t_begin = get_current_time();
    double t_total, t_can, t_quad;

    Candidates candidates;
	ullong global_weight = 0ull;
    ullong num_candidates = 0ull;

	num_candidates = CUDA_Utils_candidate_generate_wrapper(&local_mbr_arr_1, &local_mbr_arr_2, DEFAULT_MBRS_PER_BLOCK, &candidates, &global_weight);

    auto t_can_end = get_current_time();

    CUDA_Utils_quadtree_partition_wrapper(&candidates, num_candidates, global_mbr, global_weight, 
            target_num_cells, partition_results);
    
    auto t_quad_end = get_current_time();
    t_total = t_quad_end - t_begin;
    t_can = t_can_end - t_begin;
    t_quad = t_quad_end - t_can_end;

    printf("Muiltiple GPUs total time %f, candidate generate %f, quadtree %f.\n", 
        t_total, t_can, t_quad);

    return num_candidates;
}


/* adp in gpu approach */
ullong CUDA_Utils_adp_wrapper(MBR_Array *mbr_arr_1, MBR_Array *mbr_arr_2, MBR global_mbr, ullong target_num_cells, MBR_Array *partition_results)
{
    //size_t free, total;
    //printf("\n");
    //cudamemgetinfo(&free,&total);   
    //printf("%d kb free of total %d kb\n",free/1024,total/1024);


    //cudaSetDevice(1);
    auto t_begin = get_current_time();
    double t_total, t_can, t_quad;

    Candidates candidates;
	ullong global_weight = 0ull;
    ullong num_candidates = 0ull;

	num_candidates = CUDA_Utils_candidate_generate_wrapper(mbr_arr_1, mbr_arr_2, DEFAULT_MBRS_PER_BLOCK, &candidates, &global_weight);

    auto t_can_end = get_current_time();

    CUDA_Utils_quadtree_partition_wrapper(&candidates, num_candidates, global_mbr, global_weight, 
            target_num_cells, partition_results);
    
    auto t_quad_end = get_current_time();
    t_total = t_quad_end - t_begin;
    t_can = t_can_end - t_begin;
    t_quad = t_quad_end - t_can_end;

    printf("cuda_utils_adp_wrapper total time %f, candidate generate %f, quadtree %f.\n", 
        t_total, t_can, t_quad);

    return num_candidates;
}

/* wrapper to generate candidates */
ullong CUDA_Utils_candidate_generate_wrapper(MBR_Array *mbr_arr_1, MBR_Array *mbr_arr_2, ullong num_mbrs_per_cell, 
                                            Candidates *candidates, ullong *global_weight)
{
#ifdef DEBUG
    //printf("CUDA_Utils_candidate_generate_wrapper: nums rects %lu : %lu.\n", mbr_arr_1->size, mbr_arr_2->size);
    /*
    //used to verify device info
    int warp_size = 0;
    int cuda_device = findCudaDevice(0, nullptr);
    cudaDeviceProp deviceProps;
    checkCudaErrors(cudaGetDeviceProperties(&deviceProps, cuda_device));
    warp_size = deviceProps.warpSize;
    */
#endif //ifdef DEBUG

    auto t_begin = get_current_time();
    double t_total, t_sort;

    /* pre-processing layer 1 */
    thrust::device_vector<double4> rects_1((double4*)(mbr_arr_1->mbrs), (double4*)(mbr_arr_1->mbrs) + mbr_arr_1->size);
    thrust::device_vector<ullong> rects_1_weights(mbr_arr_1->weights, mbr_arr_1->weights + mbr_arr_1->size);

    thrust::stable_sort_by_key(rects_1.begin(), rects_1.end(), rects_1_weights.begin(), sort_d4_x());

    double4* local_rects_1 = thrust::raw_pointer_cast(rects_1.data());
    ullong* local_rects_1_weights = thrust::raw_pointer_cast(rects_1_weights.data());

    /* pre-processing layer 2 */
    thrust::device_vector<double4> rects_2((double4*)(mbr_arr_2->mbrs), (double4*)(mbr_arr_2->mbrs) + mbr_arr_2->size);
    thrust::device_vector<ullong> rects_2_weights(mbr_arr_2->weights, mbr_arr_2->weights + mbr_arr_2->size);

    thrust::stable_sort_by_key(rects_2.begin(), rects_2.end(), rects_2_weights.begin(), sort_d4_x());

    double4* local_rects_2 = thrust::raw_pointer_cast(rects_2.data());
    ullong* local_rects_2_weights = thrust::raw_pointer_cast(rects_2_weights.data());

    /* END pre-processing */
    auto t_sort_end = get_current_time();
    t_sort = t_sort_end - t_begin;

    cudaDeviceSynchronize();

    ullong num_cells_1 = (mbr_arr_1->size / num_mbrs_per_cell) + 1;
    ullong last_mbr_pos_1 = mbr_arr_1->size - num_mbrs_per_cell * (num_cells_1 - 1);
    if(num_cells_1 == 1)
        last_mbr_pos_1 = mbr_arr_1->size;

    ullong num_cells_2 = (mbr_arr_2->size / num_mbrs_per_cell) + 1;
    ullong last_mbr_pos_2 = mbr_arr_2->size - num_mbrs_per_cell*(num_cells_2-1);
    if (num_cells_2 == 1)
        last_mbr_pos_2 = mbr_arr_2->size;

    //Results: cell MBRs
    double4 *cells_1;
    checkCudaErrors(cudaMalloc(&cells_1, num_cells_1 * sizeof(double4)));
    double4 *cells_2;
    checkCudaErrors(cudaMalloc(&cells_2, num_cells_2 * sizeof(double4)));

    ullong block_size = 512;
    ullong num_blocks_1 = 1 + (mbr_arr_1->size / block_size);
    ullong num_blocks_2 = 1 + (mbr_arr_2->size / block_size);

    //get mbrs of new cells
    candidate_grouper <<<num_blocks_1, block_size>>> (mbr_arr_1->size, local_rects_1, cells_1, num_cells_1, 
                                                last_mbr_pos_1, num_mbrs_per_cell);

    candidate_grouper <<<num_blocks_2, block_size>>> (mbr_arr_2->size, local_rects_2, cells_2, num_cells_2, 
                                                last_mbr_pos_2, num_mbrs_per_cell);
        
    ullong max_num_candidates = mbr_arr_1->size * mbr_arr_2->size > MAX_NUM_CANDIDATES ? 
                                MAX_NUM_CANDIDATES : mbr_arr_1->size * mbr_arr_2->size;

    thrust::device_vector<double> dev_candidates_m_x(max_num_candidates);
    thrust::device_vector<double> dev_candidates_m_y(max_num_candidates);
    thrust::device_vector<ullong> dev_candidates_weights(max_num_candidates);

    Candidates dev_candidates_init;
    dev_candidates_init.m_x = thrust::raw_pointer_cast(&dev_candidates_m_x[0]);
    dev_candidates_init.m_y = thrust::raw_pointer_cast(&dev_candidates_m_y[0]);
    dev_candidates_init.weights = thrust::raw_pointer_cast(&dev_candidates_weights[0]);

    Candidates *dev_candidates;
    checkCudaErrors(cudaMalloc((void **) &dev_candidates, sizeof(Candidates)));
    checkCudaErrors(cudaMemcpy(dev_candidates, &dev_candidates_init, sizeof(Candidates), cudaMemcpyHostToDevice));
 
    ullong *dev_num_candidates;
    checkCudaErrors(cudaMalloc(&dev_num_candidates, sizeof(ullong)));
    checkCudaErrors(cudaMemset(dev_num_candidates, 0, sizeof(ullong)));
 
    ullong *dev_global_weight;
    checkCudaErrors(cudaMalloc(&dev_global_weight, sizeof(ullong)));
    checkCudaErrors(cudaMemset(dev_global_weight, 0UL, sizeof(ullong)));

    cudaDeviceSynchronize();

    //generate candidates and their weights
    candidate_generater <<<num_cells_1, DEFAULT_THREADS_PER_BLOCK>>> (mbr_arr_1->size, local_rects_1, 
            local_rects_1_weights, last_mbr_pos_1, mbr_arr_2->size, local_rects_2, local_rects_2_weights, 
            last_mbr_pos_2, cells_1, num_cells_1, cells_2, num_cells_2, num_mbrs_per_cell, dev_num_candidates, 
            dev_candidates, dev_global_weight);

    cudaDeviceSynchronize();
    auto t_generate_end = get_current_time();

    checkCudaErrors(cudaFree(cells_1));
    checkCudaErrors(cudaFree(cells_2));

    ullong num_candidates;	
    checkCudaErrors(cudaMemcpy(&num_candidates, dev_num_candidates, sizeof(ullong), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(global_weight, dev_global_weight, sizeof(ullong), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(dev_candidates));
    checkCudaErrors(cudaFree(dev_num_candidates));
    checkCudaErrors(cudaFree(dev_global_weight));
    
    thrust::host_vector<double> x_h(num_candidates);
    thrust::host_vector<double> y_h(num_candidates);
    thrust::host_vector<double> w_h(num_candidates);
    thrust::copy_n(dev_candidates_m_x.begin(), num_candidates, x_h.begin());
    thrust::copy_n(dev_candidates_m_y.begin(), num_candidates, y_h.begin());
    thrust::copy_n(dev_candidates_weights.begin(), num_candidates, w_h.begin());

    // copy candidates into STL vectors
    std::vector<double> stl_vect_x(num_candidates);
    std::vector<double> stl_vect_y(num_candidates);
    std::vector<ullong> stl_vect_weights(num_candidates);

    thrust::copy(x_h.begin(), x_h.end(), stl_vect_x.begin());
    thrust::copy(y_h.begin(), y_h.end(), stl_vect_y.begin());
    thrust::copy(w_h.begin(), w_h.end(), stl_vect_weights.begin());

    candidates->m_x = (double*)malloc(num_candidates*sizeof(double));
    candidates->m_y = (double*)malloc(num_candidates*sizeof(double));
    candidates->weights = (ullong*)malloc(num_candidates*sizeof(ullong));

    std::copy(stl_vect_x.begin(), stl_vect_x.end(), candidates->m_x);
    std::copy(stl_vect_y.begin(), stl_vect_y.end(), candidates->m_y);
    std::copy(stl_vect_weights.begin(), stl_vect_weights.end(), candidates->weights);

    auto t_end = get_current_time();
    t_total = t_end - t_begin;
    double t_quad = t_generate_end - t_sort_end;
    
#ifdef DEBUG
    printf("CUDA_Utils_candidate_generate_wrapper total time %f, sorting %f, can %f. Num candidates %lu :: Global weight %lu\n", 
            t_total, t_sort, t_quad, num_candidates, *global_weight);
#endif //ifdef DEBUG

    return num_candidates;
}


/* Wrapper to quadtree partition */
ullong* CUDA_Utils_quadtree_partition_wrapper( Candidates *candidates, ullong num_candidates, MBR global_mbr, ullong global_weight, 
            ullong target_partitions, MBR_Array *partitioned_cells)
{
    //for (ullong i = 0; i < num_candidates; ++i)
        //printf("%llu: %f %f, %llu\n", i, candidates->m_x[i], candidates->m_y[i], candidates->weights[i]);
    // Allocate memory for points.
    thrust::device_vector<double> x_d0(candidates->m_x, candidates->m_x + num_candidates);
    thrust::device_vector<double> x_d1(num_candidates, 0.0);
    thrust::device_vector<double> y_d0(candidates->m_y, candidates->m_y + num_candidates);
    thrust::device_vector<double> y_d1(num_candidates, 0.0);
    thrust::device_vector<ullong> w_d0(candidates->weights, candidates->weights + num_candidates);
    thrust::device_vector<ullong> w_d1(num_candidates, 0ULL);

    // Host structures to analyze the device ones.
    Q_Candidates points_init[2];
    points_init[0].set(thrust::raw_pointer_cast(&x_d0[0]), thrust::raw_pointer_cast(&y_d0[0]), thrust::raw_pointer_cast(&w_d0[0]));
    points_init[1].set(thrust::raw_pointer_cast(&x_d1[0]), thrust::raw_pointer_cast(&y_d1[0]), thrust::raw_pointer_cast(&w_d1[0]));

    // Allocate memory to store points.
    Q_Candidates *points;
    checkCudaErrors(cudaMalloc((void **) &points, 2*sizeof(Q_Candidates)));
    checkCudaErrors(cudaMemcpy(points, points_init, 2*sizeof(Q_Candidates), cudaMemcpyHostToDevice));

    const uint max_depth  = MAX_QUADTREE_DEPTH;

    const ullong min_weight_per_node = (global_weight/target_partitions) * CELL_WEIGTH_COEF;

    // Allocate memory to store the tree.
    Quadtree_node root;
    root.set_range(0ULL, num_candidates);
    root.set_bounding_box(global_mbr.min_x, global_mbr.min_y, global_mbr.max_x, global_mbr.max_y);
    Quadtree_node *nodes;
    checkCudaErrors(cudaMalloc((void **) &nodes, sizeof(Quadtree_node)));
    checkCudaErrors(cudaMemcpy(nodes, &root, sizeof(Quadtree_node), cudaMemcpyHostToDevice));

    // We set the recursion limit for CDP to max_depth.
    cudaDeviceSetLimit(cudaLimitDevRuntimeSyncDepth, max_depth+1);
    cudaDeviceSetLimit(cudaLimitDevRuntimePendingLaunchCount, target_partitions*4);

    // Allocate memory for partition limits.
    ullong *total_partitions;
    double4 *outputs;
    ulonglong2 *output_pos;
    ullong *output_weights;

    checkCudaErrors(cudaMalloc((void **) &total_partitions, sizeof(ullong)));
    checkCudaErrors(cudaMemset(total_partitions, 0, sizeof(ullong)));

    checkCudaErrors(cudaMalloc((void **) &outputs, target_partitions * 2 * sizeof(double4)));
    checkCudaErrors(cudaMalloc((void **) &output_pos, target_partitions * 2 * sizeof(ulonglong2)));
    checkCudaErrors(cudaMalloc((void **) &output_weights, target_partitions * 2 * sizeof(ullong)));
 
    // Build the quadtree.
    Parameters params(max_depth, min_weight_per_node, target_partitions);
    //printf("Launching CUDA kernel to build the quadtree. Warp size %d \n", warp_size);

    // QUADTREE_THREADS_PER_BLOCK = 256; Do not use less than 128 threads.
    const int NUM_WARPS_PER_BLOCK = QUADTREE_THREADS_PER_BLOCK / WARP_SIZE;
    const size_t smem_size = 4*NUM_WARPS_PER_BLOCK*sizeof(int);

    build_quadtree_kernel<<<1, QUADTREE_THREADS_PER_BLOCK, smem_size>>>(nodes, total_partitions, outputs, output_pos, output_weights, points, params);
    checkCudaErrors(cudaGetLastError());

    ullong host_total_partitions;	
    checkCudaErrors(cudaMemcpy(&host_total_partitions, total_partitions, sizeof(ullong), cudaMemcpyDeviceToHost));

    partitioned_cells->size = host_total_partitions;
    partitioned_cells->mbrs = (MBR*)malloc(host_total_partitions * sizeof(MBR));
    partitioned_cells->weights = (ullong*)malloc(host_total_partitions * sizeof(ullong));

    //ulonglong2* host_temp_pos;
    //host_temp_pos = (ulonglong2*)malloc(host_total_partitions * sizeof(ulonglong2));

    checkCudaErrors(cudaMemcpy(partitioned_cells->mbrs, outputs, host_total_partitions * sizeof(double4), cudaMemcpyDeviceToHost));
    //checkCudaErrors(cudaMemcpy(host_temp_pos, output_pos, host_total_partitions * sizeof(ulonglong2), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(partitioned_cells->weights, output_weights, host_total_partitions * sizeof(ullong), cudaMemcpyDeviceToHost));

    // Copy points to CPU.
    //thrust::host_vector<double> x_h(x_d0);
    //thrust::host_vector<double> y_h(y_d0);
    //thrust::host_vector<ullong> w_h(w_d0);

    //return (ullong*)host_temp_pos;
    //TODO: return null to save memory copy. Replace NULL If needed in future.
    return NULL;
}


/* CUDA sorting multiple arryas using thrust gather*/

/*
    // allocate space for the output
    thrust::device_vector<double4> sortedMBRs(local_total_partitions);
    thrust::device_vector<int2> sortedPos(local_total_partitions);

    // initialize indices vector to [0,1,2,..]
    thrust::counting_iterator<int> iter(0);
    thrust::device_vector<int> indices(local_total_partitions);
    thrust::copy(iter, iter + indices.size(), indices.begin());
    //printf("%d \n", 2);

    // first sort the keys and indices by the keys
    thrust::sort_by_key(dev_weight.begin(), dev_weight.end(), indices.begin());

    // Now reorder the ID arrays using the sorted indices
    thrust::gather(indices.begin(), indices.end(), dev_mbrs.begin(), sortedMBRs.begin());
    thrust::gather(indices.begin(), indices.end(), dev_pos.begin(), sortedPos.begin());
    //printf("%d \n", 3);
    
    ullong* local_temp_weight;
    local_temp_weight = (ullong*)malloc(local_total_partitions * sizeof(ullong));
    ullong* dummy_weight = thrust::raw_pointer_cast(dev_weight.data());
    checkCudaErrors(cudaMemcpy(local_temp_weight, dummy_weight, local_total_partitions * sizeof(ullong), cudaMemcpyDeviceToHost));

    double4* local_temp_mbrs;
    local_temp_mbrs = (double4*)malloc(local_total_partitions * sizeof(double4));
    double4* dummy_mbrs = thrust::raw_pointer_cast(sortedMBRs.data());
    checkCudaErrors(cudaMemcpy(local_temp_mbrs, dummy_mbrs, local_total_partitions * sizeof(double4), cudaMemcpyDeviceToHost));

    int2* local_temp_pos;
    local_temp_pos = (int2*)malloc(local_total_partitions * sizeof(int2));
    int2* dummy_pos = thrust::raw_pointer_cast(sortedPos.data());
    checkCudaErrors(cudaMemcpy(local_temp_pos, dummy_pos, local_total_partitions * sizeof(int2), cudaMemcpyDeviceToHost));
*/
