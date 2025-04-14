#include<cuda_quadtree.cuh>
#include <cooperative_groups.h>
namespace cg = cooperative_groups;

////////////////////////////////////////////////////////////////////////////////
// Build a quadtree on the GPU. Use CUDA Dynamic Parallelism.
// This algorithm takes points weights into consideration.
// A candidate is a point with weight.
// This algorithm will transfer to regular one if every candidate with weight = 1.
//
// The algorithm works as follows. The host (CPU) launches one block of
// QUADTREE_THREADS_PER_BLOCK threads. That block will do the following steps:
//
// 1- Check the total weight of candidates and its depth.
//
// We impose a maximum depth to the tree and a minimum weight of candidates per
// node. If the maximum depth is exceeded or the minimum weight of candidates is
// reached. The threads in the block exit.
//
// Before exiting, they perform a buffer swap if it is needed. Indeed, the
// algorithm uses two buffers to permute the points and make sure they are
// properly distributed in the quadtree. By design we want all points to be
// in the first buffer of points at the end of the algorithm. It is the reason
// why we may have to swap the buffer before leavin (if the points are in the
// 2nd buffer).
//
// 2- Computing the weight of candidates in each child.
//
// If the depth is not too high and the weight of candidates is sufficient, the
// block has to dispatch the candidates into four geometrical buckets: Its
// children. For that purpose, we compute the center of the bounding box and
// count the weights of candidates in each quadrant.
//
// The set of candidates is divided into sections. Each section is given to a
// warp of threads (32 threads).
//
// 3- Scan the warps' results to know the "global" numbers.
//
// Warps work independently from each other. At the end, each warp knows the
// weight of candidates in its section. To know the weight for the block, the
// block has to run a scan/reduce at the block level. It's a traditional
// approach. The implementation in that sample is not as optimized as what
// could be found in fast radix sorts, for example, but it relies on the same
// idea.
//
// 4- Move candidates.
//
// Now that the block knows how many candidates go in each of its 4 children, it
// remains to dispatch the candidates. It is straightforward.
//
// 5- Launch new blocks.
//
// The block launches four new blocks: One per children. Each of the four blocks
// will apply the same algorithm.
////////////////////////////////////////////////////////////////////////////////

__global__ void build_quadtree_kernel(Quadtree_node *nodes, ullong *total_partitions, double4 *outputs, 
        ulonglong2 *output_pos, ullong *output_weight, Q_Candidates *candidates, Parameters params)
{
    // Handle to thread block group
    cg::thread_block cta = cg::this_thread_block();

    // The number of warps in a block.   
    const int warpSize = 32;
    const int NUM_WARPS_PER_BLOCK = QUADTREE_THREADS_PER_BLOCK / warpSize;

    // Shared memory to store the number of points.
    extern __shared__ int smem[];

    // s_num_pts[4][NUM_WARPS_PER_BLOCK];
    // Addresses of shared memory.
    volatile int *s_num_pts[4];

    for (int i = 0 ; i < 4 ; ++i)
        s_num_pts[i] = (volatile int *) &smem[i*NUM_WARPS_PER_BLOCK];

    // Compute the coordinates of the threads in the block.
    const int warp_id = threadIdx.x / warpSize;
    const int lane_id = threadIdx.x % warpSize;

    // Mask for compaction.
    int lane_mask_lt = (1 << lane_id) - 1; // Same as: asm( "mov.u32 %0, %%lanemask_lt;" : "=r"(lane_mask_lt) );

    //threadsPerBlock  = blockDim.x * blockDim.y
    //threadNumInBlock = threadIdx.x + blockDim.x * threadIdx.y (alternatively: threadIdx.y + blockDim.y * threadIdx.x)
    //blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y  (alternatively: blockIdx.y  + gridDim.y  * blockIdx.x)

    int blockNumInGrid   = blockIdx.x  + gridDim.x  * blockIdx.y;
    // The current node.
    Quadtree_node &node = nodes[blockNumInGrid];

    // The number of points in the node.
    ullong num_candidates = node.num_candidates();

    // The weight of the node.
    ullong temp_weight  = 0ULL;  
    //printf("%d\n",params.min_weight_per_node);
    {
        ullong it = node.candidates_begin(), end = node.candidates_end();
        //for (it += threadIdx.x ; it < end ; it += QUADTREE_THREADS_PER_BLOCK)
        for (; it < end ; it++)
            temp_weight += candidates[params.candidate_selector].get_weight(it);


    }
    __syncthreads();

    //
    // 1- Check the weight of candidates and its depth.

    // Stop the recursion here. Make sure candidates[0] contains all the candidates.
    if (params.depth >= params.max_depth || temp_weight <= params.min_weight_per_node ||   *total_partitions >= params.target_partitions 
    || num_candidates <= 100)
    {
        if (params.candidate_selector == 1)
        {
            ullong it = node.candidates_begin(), end = node.candidates_end();

            for (it += threadIdx.x ; it < end ; it += QUADTREE_THREADS_PER_BLOCK)
                if (it < end)
                    candidates[0].set_point(it, candidates[1].get_point(it), candidates[1].get_weight(it));
        }

        //Thread 0 records the current Bounding_box
        
        if (threadIdx.x == QUADTREE_THREADS_PER_BLOCK-1 && temp_weight != 0UL)
        {
            int write_index = atomicAdd(total_partitions, 1);

            const Bounding_box &bbox = node.bounding_box();
            outputs[write_index] = make_double4(bbox.get_min().x, bbox.get_min().y, bbox.get_max().x, bbox.get_max().y);
            output_pos[write_index] = make_ulonglong2(node.candidates_begin(), node.candidates_end());
            output_weight[write_index] = temp_weight;

        }
        return;
    }

    // Compute the center of the bounding box of the points.
    const Bounding_box &bbox = node.bounding_box();
    double2 center;
    bbox.compute_center(center);

    // Find how many points to give to each warp.
    int num_candidates_per_warp = max(warpSize, (int)((num_candidates + NUM_WARPS_PER_BLOCK-1) / NUM_WARPS_PER_BLOCK));

    // Each warp of threads will compute the number of points to move to each quadrant.
    ullong range_begin = node.candidates_begin() + warp_id * num_candidates_per_warp;
    ullong range_end   = min(range_begin + num_candidates_per_warp, node.candidates_end());
    int warp_cnts[4] = {0, 0, 0, 0};
    //printf("%d\n",4);
    //
    // 2- Count the number of points in each child.
    //

    // Reset the counts of points per child.
    if (lane_id == 0)
    {
        s_num_pts[0][warp_id] = 0;
        s_num_pts[1][warp_id] = 0;
        s_num_pts[2][warp_id] = 0;
        s_num_pts[3][warp_id] = 0;
    }

    // Input points.
    const Q_Candidates &in_candidates = candidates[params.candidate_selector];

    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cta);
    //printf("%d\n",5);
    // Compute the number of points.
    for (ullong range_it = range_begin + tile32.thread_rank();
        tile32.any(range_it < range_end); range_it += warpSize) {
        // Is it still an active thread?
        bool is_active = range_it < range_end;

        // Load the coordinates of the point.
        double2 p =
            is_active ? in_candidates.get_point(range_it) : make_double2(0.0f, 0.0f);

        // Count top-left points.
        int num_pts =
            __popc(tile32.ballot(is_active && p.x < center.x && p.y >= center.y));
        warp_cnts[0] += tile32.shfl(num_pts, 0);

        // Count top-right points.
        num_pts =
            __popc(tile32.ballot(is_active && p.x >= center.x && p.y >= center.y));
        warp_cnts[1] += tile32.shfl(num_pts, 0);

        // Count bottom-left points.
        num_pts =
            __popc(tile32.ballot(is_active && p.x < center.x && p.y < center.y));
        warp_cnts[2] += tile32.shfl(num_pts, 0);

        // Count bottom-right points.
        num_pts =
            __popc(tile32.ballot(is_active && p.x >= center.x && p.y < center.y));
        warp_cnts[3] += tile32.shfl(num_pts, 0);
    }


    __syncthreads();
    if (tile32.thread_rank() == 0) {
        s_num_pts[0][warp_id] = warp_cnts[0];
        s_num_pts[1][warp_id] = warp_cnts[1];
        s_num_pts[2][warp_id] = warp_cnts[2];
        s_num_pts[3][warp_id] = warp_cnts[3];
    }

    // Make sure warps have finished counting.
    __syncthreads();

    //
    // 3- Scan the warps' results to know the "global" numbers.
    //

    // First 4 warps scan the numbers of points per child (inclusive scan).
    if (warp_id < 4)
    {
        int num_pts = lane_id < NUM_WARPS_PER_BLOCK ? s_num_pts[warp_id][tile32.thread_rank()] : 0;
#pragma unroll

        for (int offset = 1 ; offset < NUM_WARPS_PER_BLOCK ; offset *= 2)
        {
            int n = tile32.shfl_up(num_pts, offset);

            if (tile32.thread_rank() >= offset)
            num_pts += n;
        }

        if (tile32.thread_rank() < NUM_WARPS_PER_BLOCK)
        {
            s_num_pts[warp_id][tile32.thread_rank()] = num_pts;
        }
            
    }

    __syncthreads();

    // Compute global offsets.
    if (warp_id == 0)
    {
        int sum = s_num_pts[0][NUM_WARPS_PER_BLOCK-1];

        for (int row = 1 ; row < 4 ; ++row)
        {
            int tmp = s_num_pts[row][NUM_WARPS_PER_BLOCK-1];

            if (tile32.thread_rank() < NUM_WARPS_PER_BLOCK)
                s_num_pts[row][tile32.thread_rank()] += sum;

            cg::sync(tile32);
            sum += tmp;
        }
    }

    __syncthreads();
    int val = 0;
    // Make the scan exclusive.
    if (threadIdx.x < 4*NUM_WARPS_PER_BLOCK)
    {
        val = threadIdx.x == 0 ? 0 : smem[threadIdx.x-1];
        val += node.candidates_begin();
    }

    __syncthreads();

    if (threadIdx.x < 4*NUM_WARPS_PER_BLOCK)
    {
        smem[threadIdx.x] = val;
    }

    __syncthreads();

    //
    // 4- Move candidates.
    //

    // Output points.
    Q_Candidates &out_candidates = candidates[(params.candidate_selector+1) % 2];

    warp_cnts[0] = s_num_pts[0][warp_id];
    warp_cnts[1] = s_num_pts[1][warp_id];
    warp_cnts[2] = s_num_pts[2][warp_id];
    warp_cnts[3] = s_num_pts[3][warp_id];


    // Reorder points.
    for (ullong range_it = range_begin + tile32.thread_rank();
            tile32.any(range_it < range_end); range_it += warpSize) {
        // Is it still an active thread?
        bool is_active = range_it < range_end;

        // Load the coordinates of the point.
        double2 p =
            is_active ? in_candidates.get_point(range_it) : make_double2(0.0f, 0.0f);

        ullong w = is_active ? in_candidates.get_weight(range_it) : 0UL;

        // Count top-left points.
        bool pred = is_active && p.x < center.x && p.y >= center.y;
        int vote = tile32.ballot(pred);
        int dest = warp_cnts[0] + __popc(vote & lane_mask_lt);

        if (pred) out_candidates.set_point(dest, p, w);

        warp_cnts[0] += tile32.shfl(__popc(vote), 0);

        // Count top-right points.
        pred = is_active && p.x >= center.x && p.y >= center.y;
        vote = tile32.ballot(pred);
        dest = warp_cnts[1] + __popc(vote & lane_mask_lt);

        if (pred) out_candidates.set_point(dest, p, w);

        warp_cnts[1] += tile32.shfl(__popc(vote), 0);

        // Count bottom-left points.
        pred = is_active && p.x < center.x && p.y < center.y;
        vote = tile32.ballot(pred);
        dest = warp_cnts[2] + __popc(vote & lane_mask_lt);

        if (pred) out_candidates.set_point(dest, p, w);

        warp_cnts[2] += tile32.shfl(__popc(vote), 0);

        // Count bottom-right points.
        pred = is_active && p.x >= center.x && p.y < center.y;
        vote = tile32.ballot(pred);
        dest = warp_cnts[3] + __popc(vote & lane_mask_lt);

        if (pred) out_candidates.set_point(dest, p, w);

        warp_cnts[3] += tile32.shfl(__popc(vote), 0);
    }
    __syncthreads();

    if (tile32.thread_rank() == 0) {
        s_num_pts[0][warp_id] = warp_cnts[0];
        s_num_pts[1][warp_id] = warp_cnts[1];
        s_num_pts[2][warp_id] = warp_cnts[2];
        s_num_pts[3][warp_id] = warp_cnts[3];
      }

      __syncthreads();
    //
    // 5- Launch new blocks.
    //

    // The last thread launches new blocks.
    if (threadIdx.x == QUADTREE_THREADS_PER_BLOCK-1)
    {
        // The children.
        Quadtree_node *children;
        cudaMalloc((void **) &children, 4*sizeof(Quadtree_node));

        // Points of the bounding-box.
        const double2 &p_min = bbox.get_min();
        const double2 &p_max = bbox.get_max();

        // Set the bounding boxes of the children.
        children[0].set_bounding_box(p_min.x , center.y, center.x, p_max.y);    // Top-left.
        children[1].set_bounding_box(center.x, center.y, p_max.x , p_max.y);    // Top-right.
        children[2].set_bounding_box(p_min.x , p_min.y , center.x, center.y);   // Bottom-left.
        children[3].set_bounding_box(center.x, p_min.y , p_max.x , center.y);   // Bottom-right.

        // Set the ranges of the children.
        children[0].set_range(node.candidates_begin(),   s_num_pts[0][warp_id]);
        children[1].set_range(s_num_pts[0][warp_id], s_num_pts[1][warp_id]);
        children[2].set_range(s_num_pts[1][warp_id], s_num_pts[2][warp_id]);
        children[3].set_range(s_num_pts[2][warp_id], s_num_pts[3][warp_id]);

        // Launch 4 children.
        build_quadtree_kernel<<<4, QUADTREE_THREADS_PER_BLOCK, 4 *NUM_WARPS_PER_BLOCK *sizeof(int)>>>
            (children, total_partitions, outputs, output_pos, output_weight, candidates, Parameters(params, true));
    }
}

