#ifndef __CUDA_QUADTREE_CUH_INCLUDED__
#define __CUDA_QUADTREE_CUH_INCLUDED__
#include <global_var.h>

/* External libs */
#include <sys/time.h>
#include <helper_cuda.h>

/* Internal libs */

/* Classes */
////////////////////////////////////////////////////////////////////////////////
// A structure of 2D points (structure of arrays).
////////////////////////////////////////////////////////////////////////////////
class Q_Candidates
{
    private:
        //candidate* candidates;
        double *m_x;
        double *m_y;
        ullong *m_weight;

    public:
        // Constructor.
        __host__ __device__ Q_Candidates() : m_x(NULL), m_y(NULL), m_weight(NULL) {}

        // Constructor.
        //__host__ __device__ Candidates(candidate* in_candidates) : candidates(in_candidates) {}
        __host__ __device__ Q_Candidates(double *x, double *y, ullong *weight) : m_x(x), m_y(y), m_weight(weight) {}

        // Get a candidate.
        //__host__ __device__ __forceinline__ candidate get_candidate(int idx) const
        //{
            //return candidates[idx];
        //}

        // Get a point.
        __host__ __device__ __forceinline__ double2 get_point(ullong idx) const
        {
            return make_double2(m_x[idx], m_y[idx]);
        }

        __host__ __device__ __forceinline__ ullong get_weight(ullong idx) const
        {
            return m_weight[idx];
        }

        // Set a point.
        __host__ __device__ __forceinline__ void set_point(ullong idx, const double2 &p, const ullong &w)
        {
            m_x[idx] = p.x;
            m_y[idx] = p.y;
            m_weight[idx] = w;
        }

        // Set the candidates.
        __host__ __device__ __forceinline__ void set(double *x, double *y, ullong * weight)
        {
            m_x = x;
            m_y = y;
            m_weight = weight;
        }
};


////////////////////////////////////////////////////////////////////////////////
// A 2D bounding box
////////////////////////////////////////////////////////////////////////////////

class Bounding_box
{
    private:
        // Extreme points of the bounding box.
        double2 m_p_min;
        double2 m_p_max;

    public:
        // Constructor. Create a unit box.
        __host__ __device__ Bounding_box()
        {
            m_p_min = make_double2(0.0f, 0.0f);
            m_p_max = make_double2(1.0f, 1.0f);
        }

        // Compute the center of the bounding-box.
        __host__ __device__ void compute_center(double2 &center) const
        {
            center.x = 0.5f * (m_p_min.x + m_p_max.x);
            center.y = 0.5f * (m_p_min.y + m_p_max.y);
        }

        // The points of the box.
        __host__ __device__ __forceinline__ const double2 &get_max() const
        {
            return m_p_max;
        }

        __host__ __device__ __forceinline__ const double2 &get_min() const
        {
            return m_p_min;
        }

        // Does a box contain a point.
        __host__ __device__ bool contains(const double2 &p) const
        {
            return p.x >= m_p_min.x && p.x < m_p_max.x && p.y >= m_p_min.y && p.y < m_p_max.y;
        }

        __host__ __device__ bool contains(const double &x, const double &y) const
        {
            return x >= m_p_min.x && x < m_p_max.x && y >= m_p_min.y && y < m_p_max.y;
        }

        // Define the bounding box.
        __host__ __device__ void set(double min_x, double min_y, double max_x, double max_y)
        {
            m_p_min.x = min_x;
            m_p_min.y = min_y;
            m_p_max.x = max_x;
            m_p_max.y = max_y;
        }
};

////////////////////////////////////////////////////////////////////////////////
// A node of a quadree.
////////////////////////////////////////////////////////////////////////////////
class Quadtree_node
{
    private:
        // The identifier of the node.
        //int m_id;
        // The bounding box of the tree.
        Bounding_box m_bounding_box;
        // The range of points.
        ullong m_begin, m_end;


    public:
        // Constructor.
        __host__ __device__ Quadtree_node() : m_begin(0), m_end(0) // m_id(0), 
        {}

        // The ID of a node at its level.
        //__host__ __device__ int id() const
        //{
        //    return m_id;
        //}

        // The ID of a node at its level.
        //__host__ __device__ void set_id(int new_id)
        //{
         //   m_id = new_id;
        //}

        // The bounding box.
        __host__ __device__ __forceinline__ const Bounding_box &bounding_box() const
        {
            return m_bounding_box;
        }

        // Set the bounding box.
        __host__ __device__ __forceinline__ void set_bounding_box(double min_x, double min_y, double max_x, double max_y)
        {
            m_bounding_box.set(min_x, min_y, max_x, max_y);
        }

        // The number of candidates in the tree.
        __host__ __device__ __forceinline__ ullong num_candidates() const
        {
            return m_end - m_begin;
        }

        // The range of candidates in the tree.
        __host__ __device__ __forceinline__ ullong candidates_begin() const
        {
            return m_begin;
        }

        __host__ __device__ __forceinline__ ullong candidates_end() const
        {
            return m_end;
        }

        // Define the range for that node.
        __host__ __device__ __forceinline__ void set_range(ullong begin, ullong end)
        {
            m_begin = begin;
            m_end = end;
        }
};


////////////////////////////////////////////////////////////////////////////////
// Algorithm parameters.
////////////////////////////////////////////////////////////////////////////////
struct Parameters
{
    // Choose the right set of candidates to use as in/out.
    int candidate_selector;
    // The number of nodes at a given level (2^k for level k).
    //int num_nodes_at_this_level;
    // The recursion depth.
    int depth;
    // The max value for depth.
    const int max_depth;
    // The minimum weight of a node to stop recursion.
    const ullong min_weight_per_node;

    // The target num partitions.
    const ullong target_partitions;

    // Constructor set to default values.
    __host__ __device__ Parameters(int max_depth, ullong min_weight_per_node, ullong target_partitions) :
        candidate_selector(0),
        //num_nodes_at_this_level(1),
        depth(0),
        max_depth(max_depth),
        min_weight_per_node(min_weight_per_node),
        target_partitions(target_partitions)
    {}

    // Copy constructor. Changes the values for next iteration.
    __host__ __device__ Parameters(const Parameters &params, bool) :
        candidate_selector((params.candidate_selector+1) % 2),
        //num_nodes_at_this_level(4*params.num_nodes_at_this_level),
        depth(params.depth+1),
        max_depth(params.max_depth),
        min_weight_per_node(params.min_weight_per_node),
        target_partitions(params.target_partitions)
    {}
};

/* Global Functions */
__global__ void build_quadtree_kernel(Quadtree_node *nodes, ullong *total_partitions, double4 *outputs, ulonglong2 *output_pos, ullong *output_weight, 
    Q_Candidates *candidates, Parameters params);

#endif //ndef __CUDA_QUADTREE_CUH_INCLUDED__
