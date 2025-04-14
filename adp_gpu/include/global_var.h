#ifndef __GSJ_GLOBAL_VAR_H_INCLUDE__
#define __GSJ_GLOBAL_VAR_H_INCLUDE__

/* Compiler instructions*/
#define USE_UNSTABLE_GEOS_CPP_API

/* Global Variables */
#define VERBOSE 1
#define DEBUG 1

/* Types */
typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned long long int ullong;
#ifndef ULLONG_MAX
#define ULLONG_MAX 18446744073709551615ULL
#endif //ifndef ULLONG_MAX

typedef unsigned __int128 uint128_t;

/* GPU */
#define WAKEUP_CYCLE 1024*1024
#define WAKEUP_VAR 555
#define WAKEUP_DUMMY 777
#define DEFAULT_THREADS_PER_BLOCK 256
#define QUADTREE_THREADS_PER_BLOCK 256
#define DEFAULT_MBRS_PER_BLOCK 256ULL
#define MAX_NUM_CANDIDATES 1024ULL*1024*75
#define MAX_PARTITIONS_SINGL_RUN 1024
#define WARP_SIZE 32
#define CELL_WEIGTH_COEF 3 

#define MAX_QUADTREE_DEPTH 17

#ifndef W_POINT_STRUCT
#define W_POINT_STRUCT
typedef struct                      
{
    double x;            
    double y;       
} M_point;
#endif //ndef M_POINT_STRUCT

#ifndef CANDIDATES_STRUCT
#define CANDIDATES_STRUCT
typedef struct                      
{
    double *m_x; 
    double *m_y; 
    ullong *weights;           
} Candidates;
#endif //ndef CANDIDATES_STRUCT

#ifndef GPU_MBRS_STRUCT
#define GPU_MBRS_STRUCT
typedef struct                      
{
    double min_x;            
    double min_y;
    double max_x;            
    double max_y;        
} MBR;
typedef struct                      
{
    ullong size;
    MBR *mbrs;
    ullong *weights;           
} MBR_Array;
#endif //ndef GPU_MBRS_STRUCT


#endif //ndef __GSJ_GLOBAL_VAR_H_INCLUDE__
