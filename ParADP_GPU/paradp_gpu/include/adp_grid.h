#ifndef __ADP_GRID_H_INCLUDED__
#define __ADP_GRID_H_INCLUDED__

#include <global_var.h>

#include <stdio.h>
#include <stdlib.h>

#include <list>
#include <cmath>
#include <cstring>
#include <limits.h>
#include <unordered_map>

#include <geos/geom/Envelope.h>
#include <geos/geom/Geometry.h>
#include <geos/index/strtree/STRtree.h>
#include <spdlog/spdlog.h>  // Log lib
#include <spdlog/cfg/env.h> // support for loading levels from the environment variable

#include <mpi.h>

#include <grid_utils.h>
#include <cuda_host_interface.cuh>
//#include<cuda_utils.cuh>

#define MAX_SIZE LLONG_MAX

class Adp_Grid
{
public:
	Adp_Grid(uint num_partitions, uint num_threads, const geos::geom::Envelope *universe, std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs_1,
			 std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs_2)
	{

#ifdef DEBUG
    spdlog::info("Enter function {0:s}", "Adp_Grid constructor");
#endif
		int num_world_nodes, my_world_rank;

		MPI_Comm_size(MPI_COMM_WORLD, &num_world_nodes);
		MPI_Comm_rank(MPI_COMM_WORLD, &my_world_rank);

		/* init class private variables */
		global_env = NULL;
		vect_env_weigt_1 = new std::vector<std::pair<const geos::geom::Envelope *, ulong> *>();//(list_envs_1->begin(), list_envs_1->end());
		vect_env_weigt_2 = new std::vector<std::pair<const geos::geom::Envelope *, ulong> *>();//(list_envs_2->begin(), list_envs_2->end());
		envs_each_stripe_1 = NULL;
		envs_each_stripe_2 = NULL;
		stripe_cells = NULL;

#ifdef DEBUG
    spdlog::info("{0:s}", "Init finished");
#endif

		// Get envs with weights //
		for (std::list<std::pair<geos::geom::Envelope *, ulong> >::iterator itr = list_envs_1->begin(); itr != list_envs_1->end(); ++itr)
		{
			std::pair<geos::geom::Envelope *, ulong> temp_env_wieght = *itr;
			vect_env_weigt_1->push_back(new std::pair<const geos::geom::Envelope *, ulong>(
				temp_env_wieght.first, temp_env_wieght.second));
		}

		for (std::list<std::pair<geos::geom::Envelope *, ulong> >::iterator itr = list_envs_2->begin(); itr != list_envs_2->end(); ++itr)
		{
			std::pair<geos::geom::Envelope *, ulong> temp_env_wieght = *itr;
			vect_env_weigt_2->push_back(new std::pair<const geos::geom::Envelope *, ulong>(
				temp_env_wieght.first, temp_env_wieght.second));
		}

#ifdef DEBUG
    //spdlog::info("list_envs_1[0]:    {0:s}, vect_env_weigt_1[0]:    {1:s}", list_envs_1->front().first->toString(), vect_env_weigt_1->front()->first->toString());
	spdlog::info("{0:s}", "Get envs finished");
#endif

		/* Build index on layer 1 */
		build_rtree_index(list_envs_1);

		/* synchronize global MBR across all MPI processes */
		sync_universe(universe);

#ifdef DEBUG
    spdlog::info("{0:s}", "Index and sync universe finished");
#endif
		/* Parallel stripping based on layer 1, output in stripe_cells */
		parallel_sorting_by_regular_sampling(list_envs_1);

		/* init for next step */
		envs_each_stripe_1 = new std::vector<std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *>(num_world_nodes);
		envs_each_stripe_2 = new std::vector<std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *>(num_world_nodes);

		for (int i = 0; i < num_world_nodes; ++i)
		{
			envs_each_stripe_1->at(i) = new std::vector<std::pair<const geos::geom::Envelope *, ulong> *>();
			envs_each_stripe_2->at(i) = new std::vector<std::pair<const geos::geom::Envelope *, ulong> *>();
		}



#ifdef DEBUG
    spdlog::info("{0:s}", "Stripe finished");
#endif

		clustering_envs(num_threads, stripe_cells, vect_env_weigt_1, envs_each_stripe_1);
		clustering_envs(num_threads, stripe_cells, vect_env_weigt_2, envs_each_stripe_2);
		
		for (int i = 0; i < num_world_nodes; ++i)
		{
			 printf("Size of strip P%d, Strip number:%d for dataset 1:%lu\n", my_world_rank,i,envs_each_stripe_1->at(i)->size());
			 printf("Size of strip P%d, Strip number:%d for dataset 2:%lu\n", my_world_rank,i,envs_each_stripe_2->at(i)->size());
		}

		// global_vect_env_weigt_1 = new std::vector<std::pair<const geos::geom::Envelope*, ulong>*>();
		// collect_data(*global_vect_env_weigt_1, *vect_env_weigt_1);

		// global_vect_env_weigt_2 = new std::vector<std::pair<const geos::geom::Envelope*, ulong>*>();
		// collect_data(*global_vect_env_weigt_2, *vect_env_weigt_2);

		// printf("Size of dataset D1:%lu\n", my_world_rank,global_vect_env_weigt_1->size());
		// printf("Size of dataset D2:%lu\n",global_vect_env_weigt_2->size());

#ifdef MEMORY_ENHANCMENT
		//Enhancement: Free up the vect_env_weigt_1 & vect_env_weigt_2 since it never be used later!
		if (vect_env_weigt_1 != NULL)
			delete vect_env_weigt_1;

		if (vect_env_weigt_2 != NULL)
			delete vect_env_weigt_2;
#endif

#ifdef DEBUG
    spdlog::info("{0:s}", "Cluster finished");
#endif

		

		std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *envs_shuffled_1 =
			new std::vector<std::pair<const geos::geom::Envelope *, ulong> *>();
		std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *envs_shuffled_2 =
			new std::vector<std::pair<const geos::geom::Envelope *, ulong> *>();

		shuffle_candidates(envs_shuffled_1, envs_each_stripe_1, 0);
		shuffle_candidates(envs_shuffled_2, envs_each_stripe_2, 1);


		// sort_and_collect_dataset1(*envs_shuffled_1, *vect_env_weigt_1);

		printf("Size of strip P%d for dataset 1:%lu\n", my_world_rank,envs_shuffled_1->size());
		printf("Size of strip P%d for dataset 2:%lu\n", my_world_rank,envs_shuffled_2->size());

#ifdef MEMORY_ENHANCMENT
		//Enhancement: Free up the envs_each_stripe_1 & envs_each_stripe_1 since it never be used later!
		// if (vect_env_weigt_1 != NULL)
		// 	delete vect_env_weigt_1;

		if (envs_each_stripe_1 != NULL)
			delete envs_each_stripe_1;

		if (envs_each_stripe_2 != NULL)
			delete envs_each_stripe_2;
#endif

#ifdef DEBUG
    spdlog::info("{0:s}", "Shuffle finished");
#endif


// TODO: Here I can adapt the gpu part since I have all envs that releated to local stripe.

	MBR_Array mbr_arr_1;
	Util_Parse_vector_envs_to_mbr_arr(envs_shuffled_1, mbr_arr_1);

	printf("Util_Parse_vector_envs_to_mbr_arr for dataset 1 completed\n");

	MBR_Array mbr_arr_2;
	Util_Parse_vector_envs_to_mbr_arr(envs_shuffled_2, mbr_arr_2);
	// Util_Parse_vector_envs_to_mbr_arr(global_vect_env_weigt_2, mbr_arr_2);

	printf("Util_Parse_vector_envs_to_mbr_arr for dataset 2\n");

	printf("Going into wakeup wrapper\n");
	GPU_Utils_wakeup_wrapper(WAKEUP_CYCLE, my_world_rank);

	ullong num_candidates = 0ULL;
#ifdef DEBUG
	spdlog::info("{0:s}", "Begin CUDA ADP");
#endif
	//MBR global_mbr = {-180.0, -90.0, 180.0, 90.0};
	MBR global_mbr = {stripe_cells->at(my_world_rank)->getMinX(), 
	 				stripe_cells->at(my_world_rank)->getMinY(),
	 				stripe_cells->at(my_world_rank)->getMaxX(),
	 				stripe_cells->at(my_world_rank)->getMaxY()};
	
	MBR_Array partition_results;

	printf("GPU works starting NOW");

	//num_candidates = CUDA_Utils_adp_wrapper(&mbr_arr_1, &mbr_arr_2, temp_mbr, num_partitions, &partition_results);
	num_candidates = adp_gpu_wrapper(&mbr_arr_1, &mbr_arr_2, global_mbr, num_partitions, &partition_results);

#ifdef DEBUG
	spdlog::info("{0:s}", "End CUDA ADP");
#endif
#ifdef MEMORY_ENHANCMENT
		//Enhancement: Free up the stripe_cells since it never be used later!
		if (stripe_cells != NULL)
			delete stripe_cells;

		//Enhancement: Free up the envs_shuffled_1 & envs_shuffled_2 since it never be used later!
		if (envs_shuffled_1 != NULL)
			delete envs_shuffled_1;

		// if (envs_shuffled_2 != NULL)
		// 	delete envs_shuffled_2;

		if (global_vect_env_weigt_2 != NULL)
			delete global_vect_env_weigt_2;

#endif

#ifdef DEBUG
    spdlog::info("Leave function {0:s}", "Adp_Grid constructor");
#endif
	}

	~Adp_Grid()
	{
		if (global_env != NULL)
			delete global_env;

		if (vect_env_weigt_1 != NULL)
			delete vect_env_weigt_1;

		if (vect_env_weigt_2 != NULL)
			delete vect_env_weigt_2;

		if (envs_each_stripe_1 != NULL)
			delete envs_each_stripe_1;

		if (envs_each_stripe_2 != NULL)
			delete envs_each_stripe_2;

		if (stripe_cells != NULL)
			delete stripe_cells;

	}

private:
	geos::index::strtree::STRtree index;
	geos::geom::Envelope *global_env;
	MPI_Offset local_stripe_weigth;
	//MPI_Offset global_weight;

	std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *vect_env_weigt_1;
	std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *vect_env_weigt_2;
	std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *global_vect_env_weigt_1;
	std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *global_vect_env_weigt_2;

	std::vector<std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *> *envs_each_stripe_1;
	std::vector<std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *> *envs_each_stripe_2;

	std::vector<const geos::geom::Envelope *> *stripe_cells;

	void sort_and_collect_dataset1(std::vector<std::pair<const geos::geom::Envelope*, ulong>*>& global_vect,
        						   std::vector<std::pair<const geos::geom::Envelope*, ulong>*>& local_vect);

	void build_rtree_index(std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs_1);
	
	void sync_universe(const geos::geom::Envelope *universe);

	void shuffle_candidates(std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *envs_shuffled,
							std::vector<std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *> *envs_each_stripe, int tag);

	void collect_data(std::vector<std::pair<const geos::geom::Envelope*, ulong>*>& global_vect,
    					const std::vector<std::pair<const geos::geom::Envelope*, ulong>*>& local_vect);

	void parallel_sorting_by_regular_sampling(std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs);

	void clustering_envs(uint num_threads, std::vector<const geos::geom::Envelope *> *c_stripe_cells,
						 std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *vect_env_weigt,
						 std::vector<std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *> *envs_each_stripe);

	ullong adp_gpu_wrapper(MBR_Array *mbr_arr_1, MBR_Array *mbr_arr_2, MBR global_mbr, ullong target_num_cells, MBR_Array *partition_results);

};

#endif //ndef __ADP_GRID_H_INCLUDED__
