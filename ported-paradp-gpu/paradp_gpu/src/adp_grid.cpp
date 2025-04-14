#include <adp_grid.h>


void Adp_Grid ::build_rtree_index(std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs)
{
	for (std::list<std::pair<geos::geom::Envelope *, ulong>>::iterator itr = list_envs->begin(); itr != list_envs->end(); ++itr)
	{
		geos::geom::Envelope *temp_env = itr->first;
		//geos::geom::Geometry *temp_geom = *itr;
		index.insert(temp_env, temp_env);
	}

	return;
}

/* Need MIN(min_x min_y), MAX(max_x, max_y) accross all nodes*/
void Adp_Grid ::sync_universe(const geos::geom::Envelope *universe)
{
	double send_buf[4];
	double recv_buf[4];

	send_buf[0] = universe->getMinX();
	send_buf[1] = universe->getMinY();
	send_buf[2] = universe->getMaxX();
	send_buf[3] = universe->getMaxY();

	MPI_Allreduce(send_buf, recv_buf, 2, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(send_buf + 2, recv_buf + 2, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	global_env = new geos::geom::Envelope(recv_buf[0], recv_buf[2],
										  recv_buf[1], recv_buf[3]);

	return;
}

/* To get stripes */
void Adp_Grid ::parallel_sorting_by_regular_sampling(std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs)
{
	//1st step, sort locally
	//mergeSort(localMaxX, 0, v1Size-1);
	double *arr_x_max = NULL;
	ulong counter = 0;
	int num_world_nodes, my_world_rank, i;

	MPI_Comm_size(MPI_COMM_WORLD, &num_world_nodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_world_rank);

	arr_x_max = (double *)malloc(sizeof(double) * list_envs->size());

	for (std::list<std::pair<geos::geom::Envelope *, ulong> >::iterator itr = list_envs->begin(); itr != list_envs->end(); ++itr)
	{
		geos::geom::Envelope *temp_geom = itr->first;

		arr_x_max[counter] = temp_geom->getMaxX();
		//BugFixing: counter was always 0 and never get incremented!
		counter++;
	}

	Util_Quick_sort_double(arr_x_max, 0, list_envs->size() - 1);

	//2nd step, pick samples
	double *send_buf_samples = NULL;
	int sample_interval;

	send_buf_samples = (double *)malloc(num_world_nodes * sizeof(double));
	sample_interval = list_envs->size() / (num_world_nodes + 1);

	for (i = 0; i < num_world_nodes; ++i)
	{
		send_buf_samples[i] = arr_x_max[i * sample_interval];
	}

	//3rd step, rank 0 gathers all samples;
	double *recv_buf_samples = NULL;

	if (0 == my_world_rank)
	{
		recv_buf_samples = (double *)malloc(num_world_nodes * num_world_nodes * sizeof(double));
	}

	MPI_Gather(send_buf_samples, num_world_nodes, MPI_DOUBLE, recv_buf_samples, num_world_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Final step, rank 0 bcast all stripes
	double *stripes = NULL;
	stripes = (double *)malloc(num_world_nodes * sizeof(double));

	if (0 == my_world_rank)
	{
		Util_Quick_sort_double(recv_buf_samples, 0, num_world_nodes * num_world_nodes - 1);

		for (i = 0; i < num_world_nodes; ++i)
		{
			stripes[i] = recv_buf_samples[i * num_world_nodes];
		}
	}

	MPI_Bcast(stripes, num_world_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// All nodes generate same striples
	stripe_cells = new std::vector<const geos::geom::Envelope *>();
	for (i = 0; i < num_world_nodes; ++i)
	{
		double temp_min_y, temp_max_y;
		temp_min_y = global_env->getMinY();
		temp_max_y = global_env->getMaxY();

		if (0 == i)
		{
			//Stripe 0 uses the world envelope's minX as its minX
			stripe_cells->push_back(new const geos::geom::Envelope(global_env->getMinX(), stripes[i + 1], temp_min_y, temp_max_y));
		}
		else if (num_world_nodes - 1 == i)
		{
			//Last stripe uses the world envelope's maxX as its maxX
			stripe_cells->push_back(new const geos::geom::Envelope(stripes[i], global_env->getMaxX(), temp_min_y, temp_max_y));
		}
		else
		{
			stripe_cells->push_back(new const geos::geom::Envelope(stripes[i], stripes[i + 1], temp_min_y, temp_max_y));
		}
	}

	if (arr_x_max != NULL)
		free(arr_x_max);

	if (send_buf_samples != NULL)
		free(send_buf_samples);

	if (recv_buf_samples != NULL)
		free(recv_buf_samples);

	if (stripes != NULL)
		free(stripes);
}
/*
// The code I am writing from starts here. I am modifying the clustering_envs function and another function computeZOrder 
// for z sorting each stripe vectors
*/ 

// Function to compute Morton (Z-order) code for a given point
ullong computeZOrder(double x, double y) {
    ullong ix = static_cast<ullong>((x + 180.0) * 1e6); // Normalize longitude
    ullong iy = static_cast<ullong>((y + 90.0) * 1e6);  // Normalize latitude
    
    ullong z = 0;
    for (int i = 0; i < 32; i++) { // Interleave bits
        z |= ((ix >> i) & 1ULL) << (2 * i);
        z |= ((iy >> i) & 1ULL) << (2 * i + 1);
    }
    return z;
}

// // Updated clustering_envs function
// void Adp_Grid::clustering_envs(uint num_threads, 
//                                std::vector<const geos::geom::Envelope *> *c_stripe_cells,
//                                std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *vect_env_weigt,
//                                std::vector<std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *> *envs_each_stripe)
// {
//     uint i;
//     std::vector<std::thread> vect_threads;
//     std::mutex push_mutex;

//     for (i = 0; i < num_threads; ++i)
//     {
//         vect_threads.push_back(
//     		std::thread([i](std::vector<const geos::geom::Envelope *> *t_stripe_cells, 
//             		              uint t_num_threads, std::mutex *t_mutex,
//                 		          std::vector<std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *> *t_envs_each_stripe)
//     		{
//         		uint num_stripes = t_stripe_cells->size();

// 				for (uint j = 0; j < num_stripes; ++j)
// 				{
//     				// Create a copy of the original stripe for sorting
//     				std::vector<std::pair<const geos::geom::Envelope *, ulong> *> sorted_stripe = *(t_envs_each_stripe->at(j));

//     				// Sort the copied vector based on Z-order keys
//     				std::sort(sorted_stripe.begin(), sorted_stripe.end(),
//               				[](const auto &a, const auto &b) {
//                   				double ax = (a->first->getMinX() + a->first->getMaxX()) / 2.0;
//                   				double ay = (a->first->getMinY() + a->first->getMaxY()) / 2.0;
//                   				double bx = (b->first->getMinX() + b->first->getMaxX()) / 2.0;
//                   				double by = (b->first->getMinY() + b->first->getMaxY()) / 2.0;
//                   				return computeZOrder(ax, ay) < computeZOrder(bx, by);
//               				});

//     				// Lock before modifying shared memory
//     				t_mutex->lock();
//     				*(t_envs_each_stripe->at(j)) = std::move(sorted_stripe); // Assign back sorted stripe
//     				t_mutex->unlock();
// 				}


//     		}, c_stripe_cells, num_threads, &push_mutex, envs_each_stripe));

//     }

//     std::for_each(vect_threads.begin(), vect_threads.end(), [](std::thread &t) { t.join(); });
// }


void Adp_Grid ::clustering_envs(uint num_threads, std::vector<const geos::geom::Envelope *> *c_stripe_cells,
								std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *vect_env_weigt,
								std::vector<std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *> *envs_each_stripe)
{
	uint i;
	std::vector<std::thread> vect_threads;
	std::mutex push_mutex;

	for (i = 0; i < num_threads; ++i)
	{
		vect_threads.push_back(
			std::thread([i](std::vector<const geos::geom::Envelope *> *t_stripe_cells, uint t_num_threads, std::mutex *t_mutex,
							std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *t_vect_env_weigt,
							std::vector<std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *> *t_envs_each_stripe)
						{
							uint j, num_stripes = t_stripe_cells->size();
							uint num_envs_per_thread = t_vect_env_weigt->size() / t_num_threads;
							uint start = i * num_envs_per_thread;
							uint end = (i + 1) * num_envs_per_thread;

							if (t_num_threads - 1 == i)
							{
								end = t_vect_env_weigt->size();
							}

							std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *local_vect_env_weigt =
								new std::vector<std::pair<const geos::geom::Envelope *, ulong> *>(t_vect_env_weigt->begin() + start, t_vect_env_weigt->begin() + end - 1);

							std::vector<std::vector<std::pair<const geos::geom::Envelope *, ulong>>> *envs_in_stripes =
								new std::vector<std::vector<std::pair<const geos::geom::Envelope *, ulong>>>(num_stripes);

							geos::index::strtree::STRtree t_index;

							for (std::vector<std::pair<const geos::geom::Envelope *, ulong> *>::iterator itr = local_vect_env_weigt->begin();
								 itr != local_vect_env_weigt->end(); ++itr)
							{
								std::pair<const geos::geom::Envelope *, ulong> *temp_pair =
									new std::pair<const geos::geom::Envelope *, ulong>((*itr)->first, (*itr)->second);

								t_index.insert(temp_pair->first, temp_pair);
							}

							for (j = 0; j < num_stripes; ++j)
							{
								std::vector<void *> results;

								t_index.query(t_stripe_cells->at(j), results);

								for (std::vector<void *>::iterator itr = results.begin(); itr != results.end(); ++itr)
								{
									void *temp_pair_itr = *itr;
									std::pair<const geos::geom::Envelope *, ulong> *temp_pair =
										(std::pair<const geos::geom::Envelope *, ulong> *)temp_pair_itr;
									envs_in_stripes->at(j).push_back(*temp_pair);
								}
							}

							std::vector<void *> results;

							t_mutex->lock();
							for (j = 0; j < num_stripes; ++j)
							{
								for (std::vector<std::pair<const geos::geom::Envelope *, ulong>>::iterator itr = envs_in_stripes->at(j).begin();
									 itr != envs_in_stripes->at(j).end(); ++itr)
								{
									std::pair<const geos::geom::Envelope *, ulong> *temp_pair =
										new std::pair<const geos::geom::Envelope *, ulong>((*itr).first, (*itr).second);

									t_envs_each_stripe->at(j)->push_back(temp_pair);
								}
							}
							t_mutex->unlock();
						},
						c_stripe_cells, num_threads, &push_mutex, vect_env_weigt, envs_each_stripe));
	}

	// Looping every thread via for_each
	// The 3rd argument assigns a task
	// It tells the compiler we're using lambda ([])
	// The lambda function takes its argument as a reference to a thread, t
	// Then, joins one by one, and this works like barrier
	std::for_each(vect_threads.begin(), vect_threads.end(), [](std::thread &t)
				  { t.join(); });
}

void Adp_Grid ::shuffle_candidates(std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *envs_shuffled,
								   std::vector<std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *> *envs_each_stripe, int tag)
{
	int num_world_nodes, my_world_rank, i, j;

	MPI_Comm_size(MPI_COMM_WORLD, &num_world_nodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_world_rank);

	assert(envs_shuffled != NULL);
	assert(envs_each_stripe != NULL);

	/* Send envs of stripe i to proces i*/
	for (i = 0; i < num_world_nodes; ++i)
	{
		if (i == my_world_rank)
		{
			/* local stripe, no need send */
			for (std::vector<std::pair<const geos::geom::Envelope *, ulong> *>::iterator itr = envs_each_stripe->at(i)->begin();
				 itr != envs_each_stripe->at(i)->end(); ++itr)
			{
				std::pair<const geos::geom::Envelope *, ulong> *temp_pair =
					new std::pair<const geos::geom::Envelope *, ulong>((*itr)->first, (*itr)->second);

				envs_shuffled->push_back(temp_pair);
			}
		}
		else
		{
            spdlog::info("Send {0:d} to P{1:d}",  envs_each_stripe->at(i)->size(), i);
			int send_size = envs_each_stripe->at(i)->size();
			double *send_buf_envs;
			ulong *send_buf_weights;
			MPI_Request request;

			send_buf_envs = (double *)malloc(send_size * 4 * sizeof(double));
			send_buf_weights = (ulong *)malloc(send_size * sizeof(ulong));

			for (j = 0; j < send_size; ++j)
			{
				send_buf_envs[4 * j] = envs_each_stripe->at(i)->at(j)->first->getMinX();
				send_buf_envs[4 * j + 1] = envs_each_stripe->at(i)->at(j)->first->getMaxX();
				send_buf_envs[4 * j + 2] = envs_each_stripe->at(i)->at(j)->first->getMinY();
				send_buf_envs[4 * j + 3] = envs_each_stripe->at(i)->at(j)->first->getMaxY();
				send_buf_weights[j] = envs_each_stripe->at(i)->at(j)->second;
			}

			/* Tell process i how much data will be sent, 3 * num_world_nodes acts as flag */
			MPI_Isend(&send_size, 1, MPI_INT, i, 3 * num_world_nodes + tag, MPI_COMM_WORLD, &request);

			if (send_size > 0)
			{
				MPI_Isend(send_buf_envs, send_size * 4, MPI_DOUBLE, i, my_world_rank + 3 * num_world_nodes + tag, MPI_COMM_WORLD, &request);
				MPI_Isend(send_buf_weights, send_size, MPI_UNSIGNED_LONG, i, my_world_rank + 4 * num_world_nodes + tag, MPI_COMM_WORLD, &request);
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/* Recv envs of stripe my_world_rank */
	while (1)
	{
		MPI_Status status;
		MPI_Request request;
		int recv_size, recv_flag;
		std::chrono::duration<double> t_diff_temp;

		MPI_Irecv(&recv_size, 1, MPI_INT, MPI_ANY_SOURCE, 3 * num_world_nodes + tag, MPI_COMM_WORLD, &request);

		auto t_recv_begin = std::chrono::steady_clock::now();
		recv_flag = 0;

		while (!recv_flag)
		{
			MPI_Test(&request, &recv_flag, &status);
			auto t_recv_temp = std::chrono::steady_clock::now();

			t_diff_temp = t_recv_temp - t_recv_begin;

			/* wait upto 0.5 seconds */
			if (t_diff_temp.count() > 0.5)
				break;
		}

		/* Nothing to receive anymore, end this loop */
		if (!recv_flag)
			break;

        spdlog::info(" Recv {0:d} from P{1:d}",  recv_size, status.MPI_SOURCE);
		/* Begin to receive */
		if (recv_size > 0)
		{
			double *recv_buf_envs;
			ulong *recv_buf_weights;

			recv_buf_envs = (double *)malloc(recv_size * 4 * sizeof(double));
			recv_buf_weights = (ulong *)malloc(recv_size * sizeof(ulong));

			int sender = status.MPI_SOURCE;

			/* Block receive */
			MPI_Recv(recv_buf_envs, recv_size * 4, MPI_DOUBLE, sender, sender + 3 * num_world_nodes + tag, MPI_COMM_WORLD, &status);
			MPI_Recv(recv_buf_weights, recv_size, MPI_UNSIGNED_LONG , sender, sender + 4 * num_world_nodes + tag, MPI_COMM_WORLD, &status);

			for (i = 0; i < recv_size; ++i)
			{
				std::pair<const geos::geom::Envelope *, ulong> *temp_pair =
					new std::pair<const geos::geom::Envelope *, ulong>(new const geos::geom::Envelope(recv_buf_envs[4 * i], recv_buf_envs[4 * i + 1],
																									  recv_buf_envs[4 * i + 2], recv_buf_envs[4 * i + 3]),
																	   recv_buf_weights[i]);

				envs_shuffled->push_back(temp_pair);
			}
		}
	}

	/* To asure every process has finished*/
	MPI_Barrier(MPI_COMM_WORLD);
}

void Adp_Grid::collect_data(
    std::vector<std::pair<const geos::geom::Envelope*, unsigned long>*>& global_vect,
    const std::vector<std::pair<const geos::geom::Envelope*, unsigned long>*>& local_vect)
{
    int num_world_nodes, my_world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &num_world_nodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_world_rank);

    // Step 1: Pack local data into contiguous buffers.
    int local_count = static_cast<int>(local_vect.size());
    // Debug: Print local count on each process.
    fprintf(stderr, "Rank %d: local_count = %d\n", my_world_rank, local_count);

    double* local_envs = (double*) malloc(local_count * 4 * sizeof(double));
    unsigned long* local_weights = (unsigned long*) malloc(local_count * sizeof(unsigned long));
    if (!local_envs || !local_weights) {
        fprintf(stderr, "Memory allocation error on rank %d\n", my_world_rank);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    for (int i = 0; i < local_count; i++) {
        auto pair_ptr = local_vect[i];
        local_envs[4 * i]     = pair_ptr->first->getMinX();
        local_envs[4 * i + 1] = pair_ptr->first->getMaxX();
        local_envs[4 * i + 2] = pair_ptr->first->getMinY();
        local_envs[4 * i + 3] = pair_ptr->first->getMaxY();
        local_weights[i] = pair_ptr->second;
    }

    // Synchronize to ensure all processes are ready.
    MPI_Barrier(MPI_COMM_WORLD);

    // Step 2: Gather local counts from all processes.
    std::vector<int> all_counts(num_world_nodes, 0);
    MPI_Allgather(&local_count, 1, MPI_INT, all_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // Print the gathered counts on each process.
    if(my_world_rank==0) {
        fprintf(stderr, "Rank %d: all_counts: ", my_world_rank);
        for (int i = 0; i < num_world_nodes; i++) {
            fprintf(stderr, "%d ", all_counts[i]);
        }
        fprintf(stderr, "\n");
    }

    // Compute recvcounts and displacements for envelope data (4 doubles per envelope)
    int total_pairs = 0;
    std::vector<int> recvcounts(num_world_nodes, 0), displs(num_world_nodes, 0);
    for (int i = 0; i < num_world_nodes; i++) {
        recvcounts[i] = all_counts[i] * 4; // Each envelope contributes 4 doubles.
        displs[i] = total_pairs * 4;         // Displacement in number of doubles.
        total_pairs += all_counts[i];
    }
    // Debug: Print total pairs and computed recvcounts/displs.
    fprintf(stderr, "Rank %d: total_pairs = %d\n", my_world_rank, total_pairs);
    for (int i = 0; i < num_world_nodes; i++) {
        fprintf(stderr, "Rank %d: recvcounts[%d] = %d, displs[%d] = %d\n",
                my_world_rank, i, recvcounts[i], i, displs[i]);
    }

    // Step 3: Allocate global buffer for envelope coordinates.
    double* global_envs = (double*) malloc(total_pairs * 4 * sizeof(double));
    if (!global_envs) {
        fprintf(stderr, "Memory allocation error (global_envs) on rank %d\n", my_world_rank);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    // Gather all envelope data from all processes.
    MPI_Allgatherv(local_envs, local_count * 4, MPI_DOUBLE,
                   global_envs, recvcounts.data(), displs.data(), MPI_DOUBLE,
                   MPI_COMM_WORLD);

    // For weights, compute recvcounts and displacements (1 unsigned long per envelope)
    int total_weights = 0;
    std::vector<int> recvcounts_weights(num_world_nodes, 0), displs_weights(num_world_nodes, 0);
    for (int i = 0; i < num_world_nodes; i++) {
        recvcounts_weights[i] = all_counts[i];
        displs_weights[i] = total_weights;
        total_weights += all_counts[i];
    }
    unsigned long* global_weights = (unsigned long*) malloc(total_weights * sizeof(unsigned long));
    if (!global_weights) {
        fprintf(stderr, "Memory allocation error (global_weights) on rank %d\n", my_world_rank);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    MPI_Allgatherv(local_weights, local_count, MPI_UNSIGNED_LONG,
                   global_weights, recvcounts_weights.data(), displs_weights.data(), MPI_UNSIGNED_LONG,
                   MPI_COMM_WORLD);

    // Step 4: Reconstruct the global vector.
    global_vect.clear();
    for (int i = 0; i < total_pairs; i++) {
        double minX = global_envs[4 * i];
        double maxX = global_envs[4 * i + 1];
        double minY = global_envs[4 * i + 2];
        double maxY = global_envs[4 * i + 3];
        const geos::geom::Envelope* env = new geos::geom::Envelope(minX, maxX, minY, maxY);
        unsigned long weight = global_weights[i];
        global_vect.push_back(new std::pair<const geos::geom::Envelope*, unsigned long>(env, weight));
    }

    // Clean up temporary buffers.
    free(local_envs);
    free(local_weights);
    free(global_envs);
    free(global_weights);
}

void Adp_Grid::sort_and_collect_dataset1(
	std::vector<std::pair<const geos::geom::Envelope*, ulong>*>& global_vect,
        						   std::vector<std::pair<const geos::geom::Envelope*, ulong>*>& local_vect)
{
	std::vector<std::pair<const geos::geom::Envelope*, ulong>*> temp_global;
	collect_data(temp_global, local_vect);


	std::sort(temp_global.begin(), temp_global.end(),
			  [](const auto &a, const auto &b) {
				return a->first->getMaxX() < b->first->getMaxX();
			  });

	int my_rank, num_procs;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	int total = static_cast<int>(temp_global.size());

	int partition_size = total/ num_procs;
	int remainder = total % num_procs;

	int start_index = my_rank *partition_size + std::min(my_rank,remainder);
	int end_index = start_index + partition_size + (my_rank < remainder ? 1 : 0);

	global_vect.assign(temp_global.begin() + start_index, temp_global.begin()+end_index);
}


ullong Adp_Grid::adp_gpu_wrapper(MBR_Array *mbr_arr_1, MBR_Array *mbr_arr_2, MBR global_mbr, ullong global_num_cells, MBR_Array *partition_results)
{
	int num_world_nodes, my_world_rank;

	MPI_Comm_size(MPI_COMM_WORLD, &num_world_nodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_world_rank);
	if(my_world_rank == 0)
	{
		spdlog::info("CELL_WEIGTH_COEF = {0:d} MIN_NUM_CANDIDATES = {1:d} MAX_QUADTREE_DEPTH = {2:d}",
				 CELL_WEIGTH_COEF, MIN_NUM_CANDIDATES, MAX_QUADTREE_DEPTH);
	}

    auto t_begin = std::chrono::steady_clock::now();//get_current_time();
    double t_total, t_can, t_quad;

    Candidates candidates;
	ullong local_weight = 0ull;
    ullong num_candidates = 0ull;

	num_candidates = CUDA_Utils_candidate_generate_wrapper(mbr_arr_1, mbr_arr_2, DEFAULT_MBRS_PER_BLOCK, &candidates, &local_weight, my_world_rank);

	std::cout<< "Number of local candidates:" << num_candidates <<" in PID:" << my_world_rank << std::endl;
    auto t_can_end = std::chrono::steady_clock::now();//get_current_time();
#ifdef DEBUG
	spdlog::info("candidate_generate_wrapper DONE!");
#endif
	//Sync the global weight with other processes. 
	ullong local_num_cells = global_num_cells;
	ullong global_weight = 0;
	MPI_Allreduce(&local_weight, &global_weight, 1, MPI_UINT64_T/*MPI_OFFSET*/, MPI_SUM, MPI_COMM_WORLD);

	
	double d_global_weight = global_weight;
	double d_local_weigth = local_weight;
	double d_local_num_cells = global_num_cells / double(d_global_weight / d_local_weigth);

	local_num_cells = (uint)std::round(d_local_num_cells);
	//local_num_cells = global_num_cells / (global_weight / local_weight);

	if (local_num_cells == 0 && local_weight != 0)
		local_num_cells = 1;

#ifdef DEBUG
	spdlog::info("local_weight = {0:d} global_weight = {1:d} global_num_cells = {2:d} local_num_cells = {3:d}",
				 local_weight, global_weight, global_num_cells, local_num_cells);
#endif

	CUDA_Utils_quadtree_partition_wrapper(&candidates, num_candidates, global_mbr, local_weight, 
            local_num_cells, partition_results);
    
    auto t_quad_end = std::chrono::steady_clock::now();//get_current_time();
	
	std::chrono::duration<double> t_diff_temp;
	t_diff_temp = t_quad_end - t_begin;
    t_total = t_diff_temp.count();
	t_diff_temp = t_can_end - t_begin;
    t_can = t_diff_temp.count();
	t_diff_temp = t_quad_end - t_can_end;
    t_quad = t_diff_temp.count();

	#ifdef DEBUG
	spdlog::info("adp_gpu_wrapper total time {0:03.6f}, candidate generate {1:03.6f}, quadtree {2:03.6f}", 
        t_total, t_can, t_quad);
		//printf("cuda_utils_adp_wrapper total time %f, candidate generate %f, quadtree %f.\n", 
        //t_total, t_can, t_quad);

	spdlog::info("num_candidates = {0:d} partition_results size = {1:d}", num_candidates, partition_results->size);
	#endif

    return num_candidates;
}
