#include <global.h>

//mpirun -np 4 ./prog 8192 ./data/lakes_data_env ./data/sports_data_env 

int main(int argc, char **argv)
{
    //MPI initial    
    MPI_Init(&argc, &argv);

	int num_world_nodes, my_world_rank;

	MPI_Comm_size(MPI_COMM_WORLD, &num_world_nodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_world_rank);


#ifdef DEBUG
    spdlog::set_pattern("P" + std::to_string(my_world_rank) + " " + "[%H:%M:%S.%e] %v"); /*Format: h:m:s + content*/
#ifdef MEMORY_ENHANCMENT
    if(0 == my_world_rank)
        std::cout << "NOTE: MEMORY_ENHANCMENT is ENABLED!" << std::endl;
#endif
    spdlog::info("Enter function {}", "main");

    double t_total, t_parse, t_steal, t_join;
    t_total = t_steal = t_join = 0.0;
    std::chrono::duration<double> t_diff_temp;

    auto t_begin = std::chrono::steady_clock::now();
#endif

    const int num_partitions = std::stoi(argv[1]);
    const std::string file_path_1 = argv[2];
    const std::string file_path_2 = argv[3];

    const uint num_threads = 36;

    std::vector<std::string> *vect_strs_1 = new std::vector<std::string>();
    std::vector<std::string> *vect_strs_2 = new std::vector<std::string>();

    Util_Read_geoms_mpi(file_path_1, vect_strs_1, num_threads, MPI_COMM_WORLD);
    Util_Read_geoms_mpi(file_path_2, vect_strs_2, num_threads, MPI_COMM_WORLD);

    //std::list<geos::geom::Geometry *> *list_geoms_1 = new std::list<geos::geom::Geometry *>;
    //std::list<geos::geom::Geometry *> *list_geoms_2 = new std::list<geos::geom::Geometry *>;
    std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs_1 = new std::list<std::pair<geos::geom::Envelope *, ulong> >();
	std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs_2 =	new std::list<std::pair<geos::geom::Envelope *, ulong> >();

    //Util_Parse_string_parallel(vect_strs_1, list_geoms_1, num_threads);
    //Util_Parse_string_parallel(vect_strs_2, list_geoms_2, num_threads);
    if(0 == my_world_rank)
    {
        //std::cout << "P" << my_world_rank << " vect_strs_1 global env (1st line): " << vect_strs_1->front() << std::endl;
        //std::cout << "P" << my_world_rank << " vect_strs_2 global env (1st line): " << vect_strs_2->front() << std::endl;
        vect_strs_1->erase(vect_strs_1->begin());
        vect_strs_2->erase(vect_strs_2->begin());
        //std::cout << "P" << my_world_rank << " vect_strs_1[0]: " << vect_strs_1->front() << std::endl;
    }

// #ifdef DEBUG
//     spdlog::info("vect_strs_1 size: {0:d}, vect_strs_2 size: {1:d}", vect_strs_1->size(), vect_strs_2->size());
//     spdlog::info("vect_strs_1[0]:    {0:s}, vect_strs_2[0]:    {1:s}", vect_strs_1->front(), vect_strs_2->front());
//     spdlog::info("vect_strs_1[last]: {0:s}, vect_strs_2[last]: {1:s}", vect_strs_1->back(), vect_strs_2->back());
// #endif

    
    Util_Parse_env_string_parallel(vect_strs_1, list_envs_1, num_threads);
    Util_Parse_env_string_parallel(vect_strs_2, list_envs_2, num_threads);

#ifdef DEBUG
    //spdlog::info("DEBUG1");
#endif
#ifdef MEMORY_ENHANCMENT
    //Enhancement: Free up the vect_strs_1 & vect_strs_2 since it never be used later!
    if (vect_strs_1 != NULL)
        delete vect_strs_1;
    if (vect_strs_2 != NULL)
        delete vect_strs_2;
#endif

#ifdef DEBUG
    auto t_parse_end = std::chrono::steady_clock::now();

    t_diff_temp = t_parse_end - t_begin;

    t_parse = t_diff_temp.count();

    spdlog::info("List R size {0:d}, S size {1:d}, parsing took {2:03.3f} seconds", list_envs_1->size(), list_envs_2->size(), t_parse);
#endif

    const geos::geom::Envelope *universe = Util_Get_global_env(list_envs_1, list_envs_2);
    //const geos::geom::Envelope * universe = new geos::geom::Envelope(-180.0,180.0,-90.0,90.0);
    std::cout << "P" << my_world_rank << " universe: " << universe->toString() << std::endl;

    Adp_Grid *grid = new Adp_Grid(num_partitions, num_threads, universe, list_envs_1, list_envs_2);

    //std::list<std::pair<const geos::geom::Envelope *, ulong> *> * grid_result = grid->get_grid_with_weight();

    //std::cout << "P" << my_world_rank << " grid_result size: " << grid_result->size() << std::endl;

#ifdef DEBUG
    auto t_end = std::chrono::steady_clock::now();

    t_diff_temp = t_end - t_parse_end;
    double parAdp_gpu_time = t_diff_temp.count();
    spdlog::info("Total parAdp_gpu_time time used: {0:03.3f}", parAdp_gpu_time);

    t_diff_temp = t_end - t_begin;
    t_total = t_diff_temp.count();

    spdlog::info("Leave function {0:s}, total time used: {1:03.3f}", "main", t_total);
#endif

    MPI_Finalize();
    return 0;
}
