#include <global.h>

void parse_list_envs_to_mbr_arr(std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs, MBR_Array &mbr_arr);
void wirte_files(std::string file_path, std::string output_file_path, MBR_Array &grids);
void test(int argc, char **argv);
void test1(int argc, char **argv);

int main(int argc, char **argv)
{
	spdlog::set_pattern("[%H:%M:%S.%e] %v");

#ifdef DEBUG
	spdlog::set_level(spdlog::level::debug); // Set global log level to debug
#else
	spdlog::set_level(spdlog::level::info); // Set global log level to info
#endif

	test(argc, argv);// Single GPU "ADP-GPU ??"
	//test1(argc, argv);// Multiple GPUs "ADP-mGPU ??"

	return 0;
}


void test1(int argc, char **argv)
{
	const int num_partitions = atoi(argv[1]);
	spdlog::info("Env File 1 {}, file 2 {}, target parts {}", argv[2], argv[3], num_partitions);
	spdlog::info("Geoms File 1 {}, file 2 {}", argv[4], argv[5]);

    const std::string geoms_file_1 = argv[4];
    const std::string geoms_file_2 = argv[5];

	time_t t_start, t_end_parsing, t_end_build_grid, t_end_index, t_end;
	t_start = time(0);

	std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs_1 =
		new std::list<std::pair<geos::geom::Envelope *, ulong> >();
	std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs_2 =
		new std::list<std::pair<geos::geom::Envelope *, ulong> >();

	//geos::geom::Envelope *universe_1 = Util_Read_envs(argv[2], list_envs_1);
	//geos::geom::Envelope *universe_2 = Util_Read_envs(argv[3], list_envs_2);

	geos::geom::Envelope *universe_1 = Util_Read_envs_parallel(argv[2], list_envs_1);
	geos::geom::Envelope *universe_2 = Util_Read_envs_parallel(argv[3], list_envs_2);

	spdlog::debug("{}", universe_1->toString());
	spdlog::debug("{}", universe_2->toString());

	MBR_Array mbr_arr_1;
	parse_list_envs_to_mbr_arr(list_envs_1, mbr_arr_1);

	spdlog::debug("MBR 1: size {}, {} {}", mbr_arr_1.size, mbr_arr_1.mbrs[0].min_x, mbr_arr_1.mbrs[0].min_y);

	MBR_Array mbr_arr_2;
	parse_list_envs_to_mbr_arr(list_envs_2, mbr_arr_2);

	spdlog::debug("MBR 2: size {}, {} {}", mbr_arr_2.size, mbr_arr_2.mbrs[0].min_x, mbr_arr_2.mbrs[0].min_y);

	GPU_Utils_wakeup_wrapper(WAKEUP_CYCLE);

    //int num_gpus = gpu_count();
	int num_gpus = atoi(argv[6]);

	spdlog::debug("Begin CUDA ADP with GPUS {}", num_gpus);
    MBR temp_mbr = {-180.0, -90.0, 180.0, 90.0};
	MBR_Array partition_results[num_gpus];
    std::vector<std::thread> vect_threads;

	ullong num_candidates = 0ULL;
   for (int i = 0; i < num_gpus; ++i)
    {
        vect_threads.push_back(std::thread(CUDA_Utils_gpus_wrapper, &mbr_arr_1, &mbr_arr_2, temp_mbr, num_partitions, &partition_results[i], i, num_gpus));
    }


    for (int i = 0; i < vect_threads.size(); ++i)
    {
        vect_threads[i].join();
    }

	spdlog::debug("End CUDA ADP");

	//wirte_files(geoms_file_1, "A", partition_results);
	//wirte_files(geoms_file_2, "B", partition_results);

	spdlog::debug("Write finished. Exit.");
}


void test(int argc, char **argv)
{
	const int num_partitions = atoi(argv[1]);
	spdlog::info("Env File 1 {}, file 2 {}, target parts {}", argv[2], argv[3], num_partitions);
	//spdlog::info("Geoms File 1 {}, file 2 {}", argv[4], argv[5]);

    //const std::string geoms_file_1 = argv[4];
    //const std::string vim  = argv[5];

	time_t t_start, t_end_parsing, t_end_build_grid, t_end_index, t_end;
	t_start = time(0);

	std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs_1 =
		new std::list<std::pair<geos::geom::Envelope *, ulong> >();
	std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs_2 =
		new std::list<std::pair<geos::geom::Envelope *, ulong> >();

	//geos::geom::Envelope *universe_1 = Util_Read_envs(argv[2], list_envs_1);
	//geos::geom::Envelope *universe_2 = Util_Read_envs(argv[3], list_envs_2);

	geos::geom::Envelope *universe_1 = Util_Read_envs_parallel(argv[2], list_envs_1);
	geos::geom::Envelope *universe_2 = Util_Read_envs_parallel(argv[3], list_envs_2);

	spdlog::debug("{}", universe_1->toString());
	spdlog::debug("{}", universe_2->toString());

	MBR_Array mbr_arr_1;
	parse_list_envs_to_mbr_arr(list_envs_1, mbr_arr_1);

	spdlog::debug("MBR 1: size {}, {} {}", mbr_arr_1.size, mbr_arr_1.mbrs[0].min_x, mbr_arr_1.mbrs[0].min_y);

	MBR_Array mbr_arr_2;
	parse_list_envs_to_mbr_arr(list_envs_2, mbr_arr_2);

	spdlog::debug("MBR 2: size {}, {} {}", mbr_arr_2.size, mbr_arr_2.mbrs[0].min_x, mbr_arr_2.mbrs[0].min_y);

	GPU_Utils_wakeup_wrapper(WAKEUP_CYCLE);

	ullong num_candidates = 0ULL;

	spdlog::debug("Begin CUDA ADP");

	MBR temp_mbr = {-180.0, -90.0, 180.0, 90.0};
	MBR_Array partition_results;
	num_candidates = CUDA_Utils_adp_wrapper(&mbr_arr_1, &mbr_arr_2, temp_mbr, num_partitions, &partition_results);

	spdlog::debug("num_candidates = {0:d} partition_results size = {1:d}", num_candidates, partition_results.size);
	spdlog::debug("End CUDA ADP");

	//wirte_files(geoms_file_1, "A", partition_results);
	//wirte_files(geoms_file_2, "B", partition_results);

	spdlog::debug("Write finished. Exit.");
}

void parse_list_envs_to_mbr_arr(std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs, MBR_Array &mbr_arr)
{
	mbr_arr.size = list_envs->size();
	mbr_arr.mbrs = (MBR*)malloc(mbr_arr.size * sizeof(MBR));
	mbr_arr.weights = (ullong*)malloc(mbr_arr.size * sizeof(ullong));

	ullong counter = 0ULL;

	for(auto &itr : *list_envs)
	{
		mbr_arr.mbrs[counter].min_x = itr.first->getMinX();
		mbr_arr.mbrs[counter].min_y = itr.first->getMinY();
		mbr_arr.mbrs[counter].max_x = itr.first->getMaxX();
		mbr_arr.mbrs[counter].max_y = itr.first->getMaxY();
		mbr_arr.weights[counter] = itr.second;
		++counter;
	}
}

void parse_mbr_arr_to_list_envs(MBR_Array &mbr_arr, std::list<geos::geom::Envelope *> *list_envs, bool is_print)
{
	assert(list_envs);
	for (ullong counter = 0; counter < mbr_arr.size; ++counter)
	{
		list_envs->push_back(new geos::geom::Envelope(mbr_arr.mbrs[counter].min_x, mbr_arr.mbrs[counter].max_x, mbr_arr.mbrs[counter].min_y, mbr_arr.mbrs[counter].max_y));

		if (is_print)
		{
			double min_x = mbr_arr.mbrs[counter].min_x;
			double min_y = mbr_arr.mbrs[counter].min_y;
			double max_x = mbr_arr.mbrs[counter].max_x;
			double max_y = mbr_arr.mbrs[counter].max_y; 

			std::cout<<"POLYGON ((" + std::to_string(min_x) + " " + std::to_string(min_y) + "," + std::to_string(min_x) + " " + std::to_string(max_y) 
				+ "," + std::to_string(max_x) + " " + std::to_string(max_y) + "," + std::to_string(max_x) + " " + std::to_string(min_y) + "," 
				+ std::to_string(min_x) + " " + std::to_string(min_y) + +"))";
			std::cout<<std::endl;
		}
	}
}

void wirte_files(std::string file_name, std::string output_file_path, MBR_Array &grids)
{
	std::list<geos::geom::Envelope *> *list_envs = new std::list<geos::geom::Envelope *>();
	if (output_file_path.compare("A") == 0 )
		parse_mbr_arr_to_list_envs(grids, list_envs, true);
	else
		parse_mbr_arr_to_list_envs(grids, list_envs, false);

	std::list<geos::geom::Geometry*> *list_geoms = new std::list<geos::geom::Geometry*>();
	Util_Read_geoms(file_name, list_geoms);

	geos::index::strtree::STRtree index;

	for (std::list<geos::geom::Geometry *>::iterator itr = list_geoms->begin(); itr != list_geoms->end(); ++itr)
	{
		geos::geom::Geometry *geom_ptr = *itr;

		index.insert(geom_ptr->getEnvelopeInternal(), geom_ptr);
	}

	//Make dirs
	char cwd_path[100];

	if (NULL == getcwd(cwd_path, sizeof(cwd_path)))
	{
		std::cout << "Fail to make directory.\n"
				  << std::endl;
		exit(1);
	}

	std::string cwd_path_str(cwd_path);
	std::string write_file_path = cwd_path_str + "/" + output_file_path;

	if (-1 == access(write_file_path.c_str(), F_OK))
	{
			mkdir(write_file_path.c_str(),S_IRWXU);
	}

	int name_counter = 0;
	for (std::list<geos::geom::Envelope *> ::iterator itr = list_envs->begin(); itr != list_envs->end(); ++itr)
	{
		geos::geom::Envelope* tmpEnv = *itr;
	
		std::string temp_file_name = write_file_path + "/" + std::to_string(name_counter);
	
		Util_Write_geoms_to_WKT(temp_file_name, tmpEnv, &index);
		
		name_counter++;
	}

}






