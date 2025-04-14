#include <grid_utils.h>

void Util_Read_geoms(std::string file_path, std::list<geos::geom::Geometry *> *list_geoms)
{
	std::ifstream file(file_path.c_str());
	std::list<std::string> *list_strs = new std::list<std::string>;
	std::string temp_str;

	while (std::getline(file, temp_str))
	{
		//omit empty strings and invalid strings
		if (temp_str.size() > 5)
			list_strs->push_back(temp_str);
	}

	file.close();

	if (VERBOSE)
		std::cout << "File size::" << list_strs->size() << std::endl;

	geos::io::WKTReader wktreader;
	geos::geom::Geometry *temp_geom = NULL;

	for (std::list<std::string>::iterator itr = list_strs->begin(); itr != list_strs->end(); ++itr)
	{
		temp_str = *itr;
		temp_geom = NULL;

		try
		{
			temp_geom = (wktreader.read(temp_str)).release();
		}
		catch (std::exception &e)
		{
			//throw;
		}

		if (temp_geom != NULL && temp_geom->isValid())
		{
			list_geoms->push_back(temp_geom);
		}
	}

	delete list_strs;
}

void tokenize(std::vector<std::string> &tokens, std::string str, std::string del = " ")
{
	int start = 0;
	int end = str.find(del);
	while (end != -1)
	{
		tokens.push_back(str.substr(start, end - start));
		start = end + del.size();
		end = str.find(del, start);
	}

	tokens.push_back(str.substr(start, end - start));
}

void thread_read_geoms_from_strs(std::vector<std::string> *vect_strs, std::list<geos::geom::Geometry *> *l_geoms,
								 geos::io::WKTReader *reader, std::mutex *push_mutex)
{
	std::string temp_str;
	geos::geom::Geometry *temp_geom = NULL;
	std::list<geos::geom::Geometry *> *thread_l_geoms = new std::list<geos::geom::Geometry *>();
	uint i;

	for (std::vector<std::string>::iterator itr = vect_strs->begin(); itr != vect_strs->end(); ++itr)
	{
		temp_str = *itr;
		temp_geom = NULL;

		try
		{
			temp_geom = (reader->read(temp_str)).release();
		}
		catch (std::exception &e)
		{
			//throw;
		}

		if (temp_geom != NULL)
		{ // && temp_geom->isValid()){//Disable validation check increase performance
			if (temp_geom->getGeometryTypeId() == geos::geom::GeometryTypeId::GEOS_MULTIPOLYGON ||
				temp_geom->getGeometryTypeId() == geos::geom::GeometryTypeId::GEOS_MULTILINESTRING ||
				temp_geom->getGeometryTypeId() == geos::geom::GeometryTypeId::GEOS_GEOMETRYCOLLECTION)
			{
				for (i = 0; i < temp_geom->getNumGeometries(); ++i)
				{
					thread_l_geoms->push_back(const_cast<geos::geom::Geometry *>(temp_geom->getGeometryN(i)));
				}
			}
			else
			{
				thread_l_geoms->push_back(temp_geom);
			}
		}
	}

	push_mutex->lock();

	while (!thread_l_geoms->empty())
	{
		l_geoms->push_back(thread_l_geoms->back());
		thread_l_geoms->pop_back();
	}

	push_mutex->unlock();

	//Do not free temp_geom, it's used in the future.
	//To avoid memory leak, we may free it after using.
}

void parallel_read_env(std::vector<std::string> *vect_strs, 
						std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs, std::mutex *push_mutex)
{
	std::list<std::pair<geos::geom::Envelope *, ulong> > *temp_list_envs = 
			new std::list<std::pair<geos::geom::Envelope *, ulong> >();

	std::vector<std::string> tokens;		
	for (auto &itr : *vect_strs)
	{	

		std::string temp_str = itr;

		//omit empty strings and invalid strings
		if (temp_str.size() > 5)
		{
			tokens.clear();
			tokenize(tokens, temp_str);

			if(tokens.size() < 6)
			{
				try
				{
					geos::geom::Envelope *temp_env = new geos::geom::Envelope(std::stod(tokens[0]), std::stod(tokens[1]),
																		std::stod(tokens[2]), std::stod(tokens[3]));

					ulong temp_weight = std::stoul(tokens[4]);
					if(!temp_env->isNull()){ 
						temp_list_envs->push_back(std::make_pair(temp_env, temp_weight));
					}
				}
				catch(const std::exception& e)
				{
					//std::cerr << e.what() << '\n'; 
				}
			}
		}
	}

	push_mutex->lock();
	for (auto &itr : *temp_list_envs)
		list_envs->push_back(itr);
	push_mutex->unlock();
}

void Util_Parse_string_parallel(std::vector<std::string> *vect_strs, std::list<geos::geom::Geometry *> *list_geoms, uint num_threads)
{
	uint i, vect_size;
	geos::io::WKTReader wkt_reader[num_threads];
	std::vector<std::vector<std::string> *> vect_vect_strs(num_threads, NULL);
	std::vector<std::string>::iterator sub_itr_begin, sub_itr_end;
	std::vector<std::thread> *local_thread_vect = new std::vector<std::thread>();

	vect_size = vect_strs->size();
	vect_size = vect_size / num_threads;

	std::mutex push_mutex;

	for (i = 0; i < num_threads; ++i)
	{
		sub_itr_begin = vect_strs->begin() + i * vect_size;

		if (i == num_threads - 1)
			sub_itr_end = vect_strs->end();
		else
			sub_itr_end = vect_strs->begin() + (i + 1) * vect_size;

		vect_vect_strs[i] = new std::vector<std::string>(sub_itr_begin, sub_itr_end);

		wkt_reader[i] = geos::io::WKTReader();

		local_thread_vect->push_back(std::thread(thread_read_geoms_from_strs,
												 vect_vect_strs[i], list_geoms, &wkt_reader[i], &push_mutex));
	}

	for (std::vector<std::thread>::iterator itr = local_thread_vect->begin(); itr != local_thread_vect->end(); ++itr)
	{
		(*itr).join();
	}

	delete local_thread_vect;

}

void Util_Parse_env_string_parallel(std::vector<std::string> *vect_strs, std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs, uint num_threads)
{
	std::mutex push_mutex;
	std::vector<std::thread> readers;
	uint vect_size = vect_strs->size()/num_threads;
	std::vector<std::string>::iterator sub_itr_begin, sub_itr_end;
	std::vector<std::vector<std::string> *> vect_vect_strs(num_threads, NULL);

	for (uint i = 0; i < num_threads; ++i)
	{
		sub_itr_begin = vect_strs->begin() + i * vect_size;

		if (i == num_threads - 1)
                sub_itr_end = vect_strs->end();
            else
                sub_itr_end = vect_strs->begin() + (i + 1) * vect_size;

		vect_vect_strs[i] = new std::vector<std::string>(sub_itr_begin, sub_itr_end);



		readers.push_back(std::thread(parallel_read_env, vect_vect_strs[i], list_envs, &push_mutex));
	}

	std::for_each(readers.begin(), readers.end(), [](std::thread &t) 
    {
        t.join();
    });

	//delete vect_strs;
}

void Util_Read_geoms_parallel(std::string file_path, std::list<geos::geom::Geometry *> *list_geoms, uint num_threads)
{
	geos::io::WKTReader wkt_reader[num_threads];
	std::vector<std::string> *vect_strs = new std::vector<std::string>;

	std::ifstream file(file_path.c_str());
	std::string temp_str;

	while (std::getline(file, temp_str))
	{
		//omit empty strings and invalid strings
		if (temp_str.size() > 5)
			vect_strs->push_back(temp_str);
	}

	file.close();

	if (vect_strs->empty())
		return;

	if (vect_strs->size() < 200) //Too tiny for parallel parsing
	{
		Util_Read_geoms(file_path, list_geoms);
		return;
	}

	Util_Parse_string_parallel(vect_strs, list_geoms, num_threads);

	vect_strs->clear();
	vect_strs->shrink_to_fit();
}

void Util_Write_geoms_to_WKT(std::string file_name, const geos::geom::Envelope *env,
							 geos::index::strtree::STRtree *index)
{
	//Finding geometries
	std::vector<void *> results;

	index->query(env, results);

	//do not produce empty files
	if (results.empty())
	{
		return;
	}

	std::ofstream temp_file;
	temp_file.open(file_name);

	//Writing cell bounding box
	std::string env_str = std::to_string(env->getMinX()) + " " +
						  std::to_string(env->getMaxX()) + " " +
						  std::to_string(env->getMinY()) + " " +
						  std::to_string(env->getMaxY()) + "\n";

	temp_file << env_str;

	for (std::vector<void *>::iterator void_itr = results.begin(); void_itr != results.end(); ++void_itr)
	{
		void *void_geom_ptr = *void_itr;

		geos::geom::Geometry *query_geom = (geos::geom::Geometry *)void_geom_ptr;
		temp_file << query_geom->toString() + "\n";
	}

	temp_file.close();
	return;
}

geos::geom::Geometry *Util_Covert_env_to_geom(const geos::geom::Envelope *env)
{
	if (env == NULL)
		return NULL;

	double min_x = env->getMinX();
	double max_x = env->getMaxX();
	double min_y = env->getMinY();
	double max_y = env->getMaxY();

	std::string temp_str = "POLYGON ((" + std::to_string(min_x) + " " + std::to_string(min_y) + "," +
						   std::to_string(min_x) + " " + std::to_string(max_y) + "," + std::to_string(max_x) + " " +
						   std::to_string(max_y) + "," + std::to_string(max_x) + " " + std::to_string(min_y) + "," +
						   std::to_string(min_x) + " " + std::to_string(min_y) + +"))";

	geos::io::WKTReader wktreader;

	geos::geom::Geometry *temp_geom = NULL;

	try
	{
		temp_geom = (wktreader.read(temp_str)).release();
	}
	catch (std::exception &e)
	{
		//throw;
	}

	return temp_geom;
}

void Util_Parse_vector_envs_to_mbr_arr(std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *vector_envs, MBR_Array &mbr_arr)
{
	mbr_arr.size = vector_envs->size();
	mbr_arr.mbrs = (MBR*)malloc(mbr_arr.size * sizeof(MBR));
	mbr_arr.weights = (ullong*)malloc(mbr_arr.size * sizeof(ullong));

	ullong counter = 0ULL;

	for(auto &itr : *vector_envs)
	{
		mbr_arr.mbrs[counter].min_x = itr->first->getMinX();
		mbr_arr.mbrs[counter].min_y = itr->first->getMinY();
		mbr_arr.mbrs[counter].max_x = itr->first->getMaxX();
		mbr_arr.mbrs[counter].max_y = itr->first->getMaxY();
		mbr_arr.weights[counter] = itr->second;
		++counter;
	}
}

void Util_Write_grid_to_WKT(std::string file_name, std::list<const geos::geom::Envelope *> *list_envs)
{
	std::ofstream temp_file;
	temp_file.open(file_name);

	geos::geom::Geometry *temp_geom = NULL;
	for (std::list<const geos::geom::Envelope *>::iterator itr = list_envs->begin();
		 itr != list_envs->end(); itr++)
	{
		const geos::geom::Envelope *temp_env = *itr;

		temp_geom = Util_Covert_env_to_geom(temp_env);

		if (temp_geom != NULL)
			temp_file << temp_geom->toString() + "\n";
	}

	temp_file.close();
}

void Util_Write_geoms_to_array(std::string file_name, const geos::geom::Envelope *env,
							   geos::index::strtree::STRtree *index)
{
	std::vector<void *> results;

	index->query(env, results);

	//not produce empty files
	if (results.empty())
	{
		return;
	}

	std::ofstream tmpfile;

	tmpfile.open(file_name);

	std::string x_array_str;
	std::string y_array_str;
	std::string index_array_str;
	std::string envs_array_str;
	int count = 0;

	for (std::vector<void *>::iterator vdItr = results.begin(); vdItr != results.end(); ++vdItr)
	{
		void *void_geom_ptr = *vdItr;
		geos::geom::Geometry *temp_geom = (geos::geom::Geometry *)void_geom_ptr;

		envs_array_str.append(std::to_string(temp_geom->getEnvelopeInternal()->getMinX()) + " ");
		envs_array_str.append(std::to_string(temp_geom->getEnvelopeInternal()->getMaxX()) + " ");
		envs_array_str.append(std::to_string(temp_geom->getEnvelopeInternal()->getMinY()) + " ");
		envs_array_str.append(std::to_string(temp_geom->getEnvelopeInternal()->getMaxY()) + " ");

		geos::geom::CoordinateSequence *coord_seq = (temp_geom->getCoordinates()).release();

		int endCount = coord_seq->getSize();
		for (int i = 0; i < endCount; i++)
		{
			x_array_str.append(std::to_string(coord_seq->getAt(i).x));

			y_array_str.append(std::to_string(coord_seq->getAt(i).y));

			x_array_str.append(" ");
			y_array_str.append(" ");
		}

		index_array_str.append(std::to_string(count));
		index_array_str.append(" ");
		count += endCount;
	}

	envs_array_str.erase(envs_array_str.end() - 1);
	x_array_str.erase(x_array_str.end() - 1);
	y_array_str.erase(y_array_str.end() - 1);

	x_array_str.append("\n");
	y_array_str.append("\n");
	index_array_str.append(std::to_string(count));
	index_array_str.append("\n");
	envs_array_str.append("\n");

	std::string envStr = std::to_string(env->getMinX()) + " " +
						 std::to_string(env->getMaxX()) + " " +
						 std::to_string(env->getMinY()) + " " +
						 std::to_string(env->getMaxY());

	tmpfile << x_array_str;
	tmpfile << y_array_str;
	tmpfile << index_array_str;
	tmpfile << envs_array_str;
	tmpfile << envStr;

	tmpfile.close();

	return;
}

const geos::geom::Envelope *Util_Get_global_env(std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs_1,
												std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs_2)
{
	double min_x = 180.0;
	double max_x = -180.0;
	double min_y = 90.0;
	double max_y = -90.0;

	for (std::list<std::pair<geos::geom::Envelope *, ulong> >::iterator itr = list_envs_1->begin(); itr != list_envs_1->end(); ++itr)
	{
		geos::geom::Envelope *temp_env = itr->first;

		if (temp_env->getMinX() < min_x)
			min_x = temp_env->getMinX();
		if (temp_env->getMaxX() > max_x)
			max_x = temp_env->getMaxX();
		if (temp_env->getMinY() < min_y)
			min_y = temp_env->getMinY();
		if (temp_env->getMaxY() > max_y)
			max_y = temp_env->getMaxY();
	}

	for (std::list<std::pair<geos::geom::Envelope *, ulong> >::iterator itr = list_envs_2->begin(); itr != list_envs_2->end(); ++itr)
	{
		geos::geom::Envelope *temp_env = itr->first;

		if (temp_env->getMinX() < min_x)
			min_x = temp_env->getMinX();
		if (temp_env->getMaxX() > max_x)
			max_x = temp_env->getMaxX();
		if (temp_env->getMinY() < min_y)
			min_y = temp_env->getMinY();
		if (temp_env->getMaxY() > max_y)
			max_y = temp_env->getMaxY();
	}

	const geos::geom::Envelope *env = new geos::geom::Envelope(min_x, max_x, min_y, max_y);

	return env;
}

/*const geos::geom::Envelope *Util_Get_global_env(std::list<geos::geom::Geometry *> *layer_1_geoms,
												std::list<geos::geom::Geometry *> *layer_2_geoms)
{
	double min_x = 180.0;
	double max_x = -180.0;
	double min_y = 90.0;
	double max_y = -90.0;

	for (std::list<geos::geom::Geometry *>::iterator itr = layer_1_geoms->begin(); itr != layer_1_geoms->end(); ++itr)
	{
		geos::geom::Geometry *temp_geom = *itr;

		if (temp_geom == NULL || !temp_geom->isValid())
			continue;

		const geos::geom::Envelope *temp_env = temp_geom->getEnvelopeInternal();

		if (temp_env->getMinX() < min_x)
			min_x = temp_env->getMinX();
		if (temp_env->getMaxX() > max_x)
			max_x = temp_env->getMaxX();
		if (temp_env->getMinY() < min_y)
			min_y = temp_env->getMinY();
		if (temp_env->getMaxY() > max_y)
			max_y = temp_env->getMaxY();
	}

	for (std::list<geos::geom::Geometry *>::iterator itr = layer_2_geoms->begin();
		 itr != layer_2_geoms->end(); itr++)
	{
		geos::geom::Geometry *temp_geom = *itr;

		if (temp_geom == NULL || !temp_geom->isValid())
			continue;

		const geos::geom::Envelope *temp_env = temp_geom->getEnvelopeInternal();

		if (temp_env->getMinX() < min_x)
			min_x = temp_env->getMinX();
		if (temp_env->getMaxX() > max_x)
			max_x = temp_env->getMaxX();
		if (temp_env->getMinY() < min_y)
			min_y = temp_env->getMinY();
		if (temp_env->getMaxY() > max_y)
			max_y = temp_env->getMaxY();
	}

	const geos::geom::Envelope *env = new geos::geom::Envelope(min_x, max_x, min_y, max_y);

	return env;
}*/

void Util_Quick_sort_double(double *arr, ulong low, ulong high)
{
	/* pi for partitioning index, put arr[pi] to right group */
	ulong pi;

	if (low < high)
	{
		pi = Util_Quick_sort_partition(arr, low, high);

		if(pi > 0)
			Util_Quick_sort_double(arr, low, pi - 1);

		Util_Quick_sort_double(arr, pi+1, high);
	}

	return;
}

ulong Util_Quick_sort_partition(double *arr, ulong low, ulong high)
{
	ulong return_index, i;
	/* pivot to be the last elements*/
	double temp_val, pivot = arr[high];

	/* index of smaller element and indicats the right position of pivot*/
	return_index = low - 1;

	for (i = low; i < high; ++i)
	{
		// If curent element is smaller than the pivot
		if (arr[i] < pivot)
		{
			++return_index;

			/* swap */
			temp_val = arr[return_index];
			arr[return_index] = arr[i];
			arr[i] = temp_val;
		}
	}

	++return_index;
	/* put pivot at return_index */
	arr[high] = arr[return_index];
	arr[return_index] = pivot;

	return return_index;
}

void Util_Merge_sort_double(double *arr, ulong low, ulong high)
{
	ulong mid;

	if (low < high)
	{
		mid = low + (high - low) / 2;

		Util_Merge_sort_double(arr, low, mid);
		Util_Merge_sort_double(arr, mid+1, high);

		Util_Merge_sort_merge(arr, low, mid, high);
	}
}

void Util_Merge_sort_merge(double *arr, ulong low, ulong mid, ulong high)
{
	ulong i, j, k, size_left, size_right;

	size_left = mid - low + 1;
	size_right = high - mid;

	double arr_left[size_left], arr_right[size_right];

	for (i = 0; i < size_left; ++i)
		arr_left[i] = arr[low+i];

	for (i = 0; i < size_right; ++i)
		arr_right[i] = arr[mid+1+i];

	i = 0;
	j = 0;
	k = low;

	while (i < size_left && j < size_right)
	{
		if(arr_left[i] <= arr_right[j]) 
		{
			arr[k] = arr_left[i];
			i++;
		}
		else
		{
			arr[k] = arr_right[j];
			j++;
		}

		k++;
	}

	//extra element in left array
	while (i < size_left) 
	{       
		arr[k] = arr_left[i];
		i++; 
		k++;
	}

	//extra element in right array
	while (j < size_right) 
	{     
		arr[k] = arr_right[j];
		j++; 
		k++;
	}
}

void Util_Read_geoms_mpi(std::string file_path, std::vector<std::string> *vect_strs, uint num_threads, MPI_Comm comm)
{
	int num_world_nodes, my_world_rank;

	MPI_Comm_size(comm, &num_world_nodes);
	MPI_Comm_rank(comm, &my_world_rank);

	MPI_Info myinfo;
	MPI_Info_create(&myinfo);
	MPI_Info_set(myinfo, "access_style", "read_once,sequential"); 
	MPI_Info_set(myinfo, "collective_buffering", "true"); 
	MPI_Info_set(myinfo, "romio_cb_read", "enable");
  
	MPI_Offset file_size;
	MPI_Offset local_size;
	MPI_Offset start;
	MPI_Offset end;
		
	MPI_File fh;

	MPI_File_open(comm, file_path.c_str(), MPI_MODE_RDONLY, myinfo, &fh);

	MPI_File_get_size(fh, &file_size);
	
	local_size = file_size/num_world_nodes;
	start = my_world_rank * local_size;
	end = start + local_size - 1;

    //last process need to reach the end of file
    if (num_world_nodes - 1 == my_world_rank)
	{
        end = file_size - 1;
    }
    
    MPI_Offset chunk_size = end - start + 1;

    //Allocate memory for main chunk
    char * chunk;
    chunk = (char *)malloc((chunk_size + 1) * sizeof(char));

    if (NULL == chunk)
	{
        std::cerr<<"Error in malloc for file chunk"<<std::endl;
        exit(1);
    }
    //To store incomplete strings
    char *front_chunk;
    char *back_chunk;

    MPI_File_read_at(fh, start, chunk, chunk_size, MPI_CHAR, MPI_STATUS_IGNORE);

    chunk[chunk_size] = '\0';
    
    //Get last string which should be incomplete
    MPI_Offset valid_end = chunk_size;

    if (num_world_nodes - 1 != my_world_rank)
	{
        for ( ; valid_end >=0; --valid_end)
		{
            if (chunk[valid_end] == '\n')
			{
                back_chunk = (char *)malloc((chunk_size - valid_end) * sizeof(char));
                strncpy(back_chunk, chunk+valid_end+1, chunk_size - valid_end);
                break;
            }
        }
    }

    //Get first string which should be incomplete
    MPI_Offset valid_start = 0;
    if (0 != my_world_rank)
	{
        for ( ; valid_start < chunk_size; ++valid_start)
		{
            if (chunk[valid_start] == '\n')
			{
                front_chunk = (char*)malloc((valid_start+1)*sizeof(char));
                strncpy(front_chunk, chunk, valid_start);
                front_chunk[valid_start] = '\0';
                ++valid_start;
                break;
            }
        }
    }

	MPI_Offset str_start = valid_start;
	MPI_Offset str_end = valid_start;

	for ( ; valid_start <= valid_end; ++valid_start)
	{
		if(chunk[valid_start] == '\n' || valid_start == valid_end){
			str_end = valid_start;
			MPI_Offset str_len = str_end-str_start;

			char *temp_chunk;

			temp_chunk = (char *)malloc((str_len+1) * sizeof(char));
			strncpy(temp_chunk, chunk+str_start, str_len);

			temp_chunk[str_len] = '\0';	

			vect_strs->push_back(temp_chunk);

			free(temp_chunk);

			str_start = valid_start+1;
		}
	}
	
	char * recv_chunk;

	if (0 == my_world_rank%2)
	{
		if (num_world_nodes - 1 != my_world_rank)
		{
			MPI_Offset send_len = strlen(back_chunk);

			MPI_Send(&send_len, 1, MPI_OFFSET, my_world_rank+1, my_world_rank, MPI_COMM_WORLD);
			MPI_Send(back_chunk, send_len, MPI_CHAR, my_world_rank+1, my_world_rank+1, MPI_COMM_WORLD);
		}

 		if (0 != my_world_rank)
		{
			MPI_Offset recv_len = 0;
			MPI_Status * status = new MPI_Status();

			MPI_Recv(&recv_len, 1, MPI_OFFSET, my_world_rank-1, my_world_rank-1,MPI_COMM_WORLD, status);
			recv_chunk = (char *)malloc((recv_len) * sizeof(char));
			MPI_Recv(recv_chunk, recv_len, MPI_CHAR, my_world_rank-1, my_world_rank,MPI_COMM_WORLD, status);
		}
	}
	else
	{
		MPI_Offset recv_len = 0;
		MPI_Status * status = new MPI_Status();

		MPI_Recv(&recv_len, 1, MPI_OFFSET, my_world_rank-1, my_world_rank-1,MPI_COMM_WORLD, status);
		recv_chunk = (char *)malloc((recv_len) * sizeof(char));
		MPI_Recv(recv_chunk, recv_len, MPI_CHAR, my_world_rank-1, my_world_rank,MPI_COMM_WORLD, status);

 		if (num_world_nodes - 1 != my_world_rank)
		{
			MPI_Offset send_len = strlen(back_chunk);
			MPI_Send(&send_len, 1, MPI_OFFSET, my_world_rank+1, my_world_rank, MPI_COMM_WORLD);
			MPI_Send(back_chunk, send_len, MPI_CHAR, my_world_rank+1, my_world_rank+1, MPI_COMM_WORLD);
		}
	}

 	if(my_world_rank != 0)
	{
		std::string temp_str = std::string(recv_chunk) + std::string(front_chunk);
		vect_strs->push_back(temp_str);
	}

    MPI_File_close(&fh);
}