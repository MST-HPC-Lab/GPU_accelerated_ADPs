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

geos::geom::Envelope *Util_Read_envs(std::string file_path, std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs)
{
	std::ifstream file(file_path.c_str());
	std::list<std::string> *list_strs = new std::list<std::string>;
	std::string temp_str;

	/* First line as global env*/
	std::getline(file, temp_str);
	std::vector<std::string> tokens;
	tokenize(tokens, temp_str);
	geos::geom::Envelope *global_env = new geos::geom::Envelope(std::stod(tokens[0]), std::stod(tokens[1]),
																std::stod(tokens[2]), std::stod(tokens[3]));
	tokens.clear();

	while (std::getline(file, temp_str))
	{
		//omit empty strings and invalid strings
		if (temp_str.size() > 5)
			list_strs->push_back(temp_str);
	}

	file.close();

	for (std::list<std::string>::iterator itr = list_strs->begin(); itr != list_strs->end(); ++itr)
	{
		temp_str = *itr;
		tokens.clear();
		tokenize(tokens, temp_str);

		geos::geom::Envelope *temp_env = new geos::geom::Envelope(std::stod(tokens[0]), std::stod(tokens[1]),
																  std::stod(tokens[2]), std::stod(tokens[3]));

		ulong temp_weight = std::stoul(tokens[4]);

		list_envs->push_back(std::make_pair(temp_env, temp_weight));
	}

	delete list_strs;
	return global_env;
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
		tokens.clear();
		tokenize(tokens, temp_str);

		geos::geom::Envelope *temp_env = new geos::geom::Envelope(std::stod(tokens[0]), std::stod(tokens[1]),
																std::stod(tokens[2]), std::stod(tokens[3]));

		ulong temp_weight = std::stoul(tokens[4]);

		temp_list_envs->push_back(std::make_pair(temp_env, temp_weight));
	}

	push_mutex->lock();
	for (auto &itr : *temp_list_envs)
		list_envs->push_back(itr);
	push_mutex->unlock();
}
geos::geom::Envelope *Util_Read_envs_parallel(std::string file_path, std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs, uint num_threads)
{
	std::ifstream file(file_path.c_str());
	std::vector<std::string> *vect_strs = new std::vector<std::string>;
	std::string temp_str;

	/* First line as global env*/
	std::getline(file, temp_str);
	std::vector<std::string> tokens;
	tokenize(tokens, temp_str);
	geos::geom::Envelope *global_env = new geos::geom::Envelope(std::stod(tokens[0]), std::stod(tokens[1]),
																std::stod(tokens[2]), std::stod(tokens[3]));
	tokens.clear();

	while (std::getline(file, temp_str))
	{
		//omit empty strings and invalid strings
		if (temp_str.size() > 5)
			vect_strs->push_back(temp_str);
	}

	file.close();

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

		//std::thread t(parallel_read_env, begin, end, vect_strs, list_envs, push_mutex);
		//readers.push_back(std::move(t));
	}

	std::for_each(readers.begin(), readers.end(), [](std::thread &t) 
    {
        t.join();
    });

	delete vect_strs;
	return global_env;
}

int Util_Write_geoms_to_WKT(std::string file_name, const geos::geom::Envelope *env,
							geos::index::strtree::STRtree *index)
{
	//Finding geometries
	std::vector<void *> results;

	index->query(env, results);

	//do not produce empty files
	if (results.empty())
	{
		printf("EMPTY %s  %s\n", file_name, env->toString());
		return 0;
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
	return 0;
}

geos::geom::Geometry *Util_Covert_env_to_geom(const geos::geom::Envelope *env)
{
	if (env == NULL)
		return NULL;

	double min_x = env->getMinX();
	double max_x = env->getMaxX();
	double min_y = env->getMinY();
	double max_y = env->getMaxY();

	std::string temp_str = "POLYGON ((" + std::to_string(min_x) + " " + std::to_string(min_y) + "," + std::to_string(min_x) + " " + std::to_string(max_y) + "," + std::to_string(max_x) + " " + std::to_string(max_y) + "," + std::to_string(max_x) + " " + std::to_string(min_y) + "," + std::to_string(min_x) + " " + std::to_string(min_y) + +"))";

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

		geos::geom::CoordinateSequence* coord_seq = (temp_geom->getCoordinates()).release();

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

const geos::geom::Envelope *Util_Get_global_env(std::list<geos::geom::Geometry *> *layer_1_geoms,
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
}

std::list<std::pair<geos::geom::Coordinate *, ulong> > *Util_Get_adp_candidates(std::list<geos::geom::Geometry *> *list_geoms_1,
																				std::list<geos::geom::Geometry *> *list_geoms_2)
{
	std::list<std::pair<geos::geom::Coordinate *, ulong> > *candidates = new std::list<std::pair<geos::geom::Coordinate *, ulong> >();

	geos::index::strtree::STRtree index;

	for (std::list<geos::geom::Geometry *>::iterator itr = list_geoms_1->begin(); itr != list_geoms_1->end(); ++itr)
	{
		geos::geom::Geometry *temp_geom = *itr;
		index.insert(temp_geom->getEnvelopeInternal(), temp_geom);
	}

	for (std::list<geos::geom::Geometry *>::iterator itr = list_geoms_2->begin(); itr != list_geoms_2->end(); ++itr)
	{
		geos::geom::Geometry *temp_geom = *itr;

		std::vector<void *> results;

		const geos::geom::Envelope *temp_env = temp_geom->getEnvelopeInternal();

		index.query(temp_env, results);

		for (std::vector<void *>::iterator results_itr = results.begin(); results_itr != results.end(); ++results_itr)
		{
			void *geom_ptr = *results_itr;
			geos::geom::Geometry *queried_geom = (geos::geom::Geometry *)geom_ptr;

			const geos::geom::Envelope *queried_geom_env = queried_geom->getEnvelopeInternal();

			geos::geom::Envelope *intersection_env = new geos::geom::Envelope();

			temp_env->intersection(*queried_geom_env, *intersection_env);

			if (intersection_env != NULL)
			{

				candidates->push_back(std::make_pair(
					new geos::geom::Coordinate((intersection_env->getMaxX() + intersection_env->getMinX()) / 2.0,
											   (intersection_env->getMaxY() + intersection_env->getMinY()) / 2.0),
					(ulong)(temp_geom->getNumPoints() * queried_geom->getNumPoints())));
			}
		}
	}

	return candidates;
}

void Util_Convert_envs_to_array(std::list<std::pair<const geos::geom::Envelope *, ulong> > *list_envs, double *new_arr)
{
	assert(list_envs != nullptr && new_arr != nullptr);
	if (list_envs->empty())
		return;

	/* array for minx, miny, maxx, maxy and one weight */
	/* use double to store weights to have consistency and higher upper bound*/
	ulong size = list_envs->size();
	std::cout << size << std::endl;

	ulong counter = 0;

	for (auto &itr : *list_envs)
	{
		new_arr[counter] = itr.first->getMinX();
		new_arr[size + counter] = itr.first->getMinY();
		new_arr[2 * size + counter] = itr.first->getMaxX();
		new_arr[3 * size + counter] = itr.first->getMaxY();
		new_arr[4 * size + counter] = (double)(itr.second);
		++counter;
	}
}
