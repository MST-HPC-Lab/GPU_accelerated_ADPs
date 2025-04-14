#include <adp_grid.h>

/**********************************************Normal Quadtree Partition****************************************************************************/

int Adp_Grid :: partition(std::list<const geos::geom::Envelope*>* basis, const geos::geom::Envelope *universe)
{
	grid_cells_with_candidates_num->push_back(std::pair<const geos::geom::Envelope*, ulong>(universe, (ulong)(basis->size())));

	while(grid_cells_with_candidates_num->size() < num_cells){
		partition_to_4cells(basis, grid_cells_with_candidates_num->back().first);		
	}
  return 0;
}

int Adp_Grid :: partition_to_4cells(std::list<const geos::geom::Envelope*>* basis, const geos::geom::Envelope * env)
{
	std::list<const geos::geom::Envelope*>* list_envs= new std::list<const geos::geom::Envelope*>;
	double min_x = env->getMinX();
	double min_y = env->getMinY();
	double max_x = env->getMaxX();
	double max_y = env->getMaxY();
	
	
	double hight = env->getHeight();
	double width = env->getWidth();
	//calculate centre
	double mid_x = min_x + width/2;
	double mid_y = min_y + hight/2;
	
	list_envs->push_back(new const geos::geom::Envelope(min_x, mid_x, min_y, mid_y));
	list_envs->push_back(new const geos::geom::Envelope(min_x, mid_x, mid_y, max_y));
	list_envs->push_back(new const geos::geom::Envelope(mid_x, max_x, min_y, mid_y));
	list_envs->push_back(new const geos::geom::Envelope(mid_x, max_x, mid_y, max_y));
	
	for(std::list<const geos::geom::Envelope*>::iterator itr = list_envs->begin() ; itr != list_envs->end(); itr++){
			const geos::geom::Envelope* temp_env = *itr;
			std::pair<const geos::geom::Envelope*, ulong> temp_pair (temp_env, 0);

			temp_pair.second = get_num_geoms_in_cell(basis, temp_env);

			if(temp_pair.second > 0)
				grid_cells_with_candidates_num->push_front(temp_pair);			
	}
	//pop the cell being parted
	grid_cells_with_candidates_num->pop_back();
	
	//sort by the nums of shapes in cells: from least to most
	grid_cells_with_candidates_num->sort([](const std::pair<const geos::geom::Envelope*,int> & a, const std::pair<const geos::geom::Envelope*,int> & b) 
			{ return a.second < b.second; });
	
	return 0;
}

ulong Adp_Grid :: get_num_geoms_in_cell(std::list<const geos::geom::Envelope*>* basis, const geos::geom::Envelope* env)
{
	ulong num_pairs = 0;
	ulong counter = 0;

	for (std::list<const geos::geom::Envelope*>::iterator itr = basis->begin() ; itr != basis->end(); itr++) {
		const geos::geom::Envelope* candidate_env = *itr;	

		if(env->covers(candidate_env) || env->intersects(candidate_env)){
			++num_pairs;
			++counter;
		}	
	}		
	if(counter < basis->size()/num_cells) 
			num_pairs=1;

	return num_pairs;
}
/**********************************************End Normal Quadtree Partition************************************************************************/

/**********************************************Quadtree Partition Consider geoms number, using coords***********************************************/

void Adp_Grid :: build_rtree_index(std::list<std::pair<geos::geom::Coordinate*, ulong>>* basis)
{
	for (std::list<std::pair<geos::geom::Coordinate*,ulong>>::iterator itr = basis->begin(); itr != basis->end(); ++itr)
	{
		std::pair<geos::geom::Coordinate*, ulong> temp_pair = *itr;
		index.insert(new const geos::geom::Envelope(*temp_pair.first),&temp_pair.second);
	}

	return;
}

int Adp_Grid :: partition_using_coords_weight(std::list<std::pair<geos::geom::Coordinate*, ulong>>* basis, const geos::geom::Envelope *universe)
{
	grid_cells_with_candidates_weight->push_back(std::pair<const geos::geom::Envelope*, ulong>(universe, MAX_SIZE));

	while(grid_cells_with_candidates_weight->size() < num_cells){
		partition_to_4cells_using_coords_weight(basis, grid_cells_with_candidates_weight->back().first);		
	}
  return 0;
}

int Adp_Grid :: partition_to_4cells_using_coords_weight(std::list<std::pair<geos::geom::Coordinate*, ulong>>* basis, const geos::geom::Envelope* env)
{
	std::list<const geos::geom::Envelope*>* list_envs= new std::list<const geos::geom::Envelope*>();

	double min_x = env->getMinX();
	double min_y = env->getMinY();
	double max_x = env->getMaxX();
	double max_y = env->getMaxY();
		
	double hight = env->getHeight();
	double width = env->getWidth();
	//calculate centre
	double mid_x = min_x + width/2;
	double mid_y = min_y + hight/2;
	
	list_envs->push_back(new const geos::geom::Envelope(min_x, mid_x, min_y, mid_y));
	list_envs->push_back(new const geos::geom::Envelope(min_x, mid_x, mid_y, max_y));
	list_envs->push_back(new const geos::geom::Envelope(mid_x, max_x, min_y, mid_y));
	list_envs->push_back(new const geos::geom::Envelope(mid_x, max_x, mid_y, max_y));
	
	for(std::list<const geos::geom::Envelope*>::iterator itr = list_envs->begin() ; itr != list_envs->end(); itr++){
			const geos::geom::Envelope* temp_env = *itr;

			std::pair<const geos::geom::Envelope*, ulong> temp_pair (temp_env, 0);
			temp_pair.second = get_complexity_for_cell(basis, temp_env);

			if(temp_pair.second!=0)
				grid_cells_with_candidates_weight->push_front(temp_pair);			
	}

	//pop the cell being parted
	grid_cells_with_candidates_weight->pop_back();
	
	grid_cells_with_candidates_weight->sort([](const std::pair<const geos::geom::Envelope*, ulong> & a, const std::pair<const geos::geom::Envelope*, ulong> & b) 
			{ return a.second < b.second; });
	
	return 0;
}

ulong Adp_Grid :: get_complexity_for_cell(std::list<std::pair<geos::geom::Coordinate*, ulong>>* basis, const geos::geom::Envelope* env)
{
	ulong complexity = 0;
	
	std::vector<void *> results;

	index.query(env, results);

	for(std::vector<void *>::iterator itr = results.begin(); itr != results.end(); ++itr){
		void* temp_void = *itr;

		ulong* temp_weight = (ulong*)temp_void;
		
		if(LLONG_MAX - complexity < *temp_weight){
		/* A cell has reach the maximum weight that a ulong can hold, usually this cell should be divided again*/
			complexity = LLONG_MAX;
			return complexity;
		}
		else
		{
			complexity += *temp_weight;
		}
	}

	return complexity;
}
/******************************************End Quadtree Partition Consider geoms number, using coords***********************************************/

/**********************************************Quadtree Partition Consider geoms number, using envs*************************************************/

void Adp_Grid :: build_rtree_index(std::list<std::pair<const geos::geom::Envelope*, ulong>>* basis)
{
	for (std::list<std::pair<const geos::geom::Envelope*,ulong> >::iterator itr = basis->begin(); itr != basis->end(); ++itr)
	{
		std::pair<const geos::geom::Envelope*, ulong> temp_pair = *itr;
		index.insert(temp_pair.first,&temp_pair.second);
	}

	return;
}

int Adp_Grid :: partition_using_envs_weight(std::list<std::pair<const geos::geom::Envelope*, ulong>>* basis, const geos::geom::Envelope *universe)
{
	grid_cells_with_candidates_weight->push_back(std::pair<const geos::geom::Envelope*, ulong>(universe, MAX_SIZE));

	while(grid_cells_with_candidates_weight->size() < num_cells){
		partition_to_4cells_using_envs_weight(basis, grid_cells_with_candidates_weight->back().first);		
	}

	return 0;
}

int Adp_Grid :: partition_to_4cells_using_envs_weight(std::list<std::pair<const geos::geom::Envelope*, ulong> >* basis, const geos::geom::Envelope * env)
{
	std::list<const geos::geom::Envelope*>* list_envs= new std::list<const geos::geom::Envelope*>;
	double min_x = env->getMinX();
	double min_y = env->getMinY();
	double max_x = env->getMaxX();
	double max_y = env->getMaxY();
	
	
	double hight = env->getHeight();
	double width = env->getWidth();
	//calculate centre
	double mid_x = min_x + width/2;
	double mid_y = min_y + hight/2;
	
	list_envs->push_back(new const geos::geom::Envelope(min_x, mid_x, min_y, mid_y));
	list_envs->push_back(new const geos::geom::Envelope(min_x, mid_x, mid_y, max_y));
	list_envs->push_back(new const geos::geom::Envelope(mid_x, max_x, min_y, mid_y));
	list_envs->push_back(new const geos::geom::Envelope(mid_x, max_x, mid_y, max_y));
	
	for(std::list<const geos::geom::Envelope*>::iterator itr = list_envs->begin() ; itr != list_envs->end(); itr++){
			const geos::geom::Envelope* temp_env = *itr;

			std::pair<const geos::geom::Envelope*, ulong> temp_pair (temp_env, 0);

			temp_pair.second = get_complexity_for_cell(basis, temp_env);

			if(temp_pair.second!=0)
				grid_cells_with_candidates_weight->push_front(temp_pair);			
	}

	//pop the cell being parted
	grid_cells_with_candidates_weight->pop_back();
	
	//sort by the nums of shapes in cells: from least to most
	grid_cells_with_candidates_weight->sort([](const std::pair<const geos::geom::Envelope*,ulong> & a, const std::pair<const geos::geom::Envelope*,ulong> & b) 
			{ return a.second < b.second; });

	return 0;
	
}

ulong Adp_Grid :: get_complexity_for_cell(std::list<std::pair<const geos::geom::Envelope*, ulong>>* basis, const geos::geom::Envelope* env)
{
	ulong complexity = 0;
	
	std::vector<void *> results;

	index.query(env, results);

	for(std::vector<void *>::iterator itr = results.begin(); itr != results.end(); ++itr){
		void* temp_void = *itr;

		ulong* temp_weight = (ulong*)temp_void;
		
		if(LLONG_MAX - complexity < *temp_weight){
		/* A cell has reach the maximum weight that a ulong can hold, usually this cell should be divided again*/
			complexity = LLONG_MAX;
			return complexity;
		}
		else
		{
			complexity += *temp_weight;
		}
	}

	return complexity;
}
/******************************************End Quadtree Partition Consider geoms number, using envs*************************************************/
