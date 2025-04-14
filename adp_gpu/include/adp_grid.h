#ifndef __ADP_GRID_H_INCLUDED__
#define __ADP_GRID_H_INCLUDED__

#include <global_var.h>

#include <stdio.h>
#include <stdlib.h>

#include <list>
#include <cmath>
#include <limits.h>

#include <geos/geom/Envelope.h>
#include <geos/geom/Geometry.h>
#include <geos/index/strtree/STRtree.h>

#define MAX_SIZE LLONG_MAX

class Adp_Grid
{
	public:
	Adp_Grid(int cells, const geos::geom::Envelope *universe, std::list<const geos::geom::Envelope*>* basis)
	{
		num_cells = cells;
		grid_cells_with_candidates_weight = NULL;

		grid_cells_with_candidates_num = new std::list<std::pair<const geos::geom::Envelope*, ulong>>; 

		partition(basis, universe);
	}

	Adp_Grid(int cells, const geos::geom::Envelope *universe, std::list<std::pair<const geos::geom::Envelope*, ulong> >* basis)
	{
		num_cells = cells;
		grid_cells_with_candidates_num = NULL;

		grid_cells_with_candidates_weight = new std::list<std::pair<const geos::geom::Envelope*, ulong>>; 

		build_rtree_index(basis); 

		partition_using_envs_weight(basis, universe);
	}

	Adp_Grid(int cells, const geos::geom::Envelope *universe, std::list<std::pair<geos::geom::Coordinate*, ulong > > * basis)
	{
		num_cells = cells;
		grid_cells_with_candidates_num = NULL;
		
		//totalComplexity = compute_total_complexity(basis);

		grid_cells_with_candidates_weight = new std::list<std::pair<const geos::geom::Envelope*, ulong>>; 
		
		build_rtree_index(basis); 

		partition_using_coords_weight(basis, universe);
	} 
 
	std::list<std::pair<const geos::geom::Envelope*, ulong> > * get_grid_with_candidates_num(void){
		return grid_cells_with_candidates_num;
	}

	std::list<std::pair<const geos::geom::Envelope*, ulong> > * get_grid_with_candidates_weight(void){
		return grid_cells_with_candidates_weight;
	}

	private:
	uint num_iterations;
	uint num_cells;
	std::list<const geos::geom::Envelope *> *grid_cells;
	geos::index::strtree::STRtree index;

	std::list<std::pair<const geos::geom::Envelope*, ulong> > * grid_cells_with_candidates_num;
	std::list<std::pair<const geos::geom::Envelope*, ulong> > * grid_cells_with_candidates_weight;
		
	void build_rtree_index(std::list<std::pair<geos::geom::Coordinate*, ulong>>* basis);

	void build_rtree_index(std::list<std::pair<const geos::geom::Envelope*, ulong> >* basis);

	ulong compute_total_complexity(std::list<std::pair<geos::geom::Coordinate*, ulong>>* basis);
	
	int partition(std::list<const geos::geom::Envelope *> * basis, const geos::geom::Envelope * universe);

	int partition_using_envs_weight(std::list<std::pair<const geos::geom::Envelope*, ulong> >* basis, 
			const geos::geom::Envelope *universe);

	int partition_using_coords_weight(std::list<std::pair<geos::geom::Coordinate*, ulong> >* basis, 
			const geos::geom::Envelope *universe);

	int partition_to_4cells(std::list<const geos::geom::Envelope*>* basis, const geos::geom::Envelope * env);

	int partition_to_4cells_using_envs_weight(std::list<std::pair<const geos::geom::Envelope*, ulong> >* basis, 
			const geos::geom::Envelope * env);

	int partition_to_4cells_using_coords_weight(std::list<std::pair<geos::geom::Coordinate*, ulong> >* basis, 
			const geos::geom::Envelope * env);

	ulong get_num_geoms_in_cell(std::list<const geos::geom::Envelope*>* basis, const geos::geom::Envelope* env);

	ulong get_complexity_for_cell(std::list<std::pair<const geos::geom::Envelope*, ulong>>* basis, const geos::geom::Envelope* env);

	ulong get_complexity_for_cell(std::list<std::pair<geos::geom::Coordinate*, ulong>>* basis, const geos::geom::Envelope* env); 
};

#endif //ndef __ADP_GRID_H_INCLUDED__
