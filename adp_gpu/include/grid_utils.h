#ifndef __GRID_UTILS_H_INCLUDED__
#define __GRID_UTILS_H_INCLUDED__

#include "global_var.h"

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <list>
#include <cmath>
#include <utility>      // std::pair
#include <mutex>
#include <thread>

#include <geos/geom/Geometry.h>
#include <geos/io/WKTReader.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/index/strtree/STRtree.h>
#include <geos/geom/PrecisionModel.h>

void Util_Read_geoms(std::string file_path, std::list<geos::geom::Geometry *> *list_geoms);

int Util_Write_geoms_to_WKT(std::string file_name, const geos::geom::Envelope *env,
							geos::index::strtree::STRtree *index);

geos::geom::Geometry *Util_Covert_env_to_geom(const geos::geom::Envelope *env);

void Util_Write_grid_to_WKT(std::string file_name, std::list<const geos::geom::Envelope *> *list_envs);

void Util_Write_geoms_to_array(std::string file_name, const geos::geom::Envelope *env,
							   geos::index::strtree::STRtree *index);

const geos::geom::Envelope *Util_Get_global_env(std::list<geos::geom::Geometry *> *layer_1_geoms,
												std::list<geos::geom::Geometry *> *layer_2_geoms);

std::list<std::pair<geos::geom::Coordinate *, ulong> > *Util_Get_adp_candidates(
	std::list<geos::geom::Geometry *> *list_geoms_1, std::list<geos::geom::Geometry *> *list_geoms_2);

void Util_Convert_envs_to_array(std::list<std::pair<const geos::geom::Envelope *, ulong> > *list_envs,
								double *new_arr);

geos::geom::Envelope *Util_Read_envs(std::string file_path,
									 std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs);

geos::geom::Envelope *Util_Read_envs_parallel(std::string file_path, std::list<std::pair<geos::geom::Envelope *, ulong> > 
								*list_envs, uint num_threads=36);

#endif //ndef __GRID_UTILS_H_INCLUDED__
