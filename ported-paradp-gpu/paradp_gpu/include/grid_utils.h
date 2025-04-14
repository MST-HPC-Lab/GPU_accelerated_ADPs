#ifndef __GRID_UTILS_H_INCLUDED__
#define __GRID_UTILS_H_INCLUDED__

/* Global variables*/
#include <global_var.h>

/* c c++ headers */
#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <list>
#include <string.h>
#include <cmath>
#include <thread> // std::thread
#include <mutex>  // std::mutex

/* External libs */
#include <geos/geom/Geometry.h>
#include <geos/io/WKTReader.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/index/strtree/STRtree.h>
#include <geos/geom/PrecisionModel.h>

#include <mpi.h>

/* Internal libs */
/* None */

void Util_Read_geoms(std::string file_path, std::list<geos::geom::Geometry *> *list_geoms);

void Util_Read_geoms_parallel(std::string file_path, std::list<geos::geom::Geometry *> *list_geoms, uint num_threads);

void Util_Write_geoms_to_WKT(std::string file_name, const geos::geom::Envelope *env,
							 geos::index::strtree::STRtree *index);

geos::geom::Geometry *Util_Covert_env_to_geom(const geos::geom::Envelope *env);

void Util_Write_grid_to_WKT(std::string file_name, std::list<const geos::geom::Envelope *> *list_envs);

void Util_Write_geoms_to_array(std::string file_name, const geos::geom::Envelope *env,
							   geos::index::strtree::STRtree *index);

//const geos::geom::Envelope *Util_Get_global_env(std::list<geos::geom::Geometry *> *layer_1_geoms,
//												std::list<geos::geom::Geometry *> *layer_2_geoms);

const geos::geom::Envelope *Util_Get_global_env(std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs_1,
												std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs_2);												

void Util_Quick_sort_double(double *arr, ulong low, ulong high);

ulong Util_Quick_sort_partition(double *arr, ulong low, ulong high);

void Util_Merge_sort_double(double *arr, ulong low, ulong high);

void Util_Merge_sort_merge(double *arr, ulong low, ulong mid, ulong high);

void Util_Read_geoms_mpi(std::string file_path, std::vector<std::string> *vect_strs, uint num_threads, MPI_Comm comm);

void Util_Parse_string_parallel(std::vector<std::string> *vect_strs, std::list<geos::geom::Geometry *> *list_geoms, uint num_threads);

void Util_Parse_env_string_parallel(std::vector<std::string> *vect_strs, std::list<std::pair<geos::geom::Envelope *, ulong> > *list_envs, uint num_threads);

void Util_Parse_vector_envs_to_mbr_arr(std::vector<std::pair<const geos::geom::Envelope *, ulong> *> *vector_envs, MBR_Array &mbr_arr);

#endif //ndef __GRID_UTILS_H_INCLUDED__