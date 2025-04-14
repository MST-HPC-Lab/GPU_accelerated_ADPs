#ifndef __GSJ_GLOBAL_HEADER_H_INCLUDE__
#define __GSJ_GLOBAL_HEADER_H_INCLUDE__

/* Global variables*/
#include <global_var.h>

/* c c++ headers */
#include <vector>
#include <fstream>
#include <list>
#include <vector>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <ctime>
#include <cmath>
#include <thread>
#include <mutex>
#include <unordered_map>
#include <unistd.h>
#include <utility>
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono> // std::chrono::seconds, std::chrono::milliseconds

/* External libs */
#include <geos/geom/Coordinate.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Point.h>
#include <geos/io/WKTReader.h>
#include <geos/index/strtree/STRtree.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/geom/LinearRing.h>
#include <spdlog/spdlog.h>
#include <spdlog/cfg/env.h> // support for loading levels from the environment variable

/* Internal headers */
#include <grid_utils.h>
#include <adp_grid.h>
#include <cuda_host_interface.cuh>

/*Global Functions*/


#endif //ndef __GSJ_GLOBAL_HEADER_H_INCLUDE__