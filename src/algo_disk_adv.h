/**
The module contains the code concerning the fully adversary algorithm on GPS point
**/
#ifndef __KCENTER_DISK_ADV_HEADER__
#define __KCENTER_DISK_ADV_HEADER__

#include "point.h"
#include "data_fully_adv.h"
#include "query.h"

#include <stdint.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

typedef struct {
	unsigned int k;			/* Maximum number of cluster allowed */
	unsigned int z;			/* Maximum number of outliers allowed */
	double radius;			/* maximum cluster radius of current level */
	void *array;			/* pointer to all points */

	std::unordered_set<unsigned int> points; /* Currently inserted points */
	std::vector<unsigned int> selected_centers;			/* S, the selected k centers */
	double true_radius;		/* true radius */
} Disk_cluster;

void disk_initialise_level(Disk_cluster * level, unsigned int k, void * array, unsigned int z);

void disk_delete_level(Disk_cluster * clusters);

void disk_compute_true_radius(Disk_cluster * level);

void disk_k_center_add(Disk_cluster * level, unsigned int index);

void disk_k_center_delete(Disk_cluster * level, unsigned int index);

void disk_k_center_solver(Disk_cluster *level, double epsilon, double dmin, double dmax);

void disk_k_center_run(Disk_cluster *level, struct query_provider * queries, double epsilon, double dmin, double dmax);
#endif
