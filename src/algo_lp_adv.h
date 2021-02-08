/**
The module contains the code concerning the fully adversary algorithm on GPS point
**/
#ifndef __KCENTER_LP_ADV_HEADER__
#define __KCENTER_LP_ADV_HEADER__

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
} Lp_cluster;

void lp_initialise_level(Lp_cluster * level, unsigned int k, void * array, unsigned int z);

void lp_delete_level(Lp_cluster * clusters);

void lp_compute_true_radius(Lp_cluster * level);

void lp_k_center_add(Lp_cluster * level, unsigned int index);

void lp_k_center_delete(Lp_cluster * level, unsigned int index);

void lp_k_center_solver(Lp_cluster *level, double epsilon, double dmin, double dmax);

void lp_k_center_run(Lp_cluster *level, struct query_provider * queries, double epsilon, double dmin, double dmax);
#endif
