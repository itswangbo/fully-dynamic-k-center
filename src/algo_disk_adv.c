/**
The module contains the code concerning the fully adversary algorithm on GPS point
**/
#include "utils.h"
#include "query.h"
#include "data_fully_adv.h"
#include "algo_disk_adv.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <float.h>

#include <algorithm>
#include <set>
#include <map>
#include <vector>
#include <unordered_set>

void
disk_initialise_level(Disk_cluster * level, unsigned int k,
			   void * array, unsigned int z)
{
	level->k = k;
	level->z = z;

	level->array = array;
}

void disk_delete_level(Disk_cluster * level)
{
	level->selected_centers.clear();
	level->points.clear();

	level->array = NULL;

	delete level;
}

/* insertion */
void disk_k_center_add(Disk_cluster * level, unsigned int index)
{
	level->points.insert(index);
}

/* deletion */
void disk_k_center_delete(Disk_cluster * level, unsigned int index)
{
	level->points.erase(index);
}

/* sy_relation: pair<j, B_j> d(i,j)<=radius */
void construct_sy_relation(Disk_cluster *level, double radius, std::unordered_map<unsigned int, std::unordered_set<unsigned int> > *sy_relation)
{
	double tmp;
	std::unordered_set<unsigned int>::iterator iter=level->points.begin();
	std::vector<unsigned int> inserted_points;
	std::unordered_set<unsigned int> cur;

	sy_relation->clear();

	while(iter != level->points.end()){
		for(int i=0; i<inserted_points.size(); i++){
			tmp=fully_adv_distance(((Geo_point *)level->array) + inserted_points.at(i),
				((Geo_point *)level->array) + *iter);
			if(tmp <= radius){
				cur.insert(inserted_points.at(i));
				(*sy_relation)[inserted_points.at(i)].insert(*iter);
			}
		}
		sy_relation->insert(std::make_pair(*iter, cur));

		cur.clear();

		inserted_points.push_back(*iter);
		iter++;
	}
}

//delete @sel_center, record the centers whose bhat value are infected and update expanded symmetric relation
void select_update_card_and_syrelation(Disk_cluster *level, std::unordered_map<unsigned int, std::unordered_set<unsigned int> >* sy_relation, 
				std::unordered_map<unsigned int, unsigned int>* modify_card, unsigned int sel_center, std::unordered_set<unsigned int>* old_points,  
				std::set<std::pair<unsigned int, unsigned int>, std::greater<std::pair<unsigned int, unsigned int> > >* cardinality,
				std::unordered_map<unsigned int, unsigned int>* card_helper)
{
	(*old_points).erase(sel_center);

	unsigned int temp_card=(*card_helper)[sel_center];
	(*cardinality).erase(std::make_pair(temp_card, sel_center));
	(*card_helper).erase(sel_center);

	if((*sy_relation).find(sel_center)!=(*sy_relation).end()){
		std::unordered_set<unsigned int>::iterator neighbors =(*sy_relation)[sel_center].begin();
		while(neighbors != (*sy_relation)[sel_center].end()){
			if((*modify_card).find(*neighbors) != (*modify_card).end()){
				(*modify_card)[*neighbors]++;
			}else (*modify_card)[*neighbors]=1;

			(*sy_relation)[*neighbors].erase(sel_center);
			if((*sy_relation)[*neighbors].empty()) (*sy_relation).erase(*neighbors);
			neighbors++;
		}
		(*sy_relation).erase(sel_center);
	}
	return;
}

/* return @1: can cover enough points; @0: cannot cover enough points */
int check_points_cover(Disk_cluster* level, std::unordered_map<unsigned int, std::unordered_set<unsigned int> > sy_relation)
{
	/* pair<cardinality, points>, in decreasing order of cardinality */
	std::set<std::pair<unsigned int, unsigned int>, std::greater<std::pair<unsigned int, unsigned int> > > cardinality;
	/* map<points, cardinality>, used as index to find element in cardinality */	
	std::unordered_map<unsigned int, unsigned int> card_helper;

	/* map<points, modified cardinality>, record the change of cardinality of each point while updating sy_relation */	
	std::unordered_map<unsigned int, unsigned int> modify_card;
	
	std::unordered_map<unsigned int, std::unordered_set<unsigned int> >::iterator iter=sy_relation.begin();
	std::unordered_set<unsigned int>::iterator iter_old_centers;
	std::unordered_set<unsigned int>::iterator iter_center;

	std::unordered_set<unsigned int> old_points;
	std::unordered_set<unsigned int> temp_old_points;

	double tmp;

	/* construct cardinality */
	while(iter != sy_relation.end()){
		/* initialize uncovered points */
		old_points.insert(iter->first);

		cardinality.insert(std::make_pair(iter->second.size()+1, iter->first));
		card_helper[iter->first] = iter->second.size()+1;

		iter++;
	}

	int num_cur_points=level->points.size();

	level->selected_centers.clear();

	/* pick centers */
	int num_picked_centers=0, num_centers=level->k, num_picked_points=0, size, cur_card;
	unsigned int picked_center;
	while(!old_points.empty() && num_picked_centers<num_centers){
		picked_center=cardinality.begin()->second;
		level->selected_centers.push_back(picked_center);
		num_picked_centers++;

		select_update_card_and_syrelation(level, &sy_relation, &modify_card, picked_center, &old_points, &cardinality, &card_helper);
		num_picked_points++;
	
		if(!old_points.empty()){
			temp_old_points=old_points;

			iter_old_centers=old_points.begin();			
			while(iter_old_centers != old_points.end()){
				tmp = fully_adv_distance(((Geo_point *)level->array) + picked_center,
					((Geo_point *)level->array) + *iter_old_centers);
				if(tmp<=3*(level->radius)){

					select_update_card_and_syrelation(level, &sy_relation, &modify_card, *iter_old_centers, 
										&temp_old_points, &cardinality, &card_helper);

					num_picked_points++;
				}

				iter_old_centers++;
			}

			old_points = temp_old_points;
			temp_old_points.clear();
		}
		if(!old_points.empty()){
			iter_center = old_points.begin();
			while(iter_center != old_points.end()){
				if(modify_card.find(*iter_center) != modify_card.end()){

					cur_card=card_helper[*iter_center];
					cardinality.erase(std::make_pair(cur_card, *iter_center));
				
					cur_card-=modify_card[*iter_center];
					cardinality.insert(std::make_pair(cur_card, *iter_center));
					card_helper[*iter_center]=cur_card;
				}
				iter_center++;
			}
		}
	}

	if(num_picked_points>=num_cur_points-level->z) return 1;
	else return 0;
}

/* select k centers: Binary search r */
void disk_k_center_solver(Disk_cluster *level, double eps, double dmin, double dmax)
{
	/* coupute x_value and y_value by solving linear program */
	int cur_points=level->points.size();

	/* if the number of current points is no larger than k plus outliers */
	/* pick any k points */
	if(cur_points<=(level->k+level->z)){
        	level->radius=0;
        	level->selected_centers.clear();
        	std::unordered_set<unsigned int>::iterator iter=level->points.begin();
        	for(int i=0; i<level->k && i<cur_points; i++){
            		level->selected_centers.push_back(*iter);
            		iter++;
        	}
        	return;
	}
	/* end of special case */

	std::unordered_map<unsigned int, std::unordered_set<unsigned int> > sy_relation; /* pair<j, B_j> d(i,j)<=radius */

	int res;

	//double radius = pow( (1+eps), ceil(log(dmin)/log(1+eps)) );
	double radius = dmin;
	unsigned int lef=0, mid, rig=(unsigned int)(1 + ceil(log(dmax / dmin) / log(1 + eps)));
	while(lef<rig){
		mid=lef+(rig-lef)/2;
		level->radius = radius*pow((1+eps),mid);

		construct_sy_relation(level, (level->radius), &sy_relation);

        	res = check_points_cover(level, sy_relation);

		if(res == 1){
			rig=mid;
		}else{
			lef=mid+1;
		}
	}

	if(res != 1){
        	level->radius = radius*pow((1+eps),lef);

        	construct_sy_relation(level, level->radius, &sy_relation);

        	check_points_cover(level, sy_relation);
	}
}

/* compute true_radius */
void disk_compute_true_radius(Disk_cluster * level)
{
    unsigned int cur_points = level->points.size(), outliers=level->z, num_centers=level->k;
    if(cur_points<=outliers+num_centers){
        level->true_radius=0;
        return;
    }

	double k_rad, tmp;
	unsigned int nb_centers;
	std::unordered_set<unsigned int>::iterator iter=level->points.begin();

	/* sort the distance between current point and the centers in increasing order */
	std::vector<double> top_dis;
	nb_centers = level->selected_centers.size();

	while(iter != level->points.end()){
		k_rad = DBL_MAX; //Maximum double
		for(int i=0; i<nb_centers; i++){
			tmp = fully_adv_distance(((Geo_point *)level->array) + *iter,
						((Geo_point *)level->array) + level->selected_centers.at(i));
			k_rad = MIN(k_rad, tmp);
		}
		top_dis.push_back(k_rad);

		iter++;
	}

	std::sort(top_dis.begin(), top_dis.end());

	/* true radius should be the-number-of-outliers-th element */
	level->true_radius = top_dis.at(cur_points-outliers-1);
}

Error_enum fully_adv_write_log(Disk_cluster* level)
{
	if (has_log()) {
		fprintf(get_log_file(), "true radius: %lf\n", level->true_radius);
	}
	return NO_ERROR;
}

Error_enum disk_apply_one_query(Disk_cluster* level, struct query * query, unsigned int& N)
{
	if (query->type == ADD) {
		/* printf("a %u\n", query->data_index); */
		N++;
		disk_k_center_add(level, query->data_index);
	} else {
		N--;
		disk_k_center_delete(level, query->data_index);
	}

	return NO_ERROR;
	//return fully_adv_write_log(level);
}

void disk_k_center_run(Disk_cluster * level, struct query_provider * queries, double epsilon, double dmin, double dmax)
{
	struct query query;

	/* record the number of active points */
	unsigned int count = 0;

	long int initial_memory_usage = memory_usage();

	struct timeval start, end;
	double timeuse;

	printf("Run starts: memory usage %8ld KB\n", initial_memory_usage);
	while (get_next_query_set_lp(queries, &query, level->points)) {
		disk_apply_one_query(level, &query, count);
	}
	gettimeofday(&start, NULL);
	disk_k_center_solver(level, epsilon, dmin, dmax);
	gettimeofday(&end, NULL);
	
	timeuse = (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0;

	disk_compute_true_radius(level);
	fully_adv_write_log(level);
}
