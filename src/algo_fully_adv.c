/**
The module contains the code concerning the fully adversary algorithm on GPS point
**/

#include "utils.h"
#include "query.h"
#include "data_fully_adv.h"
#include "algo_fully_adv.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>

#include <set>

static unsigned int NB_POINTS=0;

void
fully_adv_initialise_level(Fully_adv_cluster * level, unsigned int k,
			   double radius, void * array,
			   unsigned int nb_points, unsigned int cluster_size)
{
	level->nb = 0;
	level->covered_points=0;
	level->k = k;
	level->radius = radius;
	initialise_set_collection(&(level->clusters), k + 1, cluster_size,
				  nb_points);
	level->centers =
	    (unsigned int *)malloc_wrapper(sizeof(*(level->centers)) * (k + 1));
	level->true_rad = (double *)malloc_wrapper(sizeof(double) * k);
	level->nb_points = nb_points;
	level->array = array;
	
	return;
}

void fully_adv_delete_level(Fully_adv_cluster * level)
{
	level->sy_relation.clear();
	level->selected_centers.clear();
	level->selected_clusters.clear();	

	free_set_collection(&(level->clusters));
	free(level->centers);
	free(level->true_rad);
	level->centers = NULL;
	level->nb_points = 0;
	level->covered_points=0;
	level->array = NULL;

	return;
}

void fully_adv_initialise_level_array(Fully_adv_cluster * levels[],
				      unsigned int k, double eps, double d_min,
				      double d_max, unsigned int *nb_instances,
				      void *points,
				      unsigned int nb_points,
				      unsigned int cluster_size,
				      unsigned int *helper_array[])
{
	unsigned int i;
	unsigned int tmp;
	double radius;

	*nb_instances = tmp = (unsigned int)(1 + ceil(log(d_max / d_min) /
						      log(1 + eps)));
	radius = d_min;
	*levels = new Fully_adv_cluster[tmp];
	//*levels = (Fully_adv_cluster *) malloc_wrapper(sizeof(**levels) * tmp);
	*helper_array =
	    (unsigned int *)malloc_wrapper(sizeof(**helper_array) * nb_points);
	//fully_adv_initialise_level((*levels), k, z, 0, points, nb_points, cluster_size);

	for (i = 0; i < tmp; i++) {
		/* p=2*r */
		fully_adv_initialise_level(*(levels) + i, k, 2*radius, points,
					   nb_points, cluster_size);
		radius *= (1 + eps);
	}

	return;
}

void
fully_adv_delete_level_array(Fully_adv_cluster levels[],
			     unsigned int nb_instances,
			     unsigned int helper_array[])
{
	unsigned int i;
	for (i = 0; i < nb_instances; i++)
		fully_adv_delete_level(levels + i);
	//free(levels);
	delete[] levels;
	free(helper_array);

	return;
}

unsigned int fully_adv_get_index_smallest(Fully_adv_cluster levels[],
					  unsigned int nb_instances, unsigned int z)
{
	/* TO-DO: try binary search? */
	unsigned int i;
	for (i = 0; i < nb_instances; i++){
		/* if the number of covered points is larger than current points - outliers */
		if(levels[i].covered_points >= (NB_POINTS-z)){
			return i;
		}
	}
	return nb_instances;
}

//insert point @index to the last cluster and update the symmetric relation
void insert_last_center(Fully_adv_cluster *level, unsigned int index){
	double tmp;
	add_element_set_collection(&(level->clusters), index, level->nb);
	std::vector<unsigned int> cur;
	for(unsigned int j=0; j<level->nb; j++){
		tmp = fully_adv_distance(((Geo_point *)level->array) + index,
			 ((Geo_point *)level->array) + level->centers[j]);
		if(tmp<=3*level->radius) cur.push_back(j);
	}

	//update sy_relation
	for(unsigned int j=0; j<cur.size(); j++){
		level->sy_relation[level->nb].insert(cur.at(j));
		level->sy_relation[cur.at(j)].insert(level->nb);
	}
	//end of updating
	level->centers[level->nb] = index;
	level->true_rad[level->nb] = 0;
	level->nb++;

	return;
}

//delete all clusters after @cluster_index, and store deleted shuffled points in @helper_array, the number of points is @size
//update the symmetric relation by adding relations which still hold 
void add_syrelation_delete_current_center(Fully_adv_cluster *level, unsigned int cluster_index, unsigned int helper_array[], unsigned int *size){

	//remove all elements after current cluster
        remove_all_elements_after_set(&(level->clusters), cluster_index, helper_array, size);
	shuffle_array(helper_array, *size);

	std::unordered_set<unsigned int>::iterator iter;

	/* update sy_relation */
	if(!level->sy_relation.empty()){
		std::unordered_map<unsigned int, std::unordered_set<unsigned int> > temp_sy_relation;
		for(unsigned int j=0; j<cluster_index; j++){
			if(level->sy_relation.find(j) != level->sy_relation.end()){
				iter=level->sy_relation[j].begin();
				while(iter != level->sy_relation[j].end()){
					if(*iter < cluster_index) temp_sy_relation[j].insert(*iter);
					iter++;
				}
			}
		}
		level->sy_relation=temp_sy_relation;	
	}
	/* end of updating */

	level->nb = cluster_index;

	return;
}

//delete all clusters after @cluster_index, and store deleted shuffled points in @helper_array, the number of points is @size
//update the symmetric relation by checking the symmetric relation of the deleted centers and then delete the corresponding relation  
void check_syrelation_delete_current_center(Fully_adv_cluster *level, unsigned int cluster_index, unsigned int helper_array[], unsigned int *size){

	//remove all elements after current cluster
        remove_all_elements_after_set(&(level->clusters), cluster_index, helper_array, size);
	shuffle_array(helper_array, *size);

	std::unordered_set<unsigned int>::iterator iter;

	/* update sy_relation */
	if(!level->sy_relation.empty()){
		for(unsigned int j=cluster_index; j<level->nb; j++){
			if(level->sy_relation.find(j) != level->sy_relation.end()){
				iter=level->sy_relation[j].begin();
				while(iter != level->sy_relation[j].end()){
					level->sy_relation[*iter].erase(j);
					if(level->sy_relation[*iter].empty()) level->sy_relation.erase(*iter);
					iter++;
				}
				level->sy_relation.erase(j);
			}
		}	
	}
	/* end of updating */

	level->nb = cluster_index;

	return;
}

/* insertion */
void new_fully_adv_k_center_add(Fully_adv_cluster * level, unsigned int index, unsigned int cluster_index, unsigned int tot_points)
{
	double tmp;

	if(cluster_index == level->nb){
		//if current cluster is the last set and is not the unclustered point set, this center will be center of current cluster
		if (level->nb < level->k) insert_last_center(level, index);
		//if current cluster is the unclustered point set, insert directly
		else add_element_set_collection(&(level->clusters), index, level->nb);
		return;
	}else{
		if(1 == (rand() % tot_points)){
			unsigned int size;
			unsigned int helper_array[tot_points];

			add_syrelation_delete_current_center(level, cluster_index, helper_array, &size);
        		insert_last_center(level, index);

        		//Do what algorithm 1 does
        		for(int i=0; i<size; i++) fully_adv_k_center_add(level, helper_array[i]);

			return;
		}else{
			tmp = fully_adv_distance(((Geo_point *)level->array) + index,
					 ((Geo_point *)level->array) + level->centers[cluster_index]);
        		if(level->radius >= tmp){
            			add_element_set_collection(&(level->clusters), index, cluster_index);
                		level->true_rad[cluster_index] = MAX(tmp, level->true_rad[cluster_index]);
                		return;
			}else{
				tot_points -= level->clusters.sets[cluster_index].card;
                		new_fully_adv_k_center_add(level, index, cluster_index+1, tot_points);
				return;
			}
		}
	}

	return;
}

void fully_adv_k_center_add(Fully_adv_cluster * level, unsigned int index)
{
	unsigned int i;
	double tmp;
	std::vector<unsigned int> cur;
	for (i = 0; i < level->nb; i++) {
		tmp = fully_adv_distance(((Geo_point *)level->array) + index,
					 ((Geo_point *)level->array) + level->centers[i]);
		if (level->radius >= tmp) {
			add_element_set_collection(&(level->clusters), index, i);
			level->true_rad[i] = MAX(tmp, level->true_rad[i]);
			return;
		}
		if(tmp<=3*level->radius) cur.push_back(i);
	}
	add_element_set_collection(&(level->clusters), index, level->nb);
	if (level->nb < level->k) {
		/* update symmetric relation */
		for(unsigned int j=0; j<cur.size(); j++){
			level->sy_relation[level->nb].insert(cur.at(j));
			level->sy_relation[cur.at(j)].insert(level->nb);
		}
		/* end of updating */
		level->centers[level->nb] = index;
		level->true_rad[level->nb] = 0;
		level->nb++;
	}

	return;
}

/* deletion */
void fully_adv_k_center_delete(Fully_adv_cluster * level, unsigned int element_index,
			  unsigned int helper_array[])
{
	unsigned int i, size, cluster_index;
	cluster_index = get_set_index(&(level->clusters), element_index);
	remove_element_set_collection(&(level->clusters), element_index);
	if (cluster_index < level-> k && element_index == level->centers[cluster_index]) {
		add_syrelation_delete_current_center(level, cluster_index, helper_array, &size);
		for (i = 0; i < size; i++) fully_adv_k_center_add(level, helper_array[i]);
	}

	return;
}

//delete @sel_center, record the centers whose bhat value are infected and update symmetric relation
void select_update_bhat_and_syrelation(Fully_adv_cluster *level, std::unordered_map<unsigned int, std::unordered_set<unsigned int> >* sy_relation, 
				std::unordered_map<unsigned int, unsigned int>* modify_bhat, unsigned int sel_center, std::unordered_set<unsigned int>* old_centers, 
				std::vector<unsigned>* sel_clusters, 
				std::set<std::pair<unsigned int, unsigned int>, std::greater<std::pair<unsigned int, unsigned int> > >* b_hat,
				std::unordered_map<unsigned int, unsigned int>* b_hat_card)
{
	(*sel_clusters).push_back(sel_center);

	(*old_centers).erase(sel_center);
		
	unsigned int temp_bhat=(*b_hat_card)[sel_center];
	(*b_hat).erase(std::make_pair(temp_bhat, sel_center));
	(*b_hat_card).erase(sel_center);

	if((*sy_relation).find(sel_center)!=(*sy_relation).end()){
		std::unordered_set<unsigned int>::iterator neighbors =(*sy_relation)[sel_center].begin();

		while(neighbors != (*sy_relation)[sel_center].end()){
/*			if(*neighbors == sel_center){
				neighbors++;
				continue;
			}
*/
			if((*modify_bhat).find(*neighbors) != (*modify_bhat).end()){
				(*modify_bhat)[*neighbors]+=level->clusters.sets[sel_center].card;
			}else (*modify_bhat)[*neighbors] = level->clusters.sets[sel_center].card;

			(*sy_relation)[*neighbors].erase(sel_center);
			if((*sy_relation)[*neighbors].empty()) (*sy_relation).erase(*neighbors);

			neighbors++;
		}
		(*sy_relation).erase(sel_center);
	}

	return;
}

/* greedy clustering */
void fully_adv_k_center_greedy(Fully_adv_cluster levels[], unsigned int nb_instances, unsigned int k){

	/*
		@nb_clusters:the number of clusters in current lelel; @temp_bhat: temporary variable to record bhat value;
		@sel_center: temporary variable to record which center we select as new center in a round; @cur_nb: temporary variable to record how many new centers are selected
		@r: the parameter of the algorithm; @tmp: temporary variable to record the distance between two centers
		@sy_relation: a deep copy of the symmetric relation of current level
		@old_centers and @temp_old_centers: the old centers that we still need to consider for new centers in the next round
		@b_hat: store pair<b_hat value, center>, the element in @b_hat will be sorted in descending order by b_hat value;
		@b_hat_card: map center to its b_hat value, auxiliary variable to find the position of a center in @b_hat
		@modify_bhat: map center to a unsigned int which denotes how much its b_hat value should decrease in this round
		@sel_clusters: temporary variable to record which clusters corresponding to old center(s) are selected for @sel_center
	*/
	unsigned int nb_clusters, temp_bhat, sel_center, cur_nb;
	double r, tmp;
	std::unordered_map<unsigned int, std::unordered_set<unsigned int> > sy_relation;
	std::unordered_set<unsigned int> old_centers, temp_old_centers;
	std::set<std::pair<unsigned int, unsigned int>, std::greater<std::pair<unsigned int, unsigned int> > > b_hat;
	std::unordered_map<unsigned int, unsigned int> b_hat_card, modify_bhat;
	std::vector<unsigned> sel_clusters;

	std::unordered_set<unsigned int>::iterator iter;
	std::unordered_set<unsigned int>::iterator iter_old_centers;
	std::unordered_set<unsigned int>::iterator center;

	for(unsigned int instance=0; instance<nb_instances; instance++){

		levels[instance].selected_centers.clear();
		levels[instance].selected_clusters.clear();
	
		sy_relation = levels[instance].sy_relation;

		nb_clusters=levels[instance].nb;
		r=levels[instance].radius/2;

		old_centers.clear();
		temp_old_centers.clear();
		for(unsigned int i=0; i<nb_clusters; i++) old_centers.insert(i);

		b_hat.clear();
		b_hat_card.clear();

		/* compute initial b_hat and b_hat_card */
		for(unsigned int i=0; i<nb_clusters; i++){
			temp_bhat=levels[instance].clusters.sets[i].card;
			if(sy_relation.find(i)!=sy_relation.end()){
				iter=sy_relation[i].begin();
				while(iter!=sy_relation[i].end()){
					temp_bhat+=levels[instance].clusters.sets[*iter].card;
					iter++;
				}
			}
			b_hat.insert(std::make_pair(temp_bhat, i));
			b_hat_card[i]=temp_bhat;
		}

		cur_nb=0;

		/* initialize covered_points */
		levels[instance].covered_points=0;

		while(cur_nb<k && old_centers.size()>0){

			modify_bhat.clear();
			sel_clusters.clear();

			sel_center=b_hat.begin()->second;

			select_update_bhat_and_syrelation(levels+instance, &sy_relation, &modify_bhat, sel_center, &old_centers, &sel_clusters, &b_hat, &b_hat_card);

			/* update covered points */
			levels[instance].covered_points+=levels[instance].clusters.sets[sel_center].card;
		
			if(!old_centers.empty()){

				temp_old_centers = old_centers;

				iter_old_centers=old_centers.begin();
				while(iter_old_centers != old_centers.end()){
					tmp = fully_adv_distance(((Geo_point *)levels[instance].array) + levels[instance].centers[sel_center],
					 	((Geo_point *)levels[instance].array) + levels[instance].centers[*iter_old_centers]);
					if(tmp<=12*r){
						select_update_bhat_and_syrelation(levels+instance, &sy_relation, &modify_bhat, 
											*iter_old_centers, &temp_old_centers, &sel_clusters, &b_hat, &b_hat_card);
						/* update covered points */
						levels[instance].covered_points+=levels[instance].clusters.sets[*iter_old_centers].card;

					}

					iter_old_centers++;
				}
				old_centers = temp_old_centers;
				temp_old_centers.clear();
			}
			levels[instance].selected_centers.push_back(sel_center);
			levels[instance].selected_clusters[sel_center]=sel_clusters;
			
			//update b_hat value
			if(!old_centers.empty()){
				center= old_centers.begin();
				while(center != old_centers.end()){
					if(modify_bhat.find(*center) != modify_bhat.end()){

						temp_bhat=b_hat_card[*center];
						b_hat.erase(std::make_pair(temp_bhat, *center));
				
						temp_bhat-=modify_bhat[*center];
						b_hat.insert(std::make_pair(temp_bhat, *center));
						b_hat_card[*center]=temp_bhat;
					}
					center++;
				}
			}
			cur_nb++;
		}
	}
	return;
}

void fully_adv_compute_true_radius(Fully_adv_cluster * level)
{
	unsigned int cur, nb=level->selected_centers.size(), nb_centers, nb_points, temp_center;
	double max_rad, tmp;

	level->true_radius = 0;	

	/* compute true radius for each new clusters */
	for(unsigned int i=0; i<nb; i++){
		cur=level->selected_centers.at(i);
		max_rad = level->true_rad[cur]; /* let initial ture radius be the radius of the old clusters */
		nb_centers=level->selected_clusters[cur].size(); /* the number of old clusters */
		for(unsigned int j=0; j<nb_centers; j++){
			if(level->selected_clusters[cur].at(j)==cur) continue; /* if the old center is exactly the new center, skip */
			
			temp_center=level->selected_clusters[cur].at(j);

			nb_points=level->clusters.sets[temp_center].card; /* the number of points in the old clusters */
			for(unsigned int k=0; k<nb_points; k++){
				/* compute the distance between the new center and the point in the old clusters */
				tmp = fully_adv_distance(((Geo_point *)level->array) + level->centers[cur],
						((Geo_point *)level->array) + level->clusters.sets[temp_center].elements[k]);
				max_rad = MAX(max_rad, tmp);
			}
		}
		level->true_radius = MAX(level->true_radius, max_rad);
		//level->true_radius[cur]=max_rad; /* store the true radius in @trie_radius */
	}
	return;
}


Error_enum fully_adv_write_log(Fully_adv_cluster levels[], unsigned int nb_instances)
{
	unsigned int result, z;

	/* can try different epsilon and different z here */
	double epsilons[3]={0, 0.1, 0.2};
	double zratio[7]={0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};

	if (has_log()) {
		fprintf(get_log_file(), "Current Points: %u \n", NB_POINTS);
		for(int i=0; i<7; i++){
			for(int j=0; j<3; j++){
				/* z*=(1+eps)*z */
				z=(unsigned int)((1+epsilons[j])*(NB_POINTS*zratio[i]));

				result = fully_adv_get_index_smallest(levels, nb_instances, z);

				fprintf(get_log_file(), "eps: %lf zratio: %lf true radius: %lf \n", 
							epsilons[j], zratio[i], levels[result].true_radius);
			}
			fprintf(get_log_file(), "\n");
		}
		fprintf(get_log_file(), "\n");
	}
	return NO_ERROR;
}


Error_enum fully_adv_apply_one_query(Fully_adv_cluster levels[], unsigned int nb_instances,
			  struct query * query, unsigned int *helper_array, unsigned int k)
{
	unsigned int i;
	if (query->type == ADD) {
		NB_POINTS++;
		for (i = 0; i < nb_instances; i++){
			new_fully_adv_k_center_add(levels + i, query->data_index, 0, NB_POINTS);
		}
	} else {
		NB_POINTS--;
		for (i = 0; i < nb_instances; i++){
			fully_adv_k_center_delete(levels + i, query->data_index, helper_array);
		}
	}

	return NO_ERROR;
}

void
fully_adv_k_center_run(Fully_adv_cluster levels[], unsigned int nb_instances,
		       struct query_provider * queries, unsigned int helper_array[], unsigned int k)
{
	struct query query;
	
	long int initial_memory_usage = memory_usage();

	struct timeval start, end;
	double timeuse;

	srand(101);

	printf("Run starts: memory usage %8ld KB\n", initial_memory_usage);
	while (get_next_query_set(queries, &query, &(levels[0].clusters))) {

		gettimeofday(&start, NULL);
		fully_adv_apply_one_query(levels, nb_instances, &query,
					  helper_array, k);

		/* Compute true radius for each layer */
		for(int i=0; i<nb_instances; i++){
			fully_adv_compute_true_radius(levels+i);
		}
		gettimeofday(&end, NULL);

		timeuse = (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0;

		fully_adv_write_log(levels, nb_instances);
	}

	return;
}
