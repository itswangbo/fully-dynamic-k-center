/**
The module contains the code concerning the fully adversary algorithm on GPS point
**/
#include "utils.h"
#include "query.h"
#include "data_fully_adv.h"
#include "algo_lp_adv.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <float.h>

#include <algorithm>
#include <set>
#include <map>
#include <vector>
#include <unordered_set>

#include "gurobi_c++.h"  /* c++ library for solving linear program */
/* creating GRB environment */
static GRBEnv* env;

static const double eps=1e-6;

void
lp_initialise_level(Lp_cluster * level, unsigned int k,
			   void * array, unsigned int z)
{
	env=new GRBEnv();

	/* turn off printing */
	env->set(GRB_IntParam_OutputFlag, 0);

	level->k = k;
	level->z = z;

	level->array = array;
}

void lp_delete_level(Lp_cluster * level)
{
	level->selected_centers.clear();
	level->points.clear();

	level->array = NULL;

	delete level;

	delete env;
}

/* insertion */
void lp_k_center_add(Lp_cluster * level, unsigned int index)
{
	level->points.insert(index);
}

/* deletion */
void lp_k_center_delete(Lp_cluster * level, unsigned int index)
{
	level->points.erase(index);
}

/* linear program solver */
int lp_solver(Lp_cluster *level, std::vector<double>* x_value, std::vector<double>* y_value, std::unordered_map<unsigned int, unsigned int>* old_new,
		std::unordered_map<unsigned int, std::vector<unsigned int> > sy_relation)
{
	unsigned int count=0, cnt;
	std::unordered_map<unsigned int, std::vector<unsigned int> >::iterator iter=sy_relation.begin();
	while(iter != sy_relation.end()){
		old_new->insert(std::make_pair(iter->first, count));
		count++;
		iter++;
	}

	/* linear program */
	GRBModel* model=new GRBModel(*env);

	int cur_points=count;

	int tot_vars=cur_points*cur_points+cur_points;
	double* lb=new double[tot_vars];
	double* ub=new double[tot_vars];
	double* obj=new double[tot_vars];
	char* type=new char[tot_vars];

	/* initialize tot_vars variables */
	for(int i=0; i<tot_vars; i++){
		lb[i]=0.0;		/* lower bound */
		ub[i]=1.0;		/* upper bound */
		obj[i]=0.0;		/* objective coefficient for new variable */
		type[i]=GRB_CONTINUOUS;	/* variable type */
	}
	/* the name of each variable is set to be NULL */
	GRBVar* vars=model->addVars(lb, ub, obj, type, NULL, tot_vars);

	GRBLinExpr rsum1=0;

	iter = sy_relation.begin();
   	while(iter != sy_relation.end()){

        	cnt = iter->second.size();

        	GRBLinExpr rsum2=0;

        	for(int pos=0; pos<cnt; pos++){
            		/* constraint 1 */
            		rsum1+=vars[(*old_new)[iter->second.at(pos)]*(cur_points) + (*old_new)[iter->first]];

            		/* constraint 2 */
            		rsum2+=vars[(*old_new)[iter->second.at(pos)]*(cur_points) + (*old_new)[iter->first]];

            		/* constraint 3 */
           	 	GRBLinExpr rsum3=0;
            		rsum3+=vars[(*old_new)[iter->second.at(pos)]*(cur_points) + (*old_new)[iter->first]];
            		rsum3-=vars[cur_points*cur_points+(*old_new)[iter->second.at(pos)]];
            		model->addConstr(rsum3 <= 0);	/* add constraint 3: x(i,j) - y_i <= 0 */
        	}

		model->addConstr(rsum2 <= 1);	/* add constraint 2: sigma(i~B_j) x(i,j) <= 1 */

        	iter++;
	}

	/* at least t points are covered */
	model->addConstr(rsum1 >= (cur_points)-(level->z));	/* add constraint 1: sigma(j)sigma(i~B_j) x(i,j) >= t */

	/* constraint 6 */
	GRBLinExpr rsum6=0;
	for(unsigned int it=0; it<cur_points; it++){
		rsum6+=vars[cur_points*cur_points+it];
	}

	model->addConstr(rsum6 <= level->k);	/* add constraint 6: sigma(i) y_i <= k */

	/* objective function: 0*x(0,0) */
	GRBLinExpr obj_f=0;
	model->setObjective(obj_f, GRB_MAXIMIZE);

	try{
		/* optimize the model */
		model->optimize();
	}catch(GRBException e){
		std::cout<<"Error code: "<<e.getErrorCode() << std::endl;
		std::cout<<e.getMessage()<<std::endl;
	}catch(...){
		std::cout<<"Exception during optimization." << std::endl;
	}

	/* 2: optimal; 3: infeasible */
	int stas = model->get(GRB_IntAttr_Status);

	/* if optimal solution, get the final value of variables */
	if(stas == GRB_OPTIMAL){
        	int cur_cnt, bari=cur_points*cur_points;
        	GRBVar* vars_val=model->getVars();
    		for(cur_cnt=0; cur_cnt<bari; cur_cnt++){
        		x_value->at(cur_cnt) = vars[cur_cnt].get(GRB_DoubleAttr_X);
    		}
    		for(; cur_cnt<tot_vars; cur_cnt++){
        		y_value->at(cur_cnt-bari)=vars[cur_cnt].get(GRB_DoubleAttr_X);
    		}
    		delete vars_val;
	}

	delete [] lb;
	delete [] ub;
	delete [] obj;
	delete [] type;

	delete vars;
	delete model;

	return stas;
}

/* modified linear program solver */
int modify_lp_solver(Lp_cluster *level, std::vector<double>* x_value, std::vector<double>* y_value, std::unordered_map<unsigned int, unsigned int>* old_new,
		std::unordered_map<unsigned int, std::vector<unsigned int> > sy_relation)
{
	unsigned int count=0, cnt;
	std::unordered_map<unsigned int, std::vector<unsigned int> >::iterator iter=sy_relation.begin();
	while(iter != sy_relation.end()){
		old_new->insert(std::make_pair(iter->first, count));
		count++;
		iter++;
	}

	/* linear program */
	GRBModel* model=new GRBModel(*env);

	int cur_points=count;

	int tot_vars=cur_points*cur_points+cur_points;
	double* lb=new double[tot_vars];
	double* ub=new double[tot_vars];
	double* obj=new double[tot_vars];
	char* type=new char[tot_vars];

	/* initialize tot_vars variables */
	for(int i=0; i<tot_vars; i++){
		lb[i]=0.0;		/* lower bound */
		ub[i]=1.0;		/* upper bound */
		obj[i]=0.0;		/* objective coefficient for new variable */
		type[i]=GRB_CONTINUOUS;	/* variable type */
	}
	/* the name of each variable is set to be NULL */
	GRBVar* vars=model->addVars(lb, ub, obj, type, NULL, tot_vars);

	GRBLinExpr rsum1=0;

	iter = sy_relation.begin();
   	while(iter != sy_relation.end()){

        cnt = iter->second.size();

        GRBLinExpr rsum2=0;

        for(int pos=0; pos<cnt; pos++){
            /* constraint 1 */
            rsum1+=vars[(*old_new)[iter->second.at(pos)]*(cur_points) + (*old_new)[iter->first]];

            /* constraint 2 */
            rsum2+=vars[(*old_new)[iter->second.at(pos)]*(cur_points) + (*old_new)[iter->first]];

            /* constraint 3 */
            GRBLinExpr rsum3=0;
            rsum3+=vars[(*old_new)[iter->second.at(pos)]*cur_points + (*old_new)[iter->first]];
            rsum3-=vars[cur_points*cur_points+(*old_new)[iter->second.at(pos)]];
            model->addConstr(rsum3 <= 0);	/* add constraint 3: x(i,j) - y_i <= 0 */
        }

		model->addConstr(rsum2 <= 1);	/* add constraint 2: sigma(i~B_j) x(i,j) <= 1 */

        iter++;
	}

	/* at least t points are covered */
	//model->addConstr(rsum1 >= (cur_points)-(level->z));	/* add constraint 1: sigma(j)sigma(i~B_j) x(i,j) >= t */

	/* constraint 6 */
	GRBLinExpr rsum6=0;
	for(unsigned int it=0; it<cur_points; it++){
		rsum6+=vars[cur_points*cur_points+it];
	}

	model->addConstr(rsum6 <= level->k);	/* add constraint 6: sigma(i) y_i <= k */

	/* objective function: 0*x(0,0) */
	model->setObjective(rsum1, GRB_MAXIMIZE);

	try{
		/* optimize the model */
		model->optimize();
	}catch(GRBException e){
		std::cout<<"Error code: "<<e.getErrorCode() << std::endl;
		std::cout<<e.getMessage()<<std::endl;
	}catch(...){
		std::cout<<"Exception during optimization." << std::endl;
	}

	/* 2: optimal; 3: infeasible */
	int stas = model->get(GRB_IntAttr_Status);

	/* if optimal solution, get the final value of variables */
	if(stas == GRB_OPTIMAL){
		double objVal=model->get(GRB_DoubleAttr_ObjVal);

		if(objVal<(cur_points-level->z-eps)){
			stas = GRB_INFEASIBLE;
		}else{
	        	int cur_cnt, bari=cur_points*cur_points;
        		GRBVar* vars_val=model->getVars();
    			for(cur_cnt=0; cur_cnt<bari; cur_cnt++){
        			x_value->at(cur_cnt) = vars[cur_cnt].get(GRB_DoubleAttr_X);
    			}
    			for(; cur_cnt<tot_vars; cur_cnt++){
        			y_value->at(cur_cnt-bari)=vars[cur_cnt].get(GRB_DoubleAttr_X);
    			}
    			delete vars_val;
		}
	}

	delete [] lb;
	delete [] ub;
	delete [] obj;
	delete [] type;

	delete vars;
	delete model;

	return stas;
}

/* filtering */
void r_filtering(Lp_cluster *level, std::vector<double> x_value, std::vector<double> y_value, std::unordered_map<unsigned int, unsigned int> old_new,
	std::set<std::pair<unsigned int, unsigned int>, std::greater<std::pair<unsigned int, unsigned int> > >* c_hat,
	std::unordered_map<unsigned int, std::vector<unsigned int> > sy_relation)
{
	std::set<std::pair<double, unsigned int>, std::greater<std::pair<double, unsigned int> > > F; /* s_j value, index j */
	std::unordered_map<unsigned int, std::vector<unsigned int> >::iterator iter=sy_relation.begin();
	double tmp, ttmp;
	int cur_points=level->points.size(), cnt;

	while(iter != sy_relation.end()){
		tmp=0;
		cnt=iter->second.size();
		for(int pos=0; pos<cnt; pos++){
			ttmp=x_value.at(old_new[iter->second.at(pos)]*cur_points + old_new[iter->first]); /* s_j: x_ij, i belongs to B_j */
            		tmp+=ttmp;
		}
		F.insert(std::make_pair(tmp, iter->first));
		iter++;
	}

	unsigned int *flag = new unsigned int[cur_points](); /* to determine which clusters are already marked, all zero */
	unsigned int selected_cluster, point_in_cluster, intersect_cluster, cur_c; /* how many clusters are selected */
	int size, ssize;
	while(!F.empty()){
		while(!F.empty() && flag[ old_new[F.begin()->second] ] == 1) F.erase(F.begin());
		if(F.empty()) break;
		selected_cluster=F.begin()->second;
		size = sy_relation[selected_cluster].size();
		cur_c=0;
		/* select current cluster */
		for(int pos=0; pos<size; pos++){
			point_in_cluster=sy_relation[selected_cluster].at(pos);
			ssize = sy_relation[point_in_cluster].size();
			/* find all unmarked clusters intersect with selected cluster */
			for(int ppos=0; ppos<ssize; ppos++){
				/* if this cluster which intersect with selected cluster is unmarked */
				intersect_cluster=sy_relation[point_in_cluster].at(ppos);
				if(flag[old_new[intersect_cluster]]==0){
					flag[old_new[intersect_cluster]]=1;
					cur_c++;
				}
			}
		}
		c_hat->insert(std::make_pair(cur_c, selected_cluster));
		F.erase(F.begin());
	}
	delete [] flag;
}

/* k center */
void r_kcenter_round(Lp_cluster *level, std::vector<double> x_value, std::vector<double> y_value, std::unordered_map<unsigned int, unsigned int> old_new,
			std::unordered_map<unsigned int, std::vector<unsigned int> > sy_relation)
{
	std::set<std::pair<unsigned int, unsigned int>, std::greater<std::pair<unsigned int, unsigned int> > > c_hat;

	r_filtering(level, x_value, y_value, old_new, &c_hat, sy_relation);

	unsigned int cur_k=0;
	level->selected_centers.clear();
	while(cur_k<=level->k && !c_hat.empty()){
		level->selected_centers.push_back(c_hat.begin()->second);
		c_hat.erase(c_hat.begin());
		cur_k++;
	}
}

/* sy_relation: pair<j, B_j> d(i,j)<=radius */
void construct_sy_relation(Lp_cluster *level, double radius, std::unordered_map<unsigned int, std::vector<unsigned int> > *sy_relation)
{
	double tmp;
	std::unordered_set<unsigned int>::iterator iter=level->points.begin();
	std::vector<unsigned int> inserted_points;
	std::vector<unsigned int> cur;

	sy_relation->clear();

	while(iter != level->points.end()){
		for(int i=0; i<inserted_points.size(); i++){
			tmp=fully_adv_distance(((Geo_point *)level->array) + inserted_points.at(i),
				((Geo_point *)level->array) + *iter);
			if(tmp <= radius){
				cur.push_back(inserted_points.at(i));
				(*sy_relation)[inserted_points.at(i)].push_back(*iter);
			}
		}
		cur.push_back(*iter);
		sy_relation->insert(std::make_pair(*iter, cur));
		cur.clear();

		inserted_points.push_back(*iter);
		iter++;
	}
}

/* select k centers: Binary search r */
void lp_k_center_solver(Lp_cluster *level, double eps, double dmin, double dmax)
{
	/* coupute x_value and y_value by solving linear program */
	int cur_points=level->points.size(), res;

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

	std::vector<double> x_value(cur_points*cur_points);		/* at(i*cur_points+j): the probability that i and j are connected */
	std::vector<double> y_value(cur_points);				/* the probability that i is picked */
	std::unordered_map<unsigned int, unsigned int> old_new;		/* map old index to the new index {0, 1, 2, 3, ...} */
	std::unordered_map<unsigned int, std::vector<unsigned int> > sy_relation; /* pair<j, B_j> d(i,j)<=radius */

	//double radius = pow( (1+eps), ceil(log(dmin)/log(1+eps)) );
	double radius = dmin;
	unsigned int lef=0, mid, rig=(unsigned int)(1 + ceil(log(dmax / dmin) / log(1 + eps)));
//int num_lp=0;
	while(lef<rig){
		mid=lef+(rig-lef)/2;
		level->radius = radius*pow((1+eps),mid);

		construct_sy_relation(level, level->radius, &sy_relation);
//num_lp++;
        	res = modify_lp_solver(level, &x_value, &y_value, &old_new, sy_relation);
		if(res == GRB_OPTIMAL){
			rig=mid;
		}else{
			lef=mid+1;
		}
	}

	if(res != GRB_OPTIMAL){
//num_lp++;
        	level->radius = radius*pow((1+eps),lef);
        	construct_sy_relation(level, level->radius, &sy_relation);
        	modify_lp_solver(level, &x_value, &y_value, &old_new, sy_relation);
	}

//std::cout<<"number of lp: "<<num_lp<<std::endl;
	r_kcenter_round(level, x_value, y_value, old_new, sy_relation);
}
	
/* compute true_radius */
void lp_compute_true_radius(Lp_cluster * level)
{
    unsigned int cur_points = level->points.size(), outliers=level->z, num_centers=level->k;
    if(cur_points-outliers-num_centers <= 0){
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

Error_enum fully_adv_write_log(Lp_cluster* level)
{
	if (has_log()) {
		fprintf(get_log_file(), "true radius: %lf\n", level->true_radius);
	}
	return NO_ERROR;
}

Error_enum lp_apply_one_query(Lp_cluster* level, struct query * query, unsigned int &N)
{
	if (query->type == ADD) {
		/* printf("a %u\n", query->data_index); */
		N++;
		lp_k_center_add(level, query->data_index);
	} else {
		N--;
		lp_k_center_delete(level, query->data_index);
	}

	return NO_ERROR;
	//return fully_adv_write_log(levels, nb_instances, nb_points, query, z);
}

void lp_k_center_run(Lp_cluster * level, struct query_provider * queries, double epsilon, double dmin, double dmax)
{
	struct query query;

	unsigned int count = 0;
	unsigned int M = 0;
	long int initial_memory_usage = memory_usage();

	struct timeval start, end;
	double timeuse;

	printf("Run starts: memory usage %8ld KB\n", initial_memory_usage);
	while (get_next_query_set_lp(queries, &query, level->points)) {
		lp_apply_one_query(level, &query, count);
	}

	gettimeofday(&start, NULL);
	lp_k_center_solver(level, epsilon, dmin, dmax);
	lp_compute_true_radius(level);
	gettimeofday(&end, NULL);

	timeuse = (end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)/1000000.0;
	
	fully_adv_write_log(level);
}
