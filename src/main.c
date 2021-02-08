/**
Main module
**/

#include "utils.h"
#include "point.h"
#include "data_fully_adv.h"
#include "query.h"
#include "algo_fully_adv.h"
#include "algo_lp_adv.h"
#include "algo_disk_adv.h"

#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include <iostream>

char *prog_name = NULL;

typedef enum {
	FULLY_ADV_K_CENTER,
	LP_K_CENTER,
	DISK_K_CENTER,
	LAST_ALGO_TYPE
} Algo_type;

struct program_args {
	char *points_path;	/* path to points file */
	char *queries_path;	/* path for query file */
	unsigned int k;		/* the famous k */
	double epsilon;		/* approximation required */
	char log_file[120];	/* path of log_file */
	int long_log;		/* type of log (display true radius or not */
	unsigned int window_length;	/* size of sliding window */
	double d_min;		/* lower bound given */
	double d_max;		/* upper_bound given */
	Algo_type algo;		/* algo type asked */
	int parallel;		/* multithread asked by user */
	unsigned int nb_thread;	/* nb of thread asked by user */
	unsigned int cluster_size;	/* limit of cluster size specified by user */
	unsigned int t;		//failure probability
	double tau;
	unsigned int z;
};

void help(void)
{
	fprintf(stderr,
		"Linear program: %s -s [-l log_file] k eps d_min d_max data_file query_file outliers\n",
		prog_name);
	fprintf(stderr,
		"3-approx program: %s -p [-l log_file] k eps d_min d_max data_file query_file outliers\n",
		prog_name);
	fprintf(stderr,
		"Fully adversary: %s -m [-l log_file] k t d_min d_max data_file query_file tau\n",
		prog_name);
}

double n_log_n(double n)
{
	return n * log(n);
}

void init_prog_args(struct program_args *prog_args)
{
	prog_args->cluster_size = 0;
	prog_args->algo = LAST_ALGO_TYPE;
	prog_args->epsilon = -1;
	prog_args->long_log = 0;
	prog_args->window_length = 0;
	prog_args->parallel = 0;
	prog_args->log_file[0] = '\0';

	prog_args->tau = -1;
	prog_args->z=0;
}

int __parse_options(int argc, char *argv[], struct program_args *prog_args)
{
	Error_enum tmp;
	int opt;
	while ((opt = getopt(argc, argv, "hvl:tsmpn:bc:u:oe")) != -1) {
		switch (opt) {
		case 'u':
			enable_time_log(optarg);
			break;
		case 't':
			prog_args->long_log = 1;
			break;
		case 'h':
			help();
			return 1;
		case 'l':
			strcpy(prog_args->log_file, optarg);
			break;
		case 's':
			prog_args->algo = LP_K_CENTER;
			break;
		case 'p':
			prog_args->algo = DISK_K_CENTER;
			break;
		case 'm':
			prog_args->algo = FULLY_ADV_K_CENTER;
			break;
		case 'c':
			tmp = (Error_enum)strtoui_wrapper(optarg, &prog_args->cluster_size);
			if (tmp || 0 == prog_args->cluster_size) {
				fprintf(stderr,
					"Positive cluster size required for -c option");
				exit(EXIT_FAILURE);
			}
			break;
		default:
			fprintf(stderr, "unknown option\n");
			help();
			exit(EXIT_FAILURE);
		}
	}
	return 0;
}

int parse_options(int argc, char *argv[], struct program_args *prog_args)
{
	char *err_ptr = NULL;
	int next_arg = 0;
	prog_name = argv[0];
	init_prog_args(prog_args);
	if (__parse_options(argc, argv, prog_args))
		return 1;

	if (argc -optind != 7) {
		help();
		exit(EXIT_FAILURE);
	}
	if (strtoui_wrapper(argv[optind], &prog_args->k) || 0 >= prog_args->k) {
		fprintf(stderr, "positive k required\n");
		help();
		exit(EXIT_FAILURE);
	}
	next_arg++;
	if(FULLY_ADV_K_CENTER != prog_args->algo){
		if (strtod_wrapper(argv[optind + next_arg], &prog_args->epsilon)
		    || 0 >= prog_args->epsilon) {
			fprintf(stderr, "positive eps required\n");
			help();
			exit(EXIT_FAILURE);
		}
		next_arg++;
	}else{
		if (strtoui_wrapper(argv[optind + next_arg], &prog_args->t) || 0 >= prog_args->t) {
			fprintf(stderr, "positive t required\n");
			help();
			exit(EXIT_FAILURE);
		}
		next_arg++;
	}
	prog_args->d_min = strtod(argv[optind + next_arg], &err_ptr);
	if (prog_args->d_min <= 0) {
		fprintf(stderr, "positive d_min required\n");
		help();
		exit(EXIT_FAILURE);
	}
	next_arg++;
	prog_args->d_max = strtod(argv[optind + next_arg], &err_ptr);
	if (prog_args->d_max < prog_args->d_min) {
		fprintf(stderr, "d_max should be greater than d_min !\n");
		help();
		exit(EXIT_FAILURE);
	}
	next_arg++;

	prog_args->points_path = argv[optind + next_arg];
	next_arg++;

	prog_args->queries_path = argv[optind + next_arg];
	next_arg++;
	if(FULLY_ADV_K_CENTER == prog_args->algo){
		if (strtod_wrapper(argv[optind + next_arg], &prog_args->tau)
	    		|| 0 >= prog_args->tau) {
			fprintf(stderr, "positive tau required\n");
			help();
			exit(EXIT_FAILURE);
		}
		next_arg++;
	}else{
		if (strtoui_wrapper(argv[optind + next_arg], &prog_args->z) || 0 > prog_args->z) {
			fprintf(stderr, "positive z required\n");
			help();
			exit(EXIT_FAILURE);
		}
		next_arg++;
	}	
	if(FULLY_ADV_K_CENTER == prog_args->algo){
		printf("k: %d t: %d d_min: %lf d_max: %lf tau: %lf \n", prog_args->k,
	       		prog_args->t, prog_args->d_min, prog_args->d_max, prog_args->tau);}
	else 
		printf("k: %d eps: %lf d_min: %.12lf d_max: %lf z: %d \n", prog_args->k,
	       		prog_args->epsilon, prog_args->d_min, prog_args->d_max, prog_args->z);
	return 0;
}

void fully_adv_k_center(struct program_args *prog_args)
{
	Fully_adv_cluster *clusters_array;
	void *array;
	struct query_provider queries;
	unsigned int size, nb_instances;
	unsigned int *helper_array;
	fully_adv_import_points(&array, &size, prog_args->points_path);
	initialise_query_provider(&queries, prog_args->queries_path);
	if (0 == prog_args->cluster_size)
		prog_args->cluster_size = size;

	printf("before initializing memory %8ld KB\n", memory_usage());
	fully_adv_initialise_level_array(&clusters_array, prog_args->t,
					 prog_args->tau, prog_args->d_min,
					 prog_args->d_max, &nb_instances, array,
					 size, prog_args->cluster_size,
					 &helper_array);
	fully_adv_k_center_run(clusters_array, nb_instances, &queries,
			       helper_array, prog_args->k);
	printf("memory %8ld KB\n", memory_usage());
	free(array);
	fully_adv_delete_level_array(clusters_array, nb_instances,
				     helper_array);
	free_query_provider(&queries);
}

void lp_adv_k_center(struct program_args *prog_args)
{
	Lp_cluster *cluster=new Lp_cluster();

	void *array;
	struct query_provider queries;
	unsigned int size;
	fully_adv_import_points(&array, &size, prog_args->points_path);
	initialise_query_provider(&queries, prog_args->queries_path);
	if (0 == prog_args->cluster_size)
		prog_args->cluster_size = size;

	lp_initialise_level(cluster, prog_args->k, array, prog_args->z);
	lp_k_center_run(cluster, &queries, prog_args->epsilon, prog_args->d_min, prog_args->d_max);
	printf("memory %8ld KB\n", memory_usage());
	free(array);
	lp_delete_level(cluster);
	free_query_provider(&queries);
}

void disk_adv_k_center(struct program_args *prog_args)
{
	Disk_cluster *cluster=new Disk_cluster();

	void *array;
	struct query_provider queries;
	unsigned int size;
	fully_adv_import_points(&array, &size, prog_args->points_path);
	initialise_query_provider(&queries, prog_args->queries_path);
	if (0 == prog_args->cluster_size)
		prog_args->cluster_size = size;

	disk_initialise_level(cluster, prog_args->k, array, prog_args->z);

	disk_k_center_run(cluster, &queries, prog_args->epsilon, prog_args->d_min, prog_args->d_max);
	printf("memory %8ld KB\n", memory_usage());
	free(array);
	disk_delete_level(cluster);
	free_query_provider(&queries);
}

int main(int argc, char *argv[])
{
	struct program_args prog_args;
	srand48(time(NULL));
	srand((unsigned int)time(NULL));
	if (parse_options(argc, argv, &prog_args))
		return 0;
	if (prog_args.long_log)
		enable_long_log(prog_args.log_file);
	else
		enable_log(prog_args.log_file);
	switch (prog_args.algo) {
	case FULLY_ADV_K_CENTER:
		printf("Fully dynamic algorithm chosen\n");
		fully_adv_k_center(&prog_args);
		break;
	case LP_K_CENTER:
		printf("Linear program algorithm chosen\n");
		lp_adv_k_center(&prog_args);
		break;
	case DISK_K_CENTER:
		printf("Offline greedy algorithm chosen\n");
		disk_adv_k_center(&prog_args);
		break;
	default:
		fprintf(stderr, "Unknow algorithm\n");
		return EXIT_FAILURE;
	}
	disable_log();
	if (has_time_log())
		disable_time_log();
	return EXIT_SUCCESS;
}
