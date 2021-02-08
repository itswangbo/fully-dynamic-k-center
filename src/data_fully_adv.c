/**
The module contains the code about the IO operation related to the fully adversary algorithm on GPS point
**/
#define _POSIX_C_SOURCE 200112L
#include "point.h"
#include "data_fully_adv.h"
#include "utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define INIT_ARRAY_SIZE 100000

static Error_enum fully_adv_read_point(char *string, Geo_point * p)
{
	char *tmp = strtok(string, "\t");
	if (!tmp)
		return FILE_FORMAT_ERROR;
	tmp = strtok(NULL, " \t\n");
	if (strtod_wrapper(tmp, &(p->longitude)))
		return FILE_FORMAT_ERROR;
	tmp = strtok(NULL, " \t\n");
	if (strtod_wrapper(tmp, &(p->latitude)))
		return FILE_FORMAT_ERROR;
	return NO_ERROR;
}

Error_enum
fully_adv_import_points(void ** point_array,
			unsigned int *nb_element, char *path)
{
	char buffer[BUFSIZ];
	unsigned line = 1;
	size_t current = 0, max_array = INIT_ARRAY_SIZE;
	FILE *f = fopen_wrapper(path, "r");
	Error_enum tmp;
	*point_array =calloc_wrapper(INIT_ARRAY_SIZE, sizeof(Geo_point));
	while (fgets(buffer, BUFSIZ, f)) {
		if (current == max_array)
			*point_array = realloc_wrapper(*point_array,&max_array,
								sizeof(Geo_point));
		if ((tmp =
		     fully_adv_read_point(buffer, ((Geo_point *)*point_array) + current))) {
			fprintf(stderr,
				"Wrong point format, incident occured line %d\n",
				line);
			return tmp;
		}
		line++;
		current++;
	}
	*nb_element = (unsigned)current;
	fclose(f);
	return NO_ERROR;
}

void fully_adv_print_points(void * parray, unsigned int nb_elements)
{
	unsigned int i;
	Geo_point *array=(Geo_point *)parray;
	for (i = 0; i < nb_elements; i++)
		printf("%u %lf %lf\n", i, array[i].longitude,
		       array[i].latitude);
}

double fully_adv_distance(void * x, void * y)
{
	return euclidean_distance((Geo_point *)x, (Geo_point *)y);
}
