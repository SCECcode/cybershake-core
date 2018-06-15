#include <stdio.h>
#include <stdlib.h>

#include "ucvm.h"

/* Code to query UCVM for elevation of surface points and write to file */

struct {
  float lon;
  float lat;
  float elev;
} point;


void parse_model_coords(char* model_coords_file, point** points, int num_pts) {
	FILE* fp_in = fopen(model_coords_file, "r");
	int i, x, y;
	for (i=0; i<num_pts; i++) {
		fscanf(fp_in, "%f %f %d %d\n", &((*points)[i].lon), &((*points)[i].lat), &x, &y);
	}
	fclose(fp_in);
}

void query_points(point* points, int num_pts) {
	
}

int main(int argc, char** argv) {
	if (argc<5) {
		printf("Usage: %s <num_surface_pts> <model_coords> <format> <output file>\n", argv[0]);
		exit(1);
	}
	int num_surf_pts = int(argv[1]);
	char* model_coords_file = argv[2];
	char* format = argv[3];
	char* output_file = argv[4];

	point* points = malloc(sizeof(point)*num_surf_pts);
	parse_model_coords(model_coords_file, &points, num_surf_pts);

	query_points(points, num_surf_pts);
}
