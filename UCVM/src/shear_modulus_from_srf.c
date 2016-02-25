#include <stdio.h>
#include <stdlib.h>

#include "rupgen_api.h"
#include "math.h"

struct coord_entry {
	float lon;
	float lat;
	int grid_x;
	int grid_y;
};

const int RAD = 6371;

int main(int argc, char** argv) {
	if (argc<9) {
		printf("Usage: %s <SRF file> <model coords> <velocity mesh> <nx> <ny> <nz> <grid spacing, km> <output file>\n", argv[0]);
		exit(1);
	}

	char* srf_filename = argv[1];
	char* model_coords = argv[2];
	char* velocity_mesh = argv[3];
	int nx = atoi(argv[4]);
	int ny = atoi(argv[5]);
	int nz = atoi(argv[6]);
	float grid_spacing = atof(argv[7]);
	char* output_file = argv[8];

	int i,j,k,m;

	struct standrupformat srf;
	_read_srf(&srf, srf_filename, 0);

	//Read model coords
	printf("Reading model coords.\n");
	FILE* fp_in = fopen(model_coords, "r");
	struct coord_entry* grid_points = malloc(nx*ny*sizeof(struct coord_entry));
	for (i=0; i<nx*ny; i++) {
		fscanf(fp_in, " %f %f %d %d\n", &grid_points[i].lon, &grid_points[i].lat, &grid_points[i].grid_x, &grid_points[i].grid_y);
	}
	fclose(fp_in);
	
	//Open velocity file
	fp_in = fopen(velocity_mesh, "rb");

	//Open output
	FILE* fp_out = fopen(output_file, "w");

	int num_close_points = 0;
	int closest_indices[4];
	float closest_dists[4];
	for (i=0; i<4; i++) {
		closest_indices[i] = -1;
		closest_dists[i] = 1.0e20;
	}

	//Determine shear modulus = rho*vs^2
	for (i=0; i<srf.srf_apnts.np; i++) {
		if (i%100==0) {
			printf("SRF point %d of %d.\n", i, srf.srf_apnts.np);
		}
		float srf_lon = srf.srf_apnts.apntvals[i].lon;
		float srf_lat = srf.srf_apnts.apntvals[i].lat;
		float srf_dep = srf.srf_apnts.apntvals[i].dep;
		//Find depth just above (nearest grid point)
		int top_dep_index = ((int)(srf.srf_apnts.apntvals[i].dep/grid_spacing));
		//Calculate distance to all points, keep the closest
		for (j=0; j<nx*ny; j++) {
			struct coord_entry entry = grid_points[j];
			//Ignore any points which are more than .1 degree away
			if (fabs(entry.lon-srf_lon)>0.05) {
				continue;
			}
			if (fabs(entry.lat-srf_lat)>0.05) {
				continue;
			}
			float lat_dist = (entry.lat-srf_lat)*3.1416*RAD/180.0;
			float med_lat = (entry.lat-srf_lat)/2.0;
			float lon_dist = (entry.lon-srf_lon)*3.1416*RAD/180.0*cos(med_lat*3.1416/180.0);
			float dist = sqrt(lat_dist*lat_dist+lon_dist*lon_dist);
			//Continue if the dist is greater than 2 grid points, not one of the close ones
			if (dist>0.8) {
				continue;
			}
			for (k=0; k<4; k++) {
				if (dist<closest_dists[k]) {
					for (m=2; m>=k; m--) {
						closest_dists[m+1] = closest_dists[m];
						closest_indices[m+1] = closest_indices[m];
					}
					closest_dists[k] = dist;
					closest_indices[k] = j;
				}
				break;
			}
		}
		//Get parameters for up to 8 points, average
		float tot = 0.0;
		float dist_tot = 0.0;
		for (j=0; j<4; j++) {
			float mesh_data[3];
			float shear_mod;
			if (closest_indices[j]!=-1) {
				num_close_points++;
				//Retrieve from mesh
				struct coord_entry entry = grid_points[closest_indices[j]];
				long offset = 3*sizeof(float)*(top_dep_index*nx*ny + entry.grid_x*ny + entry.grid_y);
				fseek(fp_in, offset, SEEK_SET);
				fread(mesh_data, sizeof(float), 3, fp_in);
				shear_mod = mesh_data[2]*mesh_data[1]*mesh_data[1];
				//include depth in dist calc
				float full_dist = sqrt(closest_dists[j]*closest_dists[j] + (top_dep_index*grid_spacing - srf_dep)*(top_dep_index*grid_spacing - srf_dep));
				tot += 1.0/full_dist*shear_mod;
				dist_tot += 1.0/full_dist;
				//If there's a deeper point available, do it again
				if (top_dep_index+1 < nz) {
					num_close_points++;
					offset = 3*sizeof(float)*((top_dep_index+1)*nx*ny + entry.grid_x*ny + entry.grid_y);
	                                fseek(fp_in, offset, SEEK_SET);
	                                fread(mesh_data, sizeof(float), 3, fp_in);
	                                shear_mod = mesh_data[2]*mesh_data[1]*mesh_data[1];
					full_dist = sqrt(closest_dists[j]*closest_dists[j] + ((top_dep_index+1)*grid_spacing - srf_dep)*((top_dep_index+1)*grid_spacing - srf_dep));
	                                tot += full_dist*shear_mod;
					dist_tot += full_dist;
				}
			} else {
				//done with points
				break;
			}
		}
		fprintf(fp_out, "%f %f %f %f\n", srf_lon, srf_lat, srf_dep, (tot/dist_tot));
	}
	fflush(fp_out);
	fclose(fp_out);
	fclose(fp_in);
}
