#include "include.h"

int main(int argc, char** argv) {
	if (argc<5) {
		printf("Usage: %s <sgt header> <x> <y> <z>\n", argv[0]);
		exit(1);
	}

	char* sgthead_in = argv[1];
	int x = atoi(argv[2]);
	int y = atoi(argv[3]);
	int z = atoi(argv[4]);

	FILE* fp_in = fopen(sgthead_in, "rb");
	struct sgtmaster sgtmast;
	fread(&sgtmast, sizeof(struct sgtmaster), 1, fp_in);
	int i;
	struct sgtindex* sgtindx = malloc(sizeof(struct sgtindex)*sgtmast.globnp);
	fread(sgtindx, sizeof(struct sgtindex), sgtmast.globnp, fp_in);
	fclose(fp_in);
	for (i=0; i<sgtmast.globnp; i++) {
		if (sgtindx[i].xsgt==x && sgtindx[i].ysgt==y && sgtindx[i].zsgt==z) {
			printf("Grid point (%d, %d, %d) found at SGT index %d.\n", x, y, z, i);
			free(sgtindx);
			return 0;
		}
	}
	printf("Grid point (%d, %d, %d) not found.\n", x, y, z);
	free(sgtindx);
	return 1;
}
