#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
	if (argc<3) {
		printf("Usage: %s <velocity file> <index to extract>\n", argv[0]);
		exit(1);
	}
	char* filename = argv[1];
	long index = atol(argv[2]);
	FILE* fp_in = fopen(filename, "rb");
	fseek(fp_in, index*4, SEEK_SET);
	float value;
	fread(&value, 1, sizeof(float), fp_in);
	fclose(fp_in);
	printf("Value in file %s at index %ld is %f.\n", filename, index, value);
	return 0;
}
