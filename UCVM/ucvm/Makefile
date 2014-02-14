

TARGET = ucvm dem

all: $(TARGET)
.PHONY: ucvm dem



############################################
# Targets
############################################

ucvm:
	cd ./src/ucvm; \
	make all;


dem:
	cd ./src/dem; \
	make all;


############################################
# Special build targets
############################################


clean:
	cd ./src/ucvm; \
	make clean; \
	cd ../dem; \
	make clean;

