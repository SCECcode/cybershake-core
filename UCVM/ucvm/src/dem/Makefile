
# General compiler/linker flags
CFLAGS = -Wall -O3 -std=c99  \
	-D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE \
	-D_FILE_OFFSET_BITS=64
LDFLAGS = -static


LDFLAGS += -lm
#-lpthread

# Compilers
CC=gcc
FC=gfortran
LD=gcc


TARGET = dem_query

all: $(TARGET)


############################################
# Executables
############################################

dem_query: dem_query.o dem.o dem_config.o
	$(LD) -o $@ $^ $(LDFLAGS)
	cp dem_query ../../bin


############################################
# Special build targets
############################################

clean:
	rm -f core *.o *~ $(TARGET) ../../bin/dem_query

