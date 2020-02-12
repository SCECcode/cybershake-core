TARGETS = get_num_rvs generate_hypos

CC = gcc
CFLAGS = -O3

MEMC = ../../../../utils/libmemcached_1.0.18
INC = -I../include/ -I/../../../../software/Getpar/getpar/include -I${OLCF_FFTW_ROOT}/include -I${MEMC}/include
LIBS = -L../lib -lrupgen -L../../../../software/Getpar/getpar/lib -lget -L${OLCF_FFTW_ROOT}/lib -lfftw3f -lm -L${MEMC}/lib -lmemcached

all: $(TARGETS)

clean:
	rm -rf *.o $(TARGETS)

get_num_rvs: get_num_rvs.c
	$(CC) $(CFLAGS) -o get_num_rvs get_num_rvs.c $(INC) $(LIBS)

generate_hypos: generate_hypos.c
	$(CC) $(CFLAGS) -o generate_hypos generate_hypos.c $(INC) $(LIBS)

