TARGETS = generate_srfs get_num_rvs get_hypo

CC = gcc
CFLAGS = -O3

MEMC = ../../../../utils/libmemcached_1.0.18
INC = -I../include/ -I/../../../../software/Getpar/getpar/include -I${TACC_FFTW3_INC} -I${MEMC}/include
LIBS = -L../lib -lrupgen -L../../../../software/Getpar/getpar/lib -lget -L${TACC_FFTW3_LIB} -lfftw3f -lm -L${MEMC}/lib -lmemcached

all: $(TARGETS)

clean:
	rm -rf *.o $(TARGETS)

generate_srfs: generate_srfs.c
	$(CC) $(CFLAGS) -o generate_srfs generate_srfs.c $(INC) $(LIBS)

get_num_rvs: get_num_rvs.c
	$(CC) $(CFLAGS) -o get_num_rvs get_num_rvs.c $(INC) $(LIBS)

get_hypo: get_hypo.c
	$(CC) $(CFLAGS) -o get_hypo get_hypo.c $(INC) $(LIBS)
