GETPAR = ../../Getpar/getpar
GETPARLIB = ${GETPAR}/lib
GETPARINCL = ${GETPAR}/include

# General compiler/linker flags
CC=gcc
FC=gfortran
#CC=icc
#FC=ifort
CFLAGS = -I${GETPARINCL} -Wall -g \
	-D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE \
	-D_FILE_OFFSET_BITS=64 -DUSE_IOBUF_MACROS -O3
LIBS = -L${GETPARLIB} -lget -lgfortran -lm
LDFLAGS = -static

# Uncomment for memcached support
#MEMC = /home/rcf-104/CyberShake2007/memcached/libmemcached-0.44
#CFLAGS := $(CFLAGS) -I${MEMC}/include -D_USE_MEMCACHED
#LIBS = ${LIBS} -L${MEMC}/lib -lmemcached

TARGET = librupgen.a

all: $(TARGET)

############################################
# Executables
############################################

librupgen.a: rupgen_api.o genslip.o fourg.o slip.o ruptime.o iofunc.o srf_subs.o misc.o
	$(AR) rcs $@ $^


genslip: genslip.o iofunc.o srf_subs.o
	$(CC) -o $@ $^ $(LIBS)

############################################
# Special build targets
############################################

clean:
	rm -f core *.o *~ $(TARGET)
