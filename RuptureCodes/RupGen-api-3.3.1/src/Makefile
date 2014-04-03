GETPAR = /work/00940/tera3d/CyberShake/software/Getpar/getpar
GETPARLIB = ${GETPAR}/lib
GETPARINCL = ${GETPAR}/include

GENSLIP = GenRandV3.3.1

# General compiler/linker flags
CC=gcc
FC=gfortran
#CC=icc
#FC=ifort
#CFLAGS = -I${GETPARINCL} -Wall \
#	-D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE \
#	-D_FILE_OFFSET_BITS=64 -DUSE_IOBUF_MACROS -p -O3 -xhost -vec-report=3
CFLAGS = -I${GETPARINCL} -Wall \
       -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE \
       -D_FILE_OFFSET_BITS=64 -DUSE_IOBUF_MACROS -g
LIBS = -L${GETPARLIB} -lget -lgfortran -lm
LDFLAGS = -static

# Uncomment for memcached support
MEMC = /work/00940/tera3d/CyberShake/utils/libmemcached
CFLAGS := $(CFLAGS) -I${MEMC}/include -D_USE_MEMCACHED
LIBS := ${LIBS} -L${MEMC}/lib -lmemcached

TARGET = librupgen.a

all: $(TARGET)

############################################
# Executables
############################################

librupgen.a: rupgen_api.o ${GENSLIP}/genslip-v3.3.1.o ${GENSLIP}/fourg.o ${GENSLIP}/slip.o ${GENSLIP}/ruptime.o ${GENSLIP}/iofunc.o ${GENSLIP}/srf_subs.o ${GENSLIP}/misc.o
	cd ${GENSLIP} && $(MAKE)
	$(AR) rcs $@ $^
	cp librupgen.a ../lib/
	cp rupgen_api.h ../include/
	cp structure.h ../include/

genslip: genslip.o iofunc.o srf_subs.o
	$(CC) -o $@ $^ $(LIBS)

############################################
# Special build targets
############################################

clean:
	rm -f core *.o *~ $(TARGET) ${GENSLIP}/*.o
