GETPAR = /work/00940/tera3d/CyberShake/software/Getpar/getpar
GETPARLIB = ${GETPAR}/lib
GETPARINCL = ${GETPAR}/include

GENSLIP = GenRandV3.3.1

GENSLIP_OBJS = ${GENSLIP}/ruptime.o ${GENSLIP}/iofunc.o ${GENSLIP}/srf_subs.o ${GENSLIP}/misc.o

# General compiler/linker flags
#CC=gcc
#FC=gfortran
CC=icc
FC=ifort
CFLAGS = -I${GETPARINCL} -Wall \
	-D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE \
	-D_FILE_OFFSET_BITS=64 -DUSE_IOBUF_MACROS -O3 -xhost -g
#CFLAGS = -I${GETPARINCL} -Wall \
#       -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE \
#       -D_FILE_OFFSET_BITS=64 -DUSE_IOBUF_MACROS -O3 -p
LIBS = -L${GETPARLIB} -lget -lgfortran -lm
LDFLAGS = -static

# Uncomment for memcached support
MEMC = /work/00940/tera3d/CyberShake/utils/libmemcached
CFLAGS := $(CFLAGS) -I${MEMC}/include -D_USE_MEMCACHED
LIBS := ${LIBS} -L${MEMC}/lib -lmemcached

TARGETS = librupgen.a librupgen_fftw.a

all: $(TARGETS)


############################################
# Executables
############################################

#librupgen.a: rupgen_api.o ${GENSLIP}/genslip-v3.3.1.o ${GENSLIP}/fourg.o ${GENSLIP}/slip.o ${GENSLIP}/ruptime.o ${GENSLIP}/iofunc.o ${GENSLIP}/srf_subs.o ${GENSLIP}/misc.o
librupgen.a: rupgen_api.o
	cd ${GENSLIP} && $(MAKE)
	$(AR) rcs $@ $^ ${GENSLIP}/genslip-v3.3.1.o ${GENSLIP}/slip.o ${GENSLIP}/fourg.o ${GENSLIP_OBJS}
	cp librupgen.a ../lib/
	cp rupgen_api.h ../include/
	cp structure.h ../include/

#librupgen_fftw.a: rupgen_api.o ${GENSLIP}/genslip-v3.3.1.fftw.o ${GENSLIP}/slip.o ${GENSLIP}/ruptime.o ${GENSLIP}/iofunc.o ${GENSLIP}/srf_subs.o ${GENSLIP}/misc.o
librupgen_fftw.a: rupgen_api.o
	cd ${GENSLIP} && $(MAKE)
	$(AR) rcs $@ $^ ${GENSLIP}/genslip-v3.3.1.fftw.o ${GENSLIP}/slip.fftw.o ${GENSLIP_OBJS}
	cp librupgen_fftw.a ../lib/
	cp rupgen_api.h ../include/
	cp structure.h ../include/

genslip: genslip.o iofunc.o srf_subs.o
	$(CC) -o $@ $^ $(LIBS)

############################################
# Special build targets
############################################

clean:
	rm -f core *.o *~ $(TARGETS) ${GENSLIP}/*.o
