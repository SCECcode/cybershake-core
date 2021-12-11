GETPAR = ../../../Getpar/getpar
GETPARLIB = ${GETPAR}/lib
GETPARINCL = ${GETPAR}/include

GENSLIP = GenRandV5.2.3

GENSLIP_OBJS = ${GENSLIP}/ruptime.o ${GENSLIP}/iofunc.o ${GENSLIP}/misc.o ${GENSLIP}/gslip_sliprate_subs.o ${GENSLIP}/geoproj_subs.o ${GENSLIP}/gslip_srf_subs.o ${GENSLIP}/wafront2d-rwg.o 

# General compiler/linker flags
CC=cc
FC=ftn
CFLAGS = -I${GETPARINCL} \
	-D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE \
	-D_FILE_OFFSET_BITS=64 -DUSE_IOBUF_MACROS -O3 -g
LIBS = -L${GETPARLIB} -lget -lgfortran -lm
LDFLAGS = -static

# Uncomment for memcached support
MEMC = /projects/sciteam/jmz/CyberShake/utils/libmemcached_1.0.18
CFLAGS := $(CFLAGS) -I${MEMC}/include -D_USE_MEMCACHED
LIBS := ${LIBS} -L${MEMC}/lib -lmemcached

TARGETS = librupgen.a

all: $(TARGETS)


############################################
# Executables
############################################

#librupgen.a: rupgen_api.o ${GENSLIP}/genslip-v3.3.1.o ${GENSLIP}/fourg.o ${GENSLIP}/slip.o ${GENSLIP}/ruptime.o ${GENSLIP}/iofunc.o ${GENSLIP}/srf_subs.o ${GENSLIP}/misc.o
librupgen.a: rupgen_api.o
	cd ${GENSLIP} && $(MAKE)
	$(AR) rcs $@ $^ ${GENSLIP}/genslip-v5.2.3a.o ${GENSLIP}/slip.o ${GENSLIP_OBJS}
	cp librupgen.a ../lib/
	cp rupgen_api.h ../include/
	cp structure.h ../include/

genslip: genslip.o iofunc.o srf_subs.o
	$(CC) -o $@ $^ $(LIBS)

############################################
# Special build targets
############################################

clean:
	rm -f core *.o *~ $(TARGETS) ${GENSLIP}/*.o
