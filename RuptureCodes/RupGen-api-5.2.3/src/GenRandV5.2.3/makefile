HEADS = include.h structure.h function.h defs.h
OBJS = iofunc.o misc.o slip.o ruptime.o gslip_srf_subs.o gslip_sliprate_subs.o geoproj_subs.o
FDRT_OBJS = wafront2d-rwg.o

GETPAR = ../../../../Getpar/getpar

LIBS = -lm ${GETPAR}/lib/libget.a
LIBS := ${OBJS} ${LIBS}

#LF_FLAGS = -D_FILE_OFFSET_BITS=32
#
# use following for large file capability
LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64
#LF_FLAGS = 
UFLAGS = -O3 -g

CC = cc
FC = ftn

CFLAGS = ${UFLAGS} ${LF_FLAGS} -I${GETPAR}/include
FFLAGS = ${UFLAGS} -ffixed-line-length-none

#Uncomment for memcached support
MEMC = /projects/sciteam/jmz/CyberShake/utils/libmemcached_1.0.18
CFLAGS := $(CFLAGS) -I${MEMC}/include -D_USE_MEMCACHED
LIBS := ${LIBS} -L${MEMC}/lib -lmemcached

##### make options

#all : genslip-v5.2.2 srf_gethypo get_rvfrac
all : genslip-v5.2.3a.o

genslip-v5.2.3a.o : genslip-v5.2.3a.c ${OBJS} ${SRF_OBJS} ${FDRT_OBJS}
	${CC} ${CFLAGS} -c -o genslip-v5.2.3a.o genslip-v5.2.3a.c -lfftw3f ${LIBS} 

genslip-v5.2.2.o : genslip-v5.2.2.c ${OBJS} ${SRF_OBJS} ${FDRT_OBJS}
	${CC} ${CFLAGS} -c -o genslip-v5.2.2.o genslip-v5.2.2.c -lfftw3f ${LIBS} 
	#cp genslip-v5.2.2 ../bin/.

srf_gethypo : srf_gethypo.c
	${CC} ${CFLAGS} -o srf_gethypo srf_gethypo.c ${LDLIBS} ${INCPAR} ${SRF_OBJS}
	cp srf_gethypo ../bin/.

get_rvfrac : get_rvfrac.c ${OBJS}
	${CC} ${CFLAGS} -o get_rvfrac get_rvfrac.c ${LDLIBS} ${INCPAR} ${SRF_OBJS}

${OBJS} : ${HEADS}

clean :
	-rm -f *.o genslip-v5.2.2 srf_gethypo get_rvfrac
