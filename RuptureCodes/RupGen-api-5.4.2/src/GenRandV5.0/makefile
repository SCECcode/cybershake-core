HEADS = include.h structure.h function.h defs.h
OBJS = iofunc.o misc.o slip.o ruptime.o srf_subs.o gslip_srf_subs.o geoproj_subs.o gslip_sliprate_subs.o wafront2d-rwg.o fourg.o

GETPARLIB = ../../../../Getpar/getpar/lib
GETPARINC = ../../../../Getpar/getpar/include

##GCCLIB = /usr/lib/gcc/x86_64-redhat-linux/4.3.2
FORTLIBS = -lgfortranbegin -lgfortran -lgcc_s
FC = gfortran

##LDLIBS = ${OBJS} -L${GCCLIB} -L${GETPARLIB} -lm -lget ${FORTLIBS}
##LDLIBS = ${OBJS} -L${GCCLIB} -L${GETPARLIB} -lm -lget ${FORTLIBS} -I ${GETPARINC} -I ${FFTW_INCDIR} -L${FFTW_LIBDIR} -lfftw3f
LDLIBS = -L${GETPARLIB} -lm -lget ${FORTLIBS} -I${GETPARINC} -I${OLCF_FFTW_ROOT}/include -L${OLCF_FFTW_ROOT}/lib -lfftw3f

#LF_FLAGS = -D_FILE_OFFSET_BITS=32
#
# use following for large file capability
LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

UFLAGS = -O3 -g

CC = gcc

CFLAGS = ${UFLAGS} ${LF_FLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-none

# Uncomment for memcached support
MEMC = ../../../../../utils/libmemcached_1.0.18
CFLAGS := $(CFLAGS) -I${MEMC}/include -D_USE_MEMCACHED
LDLIBS := ${LDLIBS} -L${MEMC}/lib -lmemcached -L${OLCF_FFTW_ROOT}/lib

##### make options

genslip-v5.4.2.o : genslip-v5.4.2.c ${OBJS}
	${CC} ${CFLAGS} -c -o genslip-v5.4.2.o genslip-v5.4.2.c ${LDLIBS}

${OBJS} : ${HEADS}

clean :
	-rm -f ${OBJS}
