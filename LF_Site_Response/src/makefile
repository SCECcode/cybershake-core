HEADS = structure.h
OBJS = wcc_siteamp14.o get_vs30.o iofunc.o misc.o integ_diff.o wcc_getpeak.o fft1d.o

GETPAR_DIR = ../../Getpar/getpar
RUPGEN_DIR = /projects/sciteam/jmz/CyberShake/software/RuptureCodes/RupGen-api-3.3.1
MEMCACHED_DIR = /projects/sciteam/jmz/CyberShake/utils/libmemcached_1.0.18

UCVM_HOME = /projects/sciteam/jmz/CyberShake/software/UCVM/ucvm_14.3.0
ETREE_HOME = ${UCVM_HOME}/lib/euclid3
PROJ4_HOME = ${UCVM_HOME}/lib/proj-4
VX_HOME = ${UCVM_HOME}/model/cvmh1191
CVMS_HOME = ${UCVM_HOME}/model/cvms4
CVMSI_HOME = ${UCVM_HOME}/model/cvms426

UCVMLIBS = -L${UCVM_HOME}/lib -lucvm -lvxapi -L${VX_HOME}/lib -lgeo -L${VX_HOME}/lib -L${CVMS_HOME}/lib -lcvms -L${CVMSI_HOME}/lib -lcvmsi -L${PROJ4_HOME}/lib -lproj -L${ETREE_HOME}/lib -letree
UCVMINC = -I${UCVM_HOME}/include -I${GETPAR_HOME}/include

LIBS = -lm ${GETPAR_DIR}/lib/libget.a
IFLAGS = -I${GETPAR_DIR}/include -I${RUPGEN_DIR}/include -I${MEMCACHED_DIR}/include
LDLIBS = -L${RUPGEN_DIR}/lib -lrupgen -L${MEMCACHED_DIR}/lib -lmemcached ${LIBS} -lstdc++ 

#LF_FLAGS = -D_FILE_OFFSET_BITS=32
#
# use following for large file capability
LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

UFLAGS = -O3 -Wimplicit-function-declaration -Wl,-upthread_mutex_destroy
#UFLAGS = -g

FC = ftn
CC = cc

CFLAGS = ${UFLAGS} ${LF_FLAGS} ${IFLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132 -finit-local-zero -Wuninitialized
#FFLAGS = ${UFLAGS} -N 132 -e 0 
##### make options

all: lf_site_response

test : test.c
	${CC} ${CFLAGS} -o test test.c

lf_site_response : lf_site_response.c ${OBJS}
	${CC} ${CFLAGS} -o lf_site_response lf_site_response.c ${OBJS} ${LDLIBS} ${UCVMLIBS}
	cp lf_site_response ../bin/

get_vs30.o : get_vs30.c
	${CC} ${CFLAGS} -c -o get_vs30.o get_vs30.c ${UCVMINC}

get_vs30 : get_vs30.o
	${CC} ${CFLAGS} -o get_vs30 get_vs30.o ${UCVMLIBS}

clean :
	rm -f ${OBJS} *.o
