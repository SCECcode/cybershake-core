HEADS = include.h structure.h function.h defs.h params.h
SRF2STOCH_OBJS = iofunc.o misc.o srf_subs.c
HFSIM_OBJS = iofunc.o misc.o srf_subs.o hb_high_v5.4.3.o integ_diff.o wcc_getpeak.o wcc_rotate.o wcc_siteamp14.o fft1d.o
MERGE_OBJS = iofunc.o misc.o wcc_tfilter.o wcc_resamp_arbdt.o wcc_add.o fourg.o
HF_SYNTH_OBJS = srf2stoch.o ${HFSIM_OBJS} hfsim.o get_vs30.o
SITEAMP_OBJS = iofunc.o misc.o fft1d.o

GETPAR_DIR = ../../Getpar/getpar
RUPGEN_DIR = /projects/sciteam/jmz/CyberShake/software/RuptureCodes/RupGen-api-3.3.1
MEMCACHED_DIR = /projects/sciteam/jmz/CyberShake/utils/libmemcached_1.0.15

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

UFLAGS = -O3 -g

FC = ftn
CC = cc

CFLAGS = ${UFLAGS} ${LF_FLAGS} ${IFLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132 -finit-local-zero -Wuninitialized
#FFLAGS = ${UFLAGS} -extend-source -nofor-main -zero

##### make options

all: hf_synth

hf_synth : hf_synth.o ${HF_SYNTH_OBJS}
	${FC} ${FFLAGS} -o hf_synth hf_synth.o ${HF_SYNTH_OBJS} ${LDLIBS} ${UCVMLIBS}
	cp hf_synth ../bin/

merge : merge.c ${MERGE_OBJS}
	${CC} ${CFLAGS} -o merge merge.c ${MERGE_OBJS} ${LDLIBS}
#	cp merge ../bin/

hfsim : hfsim.o ${HFSIM_OBJS}
	${FC} ${FFLAGS} -o hfsim hfsim.o ${HFSIM_OBJS} ${LDLIBS}
#	cp hfsim ../bin/

srf2stoch : srf2stoch.c ${SRF2STOCH_OBJS}
	${CC} ${CFLAGS} -o srf2stoch srf2stoch.c ${SRF2STOCH_OBJS} ${LDLIBS}
#	cp srf2stoch ../bin/

#hb_high_v3.1.1 : hb_high_v3.1.1.f
#	${FC} ${FFLAGS} -o hb_high_v3.1.1 hb_high_v3.1.1.f
#	cp hb_high_v3.1.1 ../bin/

wcc_siteamp14 : wcc_siteamp14_standalone.c ${SITEAMP_OBJS}
	${CC} ${CFLAGS} -o wcc_siteamp14 wcc_siteamp14_standalone.c ${SITEAMP_OBJS} ${LDLIBS}

get_vs30.o : get_vs30.c
	${CC} ${CFLAGS} -c -o get_vs30.o get_vs30.c ${UCVMINC}

clean :
	rm -f ${OBJS} ${GENRAND_OBJS} ${BAILEY_OBJS} ${GEOPROJ_OBJS} *.o
