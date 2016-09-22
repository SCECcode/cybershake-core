TARGETS = direct_synth
HEADS = include.h structure.h function.h defs.h duration.h
OBJS = master.o sgt_handler.o task_manager.o worker.o \
	fp_cache.o misc.o synth.o mpi_functions.o log.o \
	iofunc.o
SYNTH_OBJS = jbsim/jbsim3d.o jbsim/geoproj_subs.o jbsim/misc_subs.o jbsim/sgt3d_subs.o \
		jbsim/srf_rupm.o jbsim/srf_subs.o jbsim/stf_subs.o jbsim/fourg.o
PSA_OBJS = surfseis_rspectra.o setparm180.o sub_bandpass.o
ROTD_OBJS = rotd/calcrsp.o rotd/fftsub.o rotd/ft_th.o rotd/sort.o rotd/spline.o rotd/splint.o rotd/rotd.o rotd/rotd100.o
DURATION_OBJS = duration.o integ_diff.o

LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

BASE_DIR = /projects/sciteam/jmz/CyberShake
GETPAR = ../../Getpar/getpar
LIBCFU = ${BASE_DIR}/utils/libcfu_0.03
LDLIBS = -L${BASE_DIR}/software/RuptureCodes/RupGen-api-5.2.3/lib -lrupgen -L${GETPAR}/lib -lget -L${BASE_DIR}/utils/libmemcached_1.0.18/lib -lmemcached -L${LIBCFU}/lib -lcfu -lm 
IFLAGS = -I${BASE_DIR}/software/RuptureCodes/RupGen-api-5.2.3/include -I${GETPAR}/include -I${BASE_DIR}/utils/libmemcached_1.0.18/include -I${LIBCFU}/include

UFLAGS = -O3 -g -Wuninitialized -dynamic

CFLAGS = ${UFLAGS} ${IFLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132
#FFLAGS = ${UFLAGS}

CC = cc
FC = ftn

all: $(TARGETS)
	for TARGET in $(TARGETS); do \
		rm -f ../bin/$$TARGET ; \
		cp $$TARGET ../bin/ ; \
	done

direct_synth : direct_synth.c ${OBJS} ${SYNTH_OBJS} ${PSA_OBJS} ${ROTD_OBJS} ${DURATION_OBJS}
	-rm -f direct_synth.o
	${CC} ${CFLAGS} -c -o direct_synth.o direct_synth.c
	${FC} ${FFLAGS} -o direct_synth direct_synth.o ${OBJS} ${SYNTH_OBJS} ${PSA_OBJS} ${ROTD_OBJS} ${DURATION_OBJS} ${LDLIBS} -lstdc++

clean:
	-rm -f $(OBJS) $(SYNTH_OBJS) $(PSA_OBJS) $(ROTD_OBJS) $(DURATION_OBJS) $(TARGETS)
