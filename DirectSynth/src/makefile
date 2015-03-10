TARGETS = direct_synth
HEADS = include.h structure.h function.h defs.h
OBJS = master.o sgt_handler.o task_manager.o worker.o \
	fp_cache.o misc.o synth.o mpi_functions.o log.o \
	iofunc.o 
SYNTH_OBJS = jbsim/jbsim3d.o jbsim/geoproj_subs.o jbsim/misc_subs.o jbsim/sgt3d_subs.o \
		jbsim/srf_rupm.o jbsim/srf_subs.o jbsim/stf_subs.o jbsim/fourg.o
PSA_OBJS = surfseis_rspectra.o setparm180.o sub_bandpass.o
ROTD_OBJS = rotd/calcrsp.o rotd/fftsub.o rotd/ft_th.o rotd/sort.o rotd/spline.o rotd/splint.o rotd/rotd.o rotd/rotd100.o

LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

GETPAR = ../../Getpar/getpar
LIBCFU = /work/00940/tera3d/CyberShake/utils/libcfu_0.03
LDLIBS = -L/work/00940/tera3d/CyberShake/software/RuptureCodes/RupGen-api-3.3.1/lib -lrupgen -L${GETPAR}/lib -lget -L/work/00940/tera3d/CyberShake/utils/libmemcached/lib -lmemcached -L${LIBCFU}/lib -lcfu -lm 
IFLAGS = -I/work/00940/tera3d/CyberShake/software/RuptureCodes/RupGen-api-3.3.1/include -I${GETPAR}/include -I/work/00940/tera3d/CyberShake/utils/libmemcached/include -I${LIBCFU}/include

UFLAGS = -O3 -g -Wuninitialized

CFLAGS = ${UFLAGS} ${IFLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132

CC = mpicc
FC = mpif90

all: $(TARGETS)
	for TARGET in $(TARGETS); do \
		rm -f ../bin/$$TARGET ; \
		cp $$TARGET ../bin/ ; \
	done

direct_synth : direct_synth.c ${OBJS} ${SYNTH_OBJS} ${PSA_OBJS} ${ROTD_OBJS}
	-rm -f direct_synth.o
	${CC} ${CFLAGS} -c -o direct_synth.o direct_synth.c
	${FC} ${FFLAGS} -o direct_synth direct_synth.o ${OBJS} ${SYNTH_OBJS} ${PSA_OBJS} ${ROTD_OBJS} ${LDLIBS} -lstdc++

clean:
	-rm -f $(OBJS) $(SYNTH_OBJS) $(PSA_OBJS) $(ROTD_OBJS) $(TARGETS)