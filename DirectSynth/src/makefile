TARGETS = direct_synth
HEADS = include.h structure.h function.h defs.h duration.h
OBJS = master.o sgt_handler.o task_manager.o worker.o \
	fp_cache.o misc.o synth.o mpi_functions.o log.o \
	iofunc.o
SYNTH_OBJS = jbsim/geoproj_subs.o jbsim/misc_subs.o jbsim/sgt3d_subs.o \
		jbsim/srf_rupm.o jbsim/srf_subs.o jbsim/stf_subs.o jbsim/fourg.o
PSA_OBJS = surfseis_rspectra.o setparm180.o sub_bandpass.o
ROTD_OBJS = rotd/calcrsp.o rotd/fftsub.o rotd/ft_th.o rotd/sort.o rotd/spline.o rotd/splint.o rotd/rotd.o rotd/rotd100.o
DURATION_OBJS = duration.o integ_diff.o
TEST_DURATION_OBJS = ${DURATION_OBJS} iofunc.o

LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

BASE_DIR = /gpfs/alpine/proj-shared/geo112/CyberShake
GETPAR = ../../Getpar/getpar
LIBCFU = ${BASE_DIR}/utils/libcfu_0.03
CUR_RG_LIB = -L${BASE_DIR}/software/RuptureCodes/RupGen-api-5.5.2/lib
V5_4_2_RG_LIB = -L${BASE_DIR}/software/RuptureCodes/RupGen-api-5.4.2/lib
V3_3_1_RG_LIB = -L${BASE_DIR}/software/RuptureCodes/RupGen-api-3.3.1/lib
LDLIBS = -lrupgen -L${GETPAR}/lib -lget -L${BASE_DIR}/utils/libmemcached_1.0.18/lib -lmemcached -L${LIBCFU}/lib -lcfu -lm -L${OLCF_FFTW_ROOT}/lib -lfftw3f
CUR_RG_INC = -I${BASE_DIR}/software/RuptureCodes/RupGen-api-5.5.2/include
V5_4_2_RG_INC = -I${BASE_DIR}/software/RuptureCodes/RupGen-api-5.4.2/include
V3_3_1_RG_INC = -I${BASE_DIR}/software/RuptureCodes/RupGen-api-3.3.1/include
IFLAGS = -I${GETPAR}/include -I${BASE_DIR}/utils/libmemcached_1.0.18/include -I${LIBCFU}/include 

#UFLAGS = -O3 -Wuninitialized -dynamic -g
UFLAGS = -O3 -g 

CFLAGS = ${UFLAGS} ${IFLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132
#FFLAGS = ${UFLAGS} -qfixed=132

CC = mpicc
FC = mpif77

direct_synth_v3.3.1 : CFLAGS = ${UFLAGS} ${IFLAGS} -fno-expensive-optimizations -D_V3_3_1

all: $(TARGETS)
	for TARGET in $(TARGETS); do \
		rm -f ../bin/$$TARGET ; \
		cp $$TARGET ../bin/ ; \
	done

direct_synth : direct_synth.c jbsim/jbsim3d.o ${OBJS} ${SYNTH_OBJS} ${PSA_OBJS} ${ROTD_OBJS} ${DURATION_OBJS}
	-rm -f direct_synth.o
	${CC} ${CFLAGS} ${CUR_RG_INC} -c -o direct_synth.o direct_synth.c
	${FC} ${FFLAGS} -o direct_synth direct_synth.o jbsim/jbsim3d.o ${OBJS} ${SYNTH_OBJS} ${PSA_OBJS} ${ROTD_OBJS} ${DURATION_OBJS} ${CUR_RG_LIB} ${LDLIBS} -lstdc++
	cp direct_synth ../bin/

direct_synth_omp: direct_synth.c jbsim/jbsim3d_omp.o ${OBJS} ${SYNTH_OBJS} ${PSA_OBJS} ${ROTD_OBJS} ${DURATION_OBJS}
	-rm -f direct_synth.o
	${CC} ${CFLAGS} ${CUR_RG_INC} -c -o direct_synth.o direct_synth.c 
	${FC} ${FFLAGS} -o direct_synth_omp direct_synth.o jbsim/jbsim3d_omp.o ${OBJS} ${SYNTH_OBJS} ${PSA_OBJS} ${ROTD_OBJS} ${DURATION_OBJS} ${CUR_RG_LIB} ${LDLIBS} -lstdc++ -fopenmp
	cp direct_synth_omp ../bin/

direct_synth_v5.4.2 : direct_synth.c jbsim/jbsim3d_v5.4.2.o ${OBJS} ${SYNTH_OBJS} ${PSA_OBJS} ${ROTD_OBJS} ${DURATION_OBJS}
	-rm -f direct_synth_v5.4.2.o
	${CC} ${CFLAGS} ${V5_4_2_RG_INC} -c -o direct_synth_v5.4.2.o direct_synth.c
	${FC} ${FFLAGS} -o direct_synth_v5.4.2 direct_synth_v5.4.2.o jbsim/jbsim3d_v5.4.2.o ${OBJS} ${SYNTH_OBJS} ${PSA_OBJS} ${ROTD_OBJS} ${DURATION_OBJS} ${V5_4_2_RG_LIB} ${LDLIBS} -lstdc++
	cp direct_synth_v5.4.2 ../bin/

direct_synth_v3.3.1 : direct_synth.c jbsim/jbsim3d_v3.3.1.o ${OBJS} ${SYNTH_OBJS} ${PSA_OBJS} ${ROTD_OBJS} ${DURATION_OBJS}
	-rm -f direct_synth_v3.3.1.o
	${CC} ${CFLAGS} ${V3_3_1_RG_INC} -D_V3_3_1 -c -o direct_synth_v3.3.1.o direct_synth.c
	${FC} ${FFLAGS} -D_V3_3_1 -o direct_synth_v3.3.1 direct_synth_v3.3.1.o jbsim/jbsim3d_v3.3.1.o ${OBJS} ${SYNTH_OBJS} ${PSA_OBJS} ${ROTD_OBJS} ${DURATION_OBJS} ${V3_3_1_RG_LIB} ${LDLIBS} -lstdc++
	cp direct_synth_v3.3.1 ../bin/

jbsim/jbsim3d.o : jbsim/jbsim3d.c
	${CC} ${CFLAGS} ${CUR_RG_INC} -c -o jbsim/jbsim3d.o jbsim/jbsim3d.c

jbsim/jbsim3d_omp.o : jbsim/jbsim3d_omp.c
	${CC} ${CFLAGS} ${CUR_RG_INC} -c -o jbsim/jbsim3d_omp.o jbsim/jbsim3d_omp.c -fopenmp

jbsim/jbsim3d_v5.4.2.o : jbsim/jbsim3d.c
	${CC} ${CFLAGS} ${V5_4_2_RG_INC} -c -o jbsim/jbsim3d_v5.4.2.o jbsim/jbsim3d.c

jbsim/jbsim3d_v3.3.1.o : jbsim/jbsim3d.c
	${CC} ${CFLAGS} ${V3_3_1_RG_INC} -D_V3_3_1 -c -o jbsim/jbsim3d_v3.3.1.o jbsim/jbsim3d.c

test_duration : test_duration.c ${TEST_DURATION_OBJS}
	${CC} ${CFLAGS} -o test_duration test_duration.c ${TEST_DURATION_OBJS} ${LDLIBS}

clean:
	-rm -f $(OBJS) $(SYNTH_OBJS) $(PSA_OBJS) $(ROTD_OBJS) $(DURATION_OBJS) jbsim/jbsim3d.o jbsim/jbsim3d_v3.3.1.o $(TARGETS) direct_synth.o direct_synth_v3.3.1.o direct_synth_v3.3.1
