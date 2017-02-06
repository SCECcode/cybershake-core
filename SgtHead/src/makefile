TARGETS = write_head compare_head strip_head reformat_awp reformat_velocity reformat_awp_mpi
HEADS = defs.h include.h structure.h function.h
#OBJS =  absorb.o diff.o fft1d.o  interface.o iofunc.o \
        memory_slices.o  misc.o tstep.o \
        genmodel.o memory_variables.o source.o \
        geoproj_subs.o subs_mpi.o fourg.o
WRITE_OBJS = head_subs.o write_head.o
CMP_OBJS = head_subs.o compare_head.o
REF_OBJS = reformat_awp.o
STR_OBJS = strip_head.o
VEL_OBJS = reformat_velocity.o

GETPAR_HOME = ../../Getpar/getpar
LDLIBS =     ${GETPAR_HOME}/lib/libget.a -lm

UFLAGS = -g

include ../../Compilers.mk
CC = $(MY_MPICC)
FC = $(MY_MPIFC)

#CC = mpicc
#FC = mpif77

##### If hardware and OS supports 64-bit use the following to
#     enable option of opening large (>2**31 bytes) files
#
#     This option should be used for solaris 2.6 and later
#
#####LF_FLAGS = -D_FILE_OFFSET_BITS=64 `getconf LFS_LIBS` `getconf LFS_LDFLAGS` `getconf LFS_CFLAGS`
LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

##### If hardware and OS DOES NOT support 64-bit, use the following to
#     disable option of opening large (>2**31 bytes) files
#
#     This option should be used for solaris 2.5 and earlier
#
###LF_FLAGS = -D_FILE_OFFSET_BITS=32

##### C options

# optimized
CFLAGS = -O3  -h msglevel_3 ${LF_FLAGS}
FFLAGS = -h msglevel_3 ${UFLAGS}

##### make options

all: $(TARGETS)
	for TARGET in $(TARGETS); do \
		rm -rf ../bin/$$TARGET ; \
		cp $$TARGET ../bin ; \
	done

write_head : ${WRITE_OBJS}
	${CC} ${CFLAGS} -o write_head ${WRITE_OBJS} ${LDLIBS}

compare_head : ${CMP_OBJS}
	${CC} ${CFLAGS} -o compare_head ${CMP_OBJS} ${LDLIBS}

reformat_awp : ${REF_OBJS}
	${CC} ${CLFAGS} -o reformat_awp ${REF_OBJS}

reformat_awp_mpi : reformat_awp_mpi.o
	${CC} ${CLFAGS} -o reformat_awp_mpi reformat_awp_mpi.o ${LDLIBS}

reformat_velocity : ${VEL_OBJS}
	${CC} ${CFLAGS} -o reformat_velocity ${VEL_OBJS}

clean :
	-rm -f ${OBJS} ${WRITE_OBJS} ${CMP_OBJS} ${REF_OBJS} ${TARGETS}
