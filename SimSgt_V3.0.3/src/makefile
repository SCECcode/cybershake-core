HEADS = defs.h include.h structure.h function.h
OBJS =  absorb.o diff.o fft1d.o  interface.o iofunc.o \
        memory_slices.o  misc.o tstep.o \
        genmodel.o memory_variables.o source.o \
        geoproj_subs.o subs_mpi.o fourg.o
MAIN_OBJ =  main_mpi.o

GETPARLIB = ../..//Getpar/getpar/lib
LDLIBS = -L${GETPARLIB} -lm -lget

UFLAGS = -O3
CC = mpicc
FC = mpif77
#CC = gcc
#FC = gfortran

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
CFLAGS = ${UFLAGS} ${LF_FLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132

##### make options

emod3d-mpi : ${MAIN_OBJ} ${OBJS}
	${CC} ${CFLAGS} -o emod3d-mpi ${MAIN_OBJ} ${OBJS} ${LDLIBS}
	cp emod3d-mpi ../bin/

emod3d-mpiX : ${MAIN_OBJ} ${OBJS}
	${CC} ${CFLAGS} -o emod3d-mpiX ${MAIN_OBJ} ${OBJS} ${LDLIBS}

get_nodeinfo : get_nodeinfo.c ${OBJS}
	${CC} ${CFLAGS} get_nodeinfo.c -o get_nodeinfo ${OBJS} ${LDLIBS}

check_n1n2 : check_n1n2.c ${OBJS}
	${CC} ${CFLAGS} check_n1n2.c -o check_n1n2 ${OBJS} ${LDLIBS}

gen_source : gen_source.c ${OBJS}
	${CC} ${CFLAGS} gen_source.c -o gen_source ${OBJS} ${LDLIBS}

${OBJS} : ${HEADS}

${MAIN_OBJ} : ${HEADS}

clean :
	-rm -f ${MAIN_OBJ} ${OBJS}
