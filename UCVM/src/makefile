TARGETS = ucvm-mpi merge_cvm
GETPAR_HOME = ../../Getpar/getpar
#LDLIBS =     ${GETPAR_HOME}/lib/libget.a -lm 


UCVM_HOME = /scratch/projects/tg/tera3d/CyberShake2007/software-9.2.0/UCVM/ucvm
VX_HOME = /scratch/projects/tg/tera3d/CyberShake2007/software-9.2.0/UCVM/cvmh
CVMS_HOME = /scratch/projects/tg/tera3d/CyberShake2007/software-9.2.0/UCVM/cvms
UCVMLIBS = -lvxapi -L${VX_HOME}/lib -lgeo -L${VX_HOME}/gctpc/lib -L${UCVM_HOME}/lib -lucvm -L${CVMS_HOME}/lib -lcvms
UCVMINC = -I${UCVM_HOME}/include


LDLIBS = -lm ${UCVMLIBS} -L${GETPAR_HOME}/lib -lget -lgfortran

include ../../Compilers.mk
CC = $(MY_MPICC)

#mpicc -O1 -I/scratch/projects/tg/tera3d/CyberShake2007/software-9.2.0/VH-WrapC/cvmh/include -L/scratch/projects/tg/tera3d/CyberShake2007/software-9.2.0/VH-WrapC/cvmh/lib -L/scratch/projects/tg/tera3d/CyberShake2007/software-9.2.0/VH-WrapC/cvmh/gctpc/lib -lvxapi -lgeo -L../../Getpar/getpar/lib -lget -o verH-mpi verH-main_mpi.c subs_mpi.o iofunc.o

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
#CFLAGS = -O1 -compile-info -link-info ${VXINC}
CFLAGS = -O1 ${UCVMINC}
FFLAGS = -O2 -ffixed-line-length-132 -finit-local-zero
#FFLAGS = -O2

##### make options

all: $(TARGETS)
	for TARGET in $(TARGETS); do \
        	rm -rf ../bin/$$TARGET ; \
        	cp $$TARGET ../bin ; \
	done

check_media : check_media.c iofunc.o
	${CC} ${CFLAGS} -o check_media check_media.c iofunc.o ${LDLIBS}
ucvm-mpi : ucvm-main_mpi.c subs_mpi.o iofunc.o
	${CC} ${CFLAGS} -o ucvm-mpi ucvm-main_mpi.c subs_mpi.o iofunc.o ${LDLIBS} ${UCVMLIBS}

merge_cvm : merge_cvm.c subs_mpi.o iofunc.o
	${CC} ${CFLAGS} -o merge_cvm merge_cvm.c subs_mpi.o iofunc.o ${LDLIBS}

clean :
	-rm -f *.o ${TARGETS}
