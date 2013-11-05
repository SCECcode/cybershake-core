TARGETS = ucvm-mpi merge_cvm ucvm-single-mpi
GETPAR_HOME = ../../Getpar/getpar
#LDLIBS =     ${GETPAR_HOME}/lib/libget.a -lm 


UCVM_HOME = /work/00940/tera3d/CyberShake/software/UCVM/ucvm_12.2.0
VX_HOME = /work/00940/tera3d/CyberShake/software/UCVM/cvmh_11.9.1
CVMS_HOME = /work/00940/tera3d/CyberShake/software/UCVM/cvms
ETREE_HOME = /work/00940/tera3d/CyberShake/software/UCVM/euclid3-1.3
PROJ4_HOME = /work/00940/tera3d/CyberShake/software/UCVM/proj_4.7.0

UCVMLIBS = -lucvm -L${UCVM_HOME}/lib -lvxapi -L${VX_HOME}/lib -lgeo -L${VX_HOME}/gctpc/lib -L${CVMS_HOME}/lib -lcvms -L${ETREE_HOME}/libsrc -letree -L${PROJ4_HOME}/lib -lproj
UCVMINC = -I${UCVM_HOME}/include


LDLIBS = -lm -L${GETPAR_HOME}/lib -lget -lgfortran

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
CFLAGS = -O3
FFLAGS = -O3 -ffixed-line-length-132 -finit-local-zero
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
	${CC} ${CFLAGS} -o ucvm-mpi ucvm-main_mpi.c subs_mpi.o iofunc.o ${UCVMINC} ${LDLIBS} ${UCVMLIBS}

ucvm-single-mpi : ucvm-single_mpi.c subs_mpi.o iofunc.o
	${CC} ${CFLAGS} -o ucvm-single-mpi ucvm-single_mpi.c subs_mpi.o iofunc.o ${UCVMINC} ${LDLIBS} ${UCVMLIBS}

merge_cvm : merge_cvm.c subs_mpi.o iofunc.o
	${CC} ${CFLAGS} -o merge_cvm merge_cvm.c subs_mpi.o iofunc.o ${LDLIBS}

clean :
	-rm -f *.o ${TARGETS}
