TARGETS = ucvm-mpi merge_cvm reformat_velocity_mpi ucvm-single-mpi determine_surface_model
GETPAR_HOME = ../../Getpar/getpar

#UCVM_HOME = /projects/sciteam/baln/CyberShake/software/UCVM/ucvm_bbp1d
UCVM_HOME = /projects/sciteam/baln/CyberShake/software/UCVM/ucvm-15.10.0
ETREE_HOME = ${UCVM_HOME}/lib/euclid3
PROJ4_HOME = ${UCVM_HOME}/lib/proj-4
VX_HOME = ${UCVM_HOME}/model/cvmh1510
CVMS_HOME = ${UCVM_HOME}/model/cvms
CVMSI_HOME = ${UCVM_HOME}/model/cvms426
CVMS5_HOME = ${UCVM_HOME}/model/cvms5
CCA_HOME = ${UCVM_HOME}/model/cca
CENCAL_HOME = ${UCVM_HOME}/model/cencal

UCVMLIBS = -lvxapi -L${VX_HOME}/lib -lgeo -L${VX_HOME}/lib -L${CVMS_HOME}/lib -lcvms -L${CVMSI_HOME}/lib -lcvmsi  -L${CENCAL_HOME}/lib -lcencalvm -L${PROJ4_HOME}/lib -lproj -L${ETREE_HOME}/lib -letree -L${UCVM_HOME}/lib -lucvm -lstdc++ -ldl
UCVMINC = -I${UCVM_HOME}/include -I${GETPAR_HOME}/include

UCVM_BBP_HOME = /projects/sciteam/baln/CyberShake/software/UCVM/ucvm_bbp1d
ETREE_BBP_HOME = ${UCVM_BBP_HOME}/lib/euclid3
PROJ4_BBP_HOME = ${UCVM_BBP_HOME}/lib/proj-4
UCVMBBPLIBS = -L${UCVM_BBP_HOME}/lib -lucvm -L${PROJ4_BBP_HOME}/lib -lproj -L${ETREE_BBP_HOME}/lib -letree -lm -L${GETPAR_HOME}/lib -lget -L/usr/lib64 -lgfortran
UCVMBBPINC = -I${UCVM_BBP_HOME}/include -I${GETPAR_HOME}/include

LDLIBS = -lm ${UCVMLIBS} -L${GETPAR_HOME}/lib -lget -L/usr/lib64 -lgfortran

include ../../Compilers.mk
CC = $(MY_MPICC)


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
CFLAGS = -O3 ${UCVMINC} -g -dynamic
#CFLAGS = -g -Wall ${UCVMINC}
FFLAGS = -O3 -ffixed-line-length-132 -finit-local-zero -dynamic
#FFLAGS = -g -ffixed-line-length-132 -finit-local-zero
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
	cp ucvm-single-mpi ../bin/

#ucvm-single-mpi-bbp1d : ucvm-single_mpi_bbp1d.c subs_mpi.o iofunc.o
#	${CC} -O3 -o ucvm-single-mpi-bbp1d ucvm-single_mpi_bbp1d.c subs_mpi.o iofunc.o ${UCVMBBPINC} ${UCVMBBPLIBS}

merge_cvm : merge_cvm.c subs_mpi.o iofunc.o
	${CC} ${CFLAGS} -o merge_cvm merge_cvm.c subs_mpi.o iofunc.o ${LDLIBS}

test_query : test_query.c
	${CC} ${CFLAGS} -dynamic -o test_query test_query.c ${UCVMINC} ${LDLIBS} ${UCVMLIBS}

determine_surface_model : determine_surface_model.c
	${CC} ${CFLAGS} -dynamic -o determine_surface_model determine_surface_model.c ${UCVMINC} ${LDLIBS} ${UCVMLIBS}
	cp determine_surface_model ../bin/

reformat_velocity_mpi : reformat_velocity_mpi.c iofunc.o
	${CC} -g -o reformat_velocity_mpi reformat_velocity_mpi.c iofunc.o

clean :
	-rm -f *.o ${TARGETS}
