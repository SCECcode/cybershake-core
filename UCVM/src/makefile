TARGETS = ucvm-mpi merge_cvm ucvm-single-mpi determine_surface_model shear_modulus_from_srf ucvm-single-mpi-18.5
GETPAR_HOME = ../../Getpar/getpar

#Switch to 18.5 as default
UCVM_DIR = /lustre/atlas/proj-shared/geo112/CyberShake/software/UCVM
UCVM_15_10_HOME = ${UCVM_DIR}/ucvm-15.10.0
VX_15_10_HOME = ${UCVM_15_10_HOME}/model/cvmh1510
CVMS_15_10_HOME = ${UCVM_15_10_HOME}/model/cvms
CVMSI_15_10_HOME = ${UCVM_15_10_HOME}/model/cvms426
CVMS5_15_10_HOME = ${UCVM_15_10_HOME}/model/cvms5
CCA_15_10_HOME = ${UCVM_15_10_HOME}/model/cca
CENCAL_15_10_HOME = ${UCVM_15_10_HOME}/model/cencal
ETREE_15_10_HOME = ${UCVM_15_10_HOME}/lib/euclid3
PROJ4_15_10_HOME =${UCVM_15_10_HOME}/lib/proj-4

#Libs for ucvm 18.5.0
#UCVM_HOME = /lustre/atlas/proj-shared/geo112/ucvm_18_5/install
UCVM_HOME = ${UCVM_DIR}/ucvm-18.5.0
VX_HOME = ${UCVM_HOME}/model/cvmh1510
CVMS_HOME = ${UCVM_HOME}/model/cvms
CVMSI_HOME = ${UCVM_HOME}/model/cvms426
CVMS5_HOME = ${UCVM_HOME}/model/cvms5
CCA_HOME = ${UCVM_HOME}/model/cca
CENCAL_HOME = ${UCVM_HOME}/model/cencal
CS173_HOME = ${UCVM_HOME}/model/cs173
ETREE_HOME = ${UCVM_HOME}/lib/euclid3
PROJ4_HOME = ${UCVM_HOME}/lib/proj-4

#UCVMLIBS = -lucvm -L${UCVM_HOME}/lib -lvxapi -L${VX_HOME}/lib -lgeo -L${VX_HOME}/gctpc/lib -L${CVMS_HOME}/lib -lcvms -L${CVMSI_HOME}/lib -lcvmsi -L${CENCAL_HOME}/lib -lcencalvm -L${ETREE_HOME}/lib -letree -L${PROJ4_HOME}/lib -lproj -lstdc++
UCVMLIBS = -lucvm -L${UCVM_HOME}/lib -lvxapi -L${VX_HOME}/lib -lgeo -L${VX_HOME}/lib -L${CVMS_HOME}/lib -lcvms -L${CVMSI_HOME}/lib -lcvmsi  -L${CENCAL_HOME}/lib -lcencalvm -L${PROJ4_HOME}/lib -L${CCA_HOME}/lib -lcca -L${CS173_HOME}/lib -lcs173 -L/${CVMS5_HOME}/lib -lcvms5 -lproj -L${ETREE_HOME}/lib -letree -lstdc++ -ldl
UCVM_15_10_LIBS = -lucvm -L${UCVM_15_10_HOME}/lib -lvxapi -L${VX_15_10_HOME}/lib -lgeo -L${VX_15_10_HOME}/lib -L${CVMS_15_10_HOME}/lib -lcvms -L${CVMSI_15_10_HOME}/lib -lcvmsi  -L${CENCAL_15_10_HOME}/lib -lcencalvm -L${CCA_15_10_HOME}/lib -lcca -L${PROJ4_15_10_HOME}/lib -lproj -L${ETREE_15_10_HOME}/lib -letree -lstdc++ -ldl

UCVMINC = -I${UCVM_HOME}/include -I${GETPAR_HOME}/include

UCVM_15_10_INC = -I${UCVM_15_10_HOME}/include -I${GETPAR_HOME}/include

LDLIBS = -lm -L${GETPAR_HOME}/lib -lget -lgfortran

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
CFLAGS = -O3 -g -static
FFLAGS = -O3 -ffixed-line-length-132 -finit-local-zero
#FFLAGS = -O2

##### make options

all: $(TARGETS)
	for TARGET in $(TARGETS); do \
        	rm -rf ../bin/$$TARGET ; \
        	cp $$TARGET ../bin/ ; \
	done

check_media : check_media.c iofunc.o
	${CC} ${CFLAGS} -o check_media check_media.c iofunc.o ${LDLIBS}
ucvm-mpi : ucvm-main_mpi.c subs_mpi.o iofunc.o
	${CC} ${CFLAGS} -o ucvm-mpi ucvm-main_mpi.c subs_mpi.o iofunc.o ${UCVMINC} ${LDLIBS} ${UCVMLIBS}

ucvm-single-mpi-15.10 : ucvm-single_mpi-15.10.c subs_mpi.o iofunc.o
	${CC} ${CFLAGS} -o ucvm-single-mpi-15.10 ucvm-single_mpi-15.10.c subs_mpi.o iofunc.o ${UCVM_15_10_INC} ${LDLIBS} ${UCVM_15_10_LIBS}
	cp ucvm-single-mpi-15.10 ../bin/
	cp ucvm-single-mpi ../bin/

#ucvm-single-mpi-bbp1d : ucvm-single_mpi_bbp1d.c subs_mpi.o iofunc.o
#	${CC} -O3 -o ucvm-single-mpi-bbp1d ucvm-single_mpi_bbp1d.c subs_mpi.o iofunc.o ${UCVMBBPINC} ${UCVMBBPLIBS}

merge_cvm : merge_cvm.c subs_mpi.o iofunc.o
	${CC} ${CFLAGS} -o merge_cvm merge_cvm.c subs_mpi.o iofunc.o ${LDLIBS}

test_query : test_query.c
	${CC} ${CFLAGS} -o test_query test_query.c ${UCVMINC} ${LDLIBS} ${UCVMLIBS}

determine_surface_model : determine_surface_model.c
	${CC} ${CFLAGS} -o determine_surface_model determine_surface_model.c ${UCVMINC} ${LDLIBS} ${UCVMLIBS}
	cp determine_surface_model ../bin/

determine_surface_model_mpi : determine_surface_model_mpi.c
	${CC} ${CFLAGS} -o determine_surface_model_mpi determine_surface_model_mpi.c ${UCVMINC} ${LDLIBS} ${UCVMLIBS}
	cp determine_surface_model_mpi ../bin/

shear_modulus_from_srf : shear_modulus_from_srf.c
	${CC} ${CFLAGS} -o shear_modulus_from_srf shear_modulus_from_srf.c ${LDLIBS} -lm -I../../RuptureCodes/RupGen-api-3.3.1/include -lrupgen -L../../RuptureCodes/RupGen-api-3.3.1/lib -lmemcached -L/lustre/atlas/proj-shared/geo112/CyberShake/utils/libmemcached_1.0.18/lib -lstdc++
	cp shear_modulus_from_srf ../bin/

ucvm-single-mpi : ucvm-single_mpi.c subs_mpi.o iofunc.o
	${CC} ${CFLAGS} -o ucvm-single-mpi ucvm-single_mpi.c subs_mpi.o iofunc.o ${UCVMINC} ${LDLIBS} ${UCVMLIBS}
	cp ucvm-single-mpi ../bin/

get_model_info_for_db: get_model_info_for_db.c
	${CC} ${CFLAGS} -o get_model_info_for_db get_model_info_for_db.c ${UCVMINC} ${LDLIBS} ${UCVMLIBS}
	cp get_model_info_for_db ../bin/

%.o : %.c
	${CC} ${CLFAGS} -c $< -o $@ -I${GETPAR_HOME}/include

clean :
	-rm -f *.o ${TARGETS}
