# Hostname of current system
HOST=`hostname -s`

# Home directory
HOME=/share/home/00940/tera3d
#HOME = /lustre/scratch/patricks
#HOME = /lustre/widow/scratch/spatrick
#HOME = /home/rcf-104/patrices

# Path to all CVM installations
CVM_DIR = /scratch/projects/tg/tera3d/CyberShake2007/software/UCVM

# Path to 3rd party libraries
#ETREE_DIR = $(HOME)/opt/$(HOST)/euclid3-1.2
P4_DIR = /scratch/projects/tg/tera3d/CyberShake2007/software/UCVM/proj4


# General compiler/linker flags
CFLAGS = -Wall -O3 -std=c99  \
	-D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE \
	-D_FILE_OFFSET_BITS=64
LDFLAGS = -static

# Default set of cvms to include
MODEL_TARGS = ucvm_model_1d.o ucvm_model_linthurber.o
#MODEL_TARGS = ucvm_model_1d.o


# CENCALVM includes
#CC_DIR = ${CVM_DIR}/cencalvm-0.6.6
#CC_INCLUDE = $(CC_DIR)/include/cencalvm/query
#CC_LIB = $(CC_DIR)/lib
#CFLAGS += -I$(CC_INCLUDE) -D_UCVM_ENABLE_CENCAL
#LDFLAGS += -L$(CC_LIB) -l cencalvm -lstdc++ -pthread
#MODEL_TARGS += ucvm_model_cencal.o

# Etree includes
#ETREE_INCLUDE = $(ETREE_DIR)/libsrc
#ETREE_LIB = $(ETREE_DIR)/libsrc
#CFLAGS += -I$(ETREE_INCLUDE)
#LDFLAGS += -L$(ETREE_LIB) -l etree

# Proj4 includes
P4_INCLUDE = $(P4_DIR)/include
P4_LIB = $(P4_DIR)/lib
CFLAGS += -I$(P4_INCLUDE)
LDFLAGS += -L$(P4_LIB) -lproj

# CVM-H includes
VX_DIR = ${CVM_DIR}/cvmh
VX_INCLUDE = $(VX_DIR)/include
VX_LIB = $(VX_DIR)/lib
GCTPC_LIB = $(VX_DIR)/gctpc/lib
CFLAGS += -I$(VX_INCLUDE) -D_UCVM_ENABLE_CVMH
LDFLAGS += -L$(VX_LIB) -lvxapi -L$(GCTPC_LIB) -lgeo
MODEL_TARGS += ucvm_model_cvmh.o

# CVM-S4 includes
S4_DIR = ${CVM_DIR}/cvms
S4_INCLUDE = $(S4_DIR)/include
S4_LIB = $(S4_DIR)/lib
CFLAGS += -I$(S4_INCLUDE) -D_UCVM_ENABLE_CVMS
LDFLAGS += -L$(S4_LIB) -lcvms -lgfortran
MODEL_TARGS += ucvm_model_cvms.o

# Common libs to include
LDFLAGS += -lm

# Kraken/Jaguar
#CC = cc
#FC=ftn
#LD = cc

# Ranger
#CC = mpicc
#FC=mpif77
#LD = mpicc

CC = gcc
FC = gfortran
LD = gcc

# HPCC
#CC = mpicc
#FC=gfortran
#LD = mpicc

# HPCC/Shock
#CC=gcc
#FC=gfortran
#LD=gcc


TARGET = ucvm_query libucvm.a

all: $(TARGET)


############################################
# Executables
############################################

test_ucvm: test_ucvm.o ucvm.o ucvm_utils.o ucvm_config.o libucvm.a
	$(LD) -o $@ $^ $(LDFLAGS)

test_grid: test_grid.o ucvm_grid.o ucvm.o ucvm_utils.o ucvm_config.o \
		ucvm_proj_bilinear.o libucvm.a
	$(LD) -o $@ $^ $(LDFLAGS)

ucvm_query: ucvm_query.o ucvm.o ucvm_config.o ucvm_utils.o libucvm.a
	$(LD) -o $@ $^ $(LDFLAGS)
	cp ucvm_query ../../bin

libucvm.a: ucvm.o ucvm_grid.o ucvm_config.o ucvm_proj_bilinear.o \
		ucvm_utils.o $(MODEL_TARGS)
	$(AR) rcs $@ $^
	cp $@ ../../lib
	cp ucvm.h ../../include
	cp ucvm_dtypes.h ../../include
	cp ucvm_grid.h ../../include
	cp ucvm_config.h ../../include


############################################
# Special build targets
############################################


clean:
	rm -f core *.o *~ $(TARGET) ../../include/*.h ../../lib/*.a ../../bin/ucvm_query

