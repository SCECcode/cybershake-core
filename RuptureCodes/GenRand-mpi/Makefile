GETPAR = ../../Getpar/getpar
GETPARLIB = ${GETPAR}/lib
GETPARINCL = ${GETPAR}/include


# General compiler/linker flags
CC=cc
FC=ftn
CFLAGS = -I${GETPARINCL} -Wall -O3 \
	-D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE \
	-D_FILE_OFFSET_BITS=64 -DUSE_IOBUF_MACROS
LIBS = -L${GETPARLIB} -lget -lgfortran -lm
#LDFLAGS = -static

TARGET = rupgen rupgen-serial

all: $(TARGET)

############################################
# Executables
############################################

rupgen: rupgen.o genslip.o func_mpi.o iofunc.o misc.o \
	ruptime.o slip.o srf_subs.o fourg.o rg_config.o
	$(CC) -o $@ $^ $(LIBS)

rupgen-serial: rupgen-serial.o genslip.o iofunc.o misc.o \
	ruptime.o slip.o srf_subs.o fourg.o
	$(CC) -o $@ $^ $(LIBS)


############################################
# Special build targets
############################################

clean:
	rm -f core *.o *~ $(TARGET)
