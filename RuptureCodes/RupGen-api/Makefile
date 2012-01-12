GETPAR = ../../Getpar/getpar
GETPARLIB = ${GETPAR}/lib
GETPARINCL = ${GETPAR}/include


# General compiler/linker flags
CC=gcc
FC=gfortran
CFLAGS = -I${GETPARINCL} -Wall -O3 \
	-D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE \
	-D_FILE_OFFSET_BITS=64 -DUSE_IOBUF_MACROS
LIBS = -L${GETPARLIB} -lget -lgfortran -lm
LDFLAGS = -static

TARGET = librupgen.a

all: $(TARGET)

############################################
# Executables
############################################

librupgen.a: rupgen_api.o genslip.o fourg.o slip.o ruptime.o
	$(AR) rcs $@ $^


genslip: genslip.o
	$(CC) -o $@ $^ $(LIBS)

############################################
# Special build targets
############################################

clean:
	rm -f core *.o *~ $(TARGET)
