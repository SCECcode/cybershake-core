COBJS = sacio.o iofunc.o fft1d.o
FOBJS = fourg.o

GETPAR = ../../Getpar/getpar/lib

LIBS = -lm ${GETPAR}/libget.a
LDLIBS = ${COBJS} ${FOBJS} ${LIBS}

#LF_FLAGS = -D_FILE_OFFSET_BITS=32
#
# use following for large file capability
LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

UFLAGS = -O3

include ../../Compilers.mk
CC = $(MY_CC)
FC = $(MY_FC)

#CC = gcc
#FC = g77

##### compile options

CFLAGS = ${UFLAGS} ${LF_FLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132

##### make options

all: dump_rawsgt
	-rm -rf ../bin/dump_rawsgt ; \
	cp dump_rawsgt ../bin 

dump_rawsgt : dump_rawsgt.c ${COBJS} ${FOBJS}
	${CC} ${CFLAGS} -o dump_rawsgt dump_rawsgt.c ${LDLIBS}

clean :
	-rm -f ${COBJS} ${FOBJS}

${COBJS} ${FOBJS} :
