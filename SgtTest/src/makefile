TARGETS = dump_rawsgt check_for_nans
COBJS = sacio.o iofunc.o fft1d.o
FOBJS = fourg.o

GETPAR = ../../Getpar/getpar/lib

LIBS = -lm ${GETPAR}/libget.a -I../../Getpar/getpar/include
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

CFLAGS = ${UFLAGS} ${LF_FLAGS} ${LIBS}
FFLAGS = ${UFLAGS}

##### make options

all: $(TARGETS)
	for TARGET in $(TARGETS); do \
		rm -rf ../bin/$$TARGET ; \
		cp $$TARGET ../bin/ ; \
	done 

dump_rawsgt : dump_rawsgt.c ${COBJS} ${FOBJS}
	${CC} ${CFLAGS} -o dump_rawsgt dump_rawsgt.c ${LDLIBS}

check_for_nans : check_for_nans.c
	${CC} ${CFLAGS} -o check_for_nans check_for_nans.c -lm

clean :
	-rm -f ${COBJS} ${FOBJS}


