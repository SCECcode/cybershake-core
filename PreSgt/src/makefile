COBJS = gen_sgtgrid.o sacio.o iofunc.o fft1d.o geoproj_subs.o
FOBJS = fourg.o

GETPAR = ../../Getpar/getpar/lib

LIBS = -lm ${GETPAR}/libget.a

LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

UFLAGS = -O3

include ../../Compilers.mk
CC = $(MY_MPICC)
FC = $(MY_FC)

BIN_LIST = gen_sgtgrid

##### compile options

CFLAGS = $(UFLAGS) $(LF_FLAGS) $(MY_CFLAGS)
FFLAGS = $(UFLAGS) $(MY_FFLAGS)

all: $(BIN_LIST)
	mkdir -p ../bin
	for TARGET in $(BIN_LIST); do \
		rm -rf ../bin/$$TARGET ; \
		cp $$TARGET ../bin/ ; \
	done

gen_sgtgrid: ${COBJS} ${FOBJS}
	${CC} -o $@ $^ ${LIBS}

clean :
	-rm -f ${COBJS} ${FOBJS} $(BIN_LIST)

