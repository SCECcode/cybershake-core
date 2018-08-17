COBJS = gen_sgtgrid.o sacio.o iofunc.o fft1d.o geoproj_subs.o
FOBJS = fourg.o

GETPAR = ../../Getpar/getpar

LIBS = -lm -L${GETPAR}/lib -lget -L../../../utils/libcfu_0.03/lib -lcfu

LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

I_FLAGS = -I${GETPAR}/include/libget.h -I../../../utils/libcfu_0.03/include

UFLAGS = -O3

include ../../Compilers.mk
CC = $(MY_MPICC)
FC = $(MY_FC)

BIN_LIST = gen_sgtgrid gen_sgtgrid_qsort

##### compile options

CFLAGS = $(UFLAGS) $(I_FLAGS) $(LF_FLAGS) $(MY_CFLAGS)
FFLAGS = $(UFLAGS) $(MY_FFLAGS)

all: $(BIN_LIST)
	mkdir -p ../bin
	for TARGET in $(BIN_LIST); do \
		rm -rf ../bin/$$TARGET ; \
		cp $$TARGET ../bin/ ; \
	done

gen_sgtgrid: ${COBJS} ${FOBJS}
	${CC} ${CFLAGS} -o $@ $^ ${LIBS}

gen_sgtgrid_qsort: gen_sgtgrid_qsort.o sacio.o iofunc.o fft1d.o geoproj_subs.o
	${CC} ${CFLAGS} -o $@ $^ ${LIBS}

clean :
	-rm -f ${COBJS} ${FOBJS} $(BIN_LIST)

