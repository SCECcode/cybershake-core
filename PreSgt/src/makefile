COBJS = sacio.o iofunc.o fft1d.o geoproj_subs.o
FOBJS = fourg.o

#GETPAR = ${HOME}/Getpar/getpar/lib
GETPAR = ../../Getpar/getpar/lib

LIBS = -lm ${GETPAR}/libget.a
LDLIBS = ${COBJS} ${FOBJS} ${LIBS}

#LF_FLAGS = -D_FILE_OFFSET_BITS=32
#
# use following for large file capability
LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

UFLAGS = -O3

CC = gcc

BIN_LIST = gen_sgtgrid

##### compile options

CFLAGS = ${UFLAGS} ${LF_FLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132

#TARG = $@

##### make options

#MAIN_TARG := ${MAKECMDGOALS}
#TARG_SRC := ${MAKECMDGOALS}.c

all: $(BIN_LIST)
	for TARGET in $(BIN_LIST); do \
		rm -f ../bin/$$TARGET ; \
		cp $$TARGET ../bin ; \
	done

%: %.c

gen_sgtgrid : ${COBJS} ${FOBJS}
	${CC} -o $@ $@.c ${LDLIBS}

${COBJS} ${FOBJS} :

clean :
	-rm -f ${COBJS} ${FOBJS} $(BIN_LIST)
	make all


