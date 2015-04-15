HEADS = include.h structure.h function.h defs.h
OBJS = fourg.o iofunc.o misc.o slip.o ruptime.o srf_subs.o

#GETPAR = ${HOME}/Getpar/getpar/lib
GETPAR = /home/rcf-104/CyberShake2007/software/Getpar/getpar/lib

LIBS = -lm ${GETPAR}/libget.a
LDLIBS = ${OBJS} ${LIBS}

#LF_FLAGS = -D_FILE_OFFSET_BITS=32
#
# use following for large file capability
LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

UFLAGS = -O3

CC = gcc
FC = g77

CFLAGS = ${UFLAGS} ${LF_FLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132

##### make options

genslip-v2.2 : genslip-v2.2.c ${OBJS}
	${CC} ${CFLAGS} -o genslip-v2.2 genslip-v2.2.c ${LDLIBS}

${OBJS} : ${HEADS}

clean :
	-rm -f ${OBJS} genslip-v2.2
