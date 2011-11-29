HEADS = include.h structure.h function.h defs.h
OBJS = fourg.o iofunc.o misc.o slip.o ruptime.o srf_subs.o

#GETPARLIB = ${HOME}/Getpar/getpar/lib
GETPARLIB = /home/rcf-104/CyberShake2007/software/Getpar/getpar/lib

GCCLIB = /usr/lib/gcc/x86_64-redhat-linux/4.3.2
FORTLIBS = -lgfortranbegin -lgfortran -lgcc_s
FC = gfortran

LDLIBS = ${OBJS} -L${GCCLIB} -L${GETPARLIB} -lm -lget ${FORTLIBS}

#LF_FLAGS = -D_FILE_OFFSET_BITS=32
#
# use following for large file capability
LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

UFLAGS = -O3

CC = gcc

CFLAGS = ${UFLAGS} ${LF_FLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132

##### make options

generf : generf.o ${OBJS}
	g77 -o generf generf.o ${LDLIBS}

fd_ruptime : fd_ruptime.c ${OBJS}
	${CC} ${CFLAGS} -o fd_ruptime fd_ruptime.c ${LDLIBS}

genslip-v3.0.1 : genslip-v3.0.1.c ${OBJS}
	${CC} -o genslip-v3.0.1 genslip-v3.0.1.c ${LDLIBS}

genslip-v3.2 : genslip-v3.2.c ${OBJS}
	${CC} -o genslip-v3.2 genslip-v3.2.c ${LDLIBS}

genslip-v3.1 : genslip-v3.1.c ${OBJS}
	${CC} -o genslip-v3.1 genslip-v3.1.c ${LDLIBS}

genslip-v3.0 : genslip-v3.0.c ${OBJS}
	${CC} -o genslip-v3.0 genslip-v3.0.c ${LDLIBS}

genslip-v2.4 : genslip-v2.4.o ${OBJS}
	g77 -o genslip-v2.4 genslip-v2.4.o ${LDLIBS}

genslip-v2.3 : genslip-v2.3.o ${OBJS}
	g77 -o genslip-v2.3 genslip-v2.3.o ${LDLIBS}

genslip-v2.2 : genslip-v2.2.o ${OBJS}
	g77 -o genslip-v2.2 genslip-v2.2.o ${LDLIBS}

genslip-v2.1 : genslip-v2.1.o ${OBJS}
	g77 -o genslip-v2.1 genslip-v2.1.o ${LDLIBS}

genslip-v2.0 : genslip-v2.0.o ${OBJS}
	g77 -o genslip-v2.0 genslip-v2.0.o ${LDLIBS}

genslip-mreal : genslip-mreal.o ${OBJS}
	g77 -o genslip-mreal genslip-mreal.o ${LDLIBS}

genslip-mreal_OLD-IO : genslip-mreal_OLD-IO.o ${OBJS}
	g77 -o genslip-mreal_OLD-IO genslip-mreal_OLD-IO.o ${LDLIBS}

genslip : genslip.o ${OBJS}
	g77 -o genslip genslip.o ${LDLIBS}

${OBJS} : ${HEADS}

clean :
	-rm -f ${OBJS} genslip-v3.2
