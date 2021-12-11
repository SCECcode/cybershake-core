TARGETS = merge_psa
HEADS = include.h structure.h function.h
OBJS =  iofunc.o misc.o wcc_tfilter.o wcc_resamp_arbdt.o wcc_add.o fourg.o surfseis_rspectra.o setparm180.o sub_bandpass.o merge.o
ROTD_OBJS = rotd/calcrsp.o rotd/fftsub.o rotd/ft_th.o rotd/sort.o rotd/spline.o rotd/splint.o rotd/rotd.o rotd/rotd100.o

LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

GETPAR = ../../Getpar/getpar
LDLIBS = -lm ${GETPAR}/lib/libget.a
IFLAGS = -I${GETPAR}/include

UFLAGS = -g

include ../../Compilers.mk

CFLAGS = ${UFLAGS} ${LF_FLAGS} ${IFLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132 -finit-local-zero

#CC = $(MY_CC)
#FC = $(MY_FC)
#CC = icc
FC = cc
CC = ftn
#FC = f77

##### make options

all: $(TARGETS)
	for TARGET in $(TARGETS); do \
		rm -f ../bin/$$TARGET ; \
		cp $$TARGET ../bin ; \
	done

merge_psa : merge_psa.o ${OBJS} ${ROTD_OBJS}
	${FC} ${LDFLAGS} -o merge_psa merge_psa.o ${OBJS} ${ROTD_OBJS} ${LDLIBS}

ray_stimes : ray_stimes.o ${OBJS}
	${CC} ${LDFLAGS} -o ray_stimes ray_stimes.o ${OBJS} ${LDLIBS}

jbsim3d : jbsim3d.o ${OBJS} ${SRF_OBJS}
	${CC} ${LDFLAGS} -o jbsim3d jbsim3d.o ${OBJS} ${LDLIBS} ${SRF_OBJS}

gen_gflist : gen_gflist.o ${OBJS} ${SRF_OBJS}
	${CC} ${LDFLAGS} -o gen_gflist gen_gflist.o ${OBJS} ${LDLIBS} ${SRF_OBJS}

bailey2srf : bailey2srf.o ${OBJS}
	${CC} ${LDFLAGS} -o bailey2srf bailey2srf.o ${OBJS} ${LDLIBS}

PARM    = setparm180.f

surfseis_rspectra : sub_bandpass.o surfseis_rspectra.o
	$(FC) $(FFLAGS) -o $@ $^ $(PARM)

${OBJS} : ${HEADS}

clean :
	-rm -f $(OBJS) $(ROTD_OBJS) $(TARGETS) *.o
