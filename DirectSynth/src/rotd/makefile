FC=gfortran
FFLAGS = -O0 -ffixed-line-length-none
CFLAGS = -O0
HEADS = baseline.h
OBJS = calcrsp.o fftsub.o ft_th.o sort.o spline.o splint.o
TARGETS = rotd

all: $(TARGETS)
	for TARGET in $(TARGETS); do \
		rm -f ../bin/$$TARGET ; \
		cp $$TARGET ../bin/ ; \
	done

rotd: ${OBJS} rotd100.o rotd.o
	${FC} ${CFLAGS} -o rotd rotd.o rotd100.o ${OBJS}

rotd50: ${OBJS} rotd50.f
	${FC} ${FFLAGS} -o rotd50 rotd50.f ${OBJS}

rotd100: ${OBJS} rotd100.f
	${FC} ${FFLAGS} -o rotd100 rotd100.f ${OBJS}

${OBJS}: ${HEADS}

clean:
	rm -f *.o rotd50 rotd100
