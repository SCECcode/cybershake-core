COBJS = geoproj_subs.o
FOBJS = geo_utm.o

GETPAR = ../../../Getpar/getpar/lib

LIBS = -lm ${GETPAR}/libget.a
LDLIBS = ${COBJS} ${FOBJS} ${LIBS}

#LF_FLAGS = -D_FILE_OFFSET_BITS=32
#
# use following for large file capability
LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

UFLAGS = -O3

include ../../../Compilers.mk
CC = $(MY_CC)
FC = $(MY_FC)

CC = gcc

BIN_LIST = gcproj
#faultcords statgrid2latlon fd2close-dist latlon2distaz \
#           latlon2xyz stat_prof_cords fdcords2latlon ll2delaz \
#	   xy2ll gen_latlon_cords latlon2statgrid ll2xy xyz2lonlat \
#           gen_model_cords latlon2topcenter llmask gcproj fault_seg2lonlat

##### compile options

CFLAGS = ${UFLAGS} ${LF_FLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132

##### make options

all : $(BIN_LIST)

gcproj : gcproj.o gcproj.c ${COBJS} ${FOBJS}
	${CC} -o gcproj gcproj.o ${LDLIBS}
	cp $(BIN_LIST) ../bin/

#$(BIN_LIST) : $$@.o $$@.c ${COBJS} ${FOBJS}
#	${CC} -o $@ $@.o ${LDLIBS}
#	cp $(BIN_LIST) ../bin/
 
${COBJS} ${FOBJS} : 

clean:
	rm -rf *.o $(BIN_LIST)

