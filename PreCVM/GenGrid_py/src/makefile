TARGETS = gen_model_cords
COBJS = geoproj_subs.o
FOBJS = geo_utm.o

GETPAR = ../../../Getpar/getpar/lib

LIBS = -lm $(GETPAR)/libget.a

LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

include ../../../Compilers.mk
CC = $(MY_CC)
FC = $(MY_FC)

#CC = gcc
#FC = f77

BIN_LIST = faultcords statgrid2latlon fd2close-dist latlon2distaz \
           latlon2xyz stat_prof_cords fdcords2latlon ll2delaz \
	              xy2ll gen_latlon_cords latlon2statgrid ll2xy xyz2lonlat \
		                 gen_model_cords latlon2topcenter llmask gcproj

##### compile options

CFLAGS = -O3 $(LF_FLAGS) $(MY_CFLAGS)
FFLAGS = -O3 $(MY_FFLAGS)

##### make options

all: $(TARGETS)
	for TARGET in $(TARGETS); do \
		rm -rf ../bin/$$TARGET ; \
		cp $$TARGET ../bin ; \
	done

# RANGER CHANGE
gen_model_cords: gen_model_cords.o $(COBJS) $(FOBJS)
	$(CC) -o $@ $^ $(LIBS)

clean:
	rm -rf *.o $(TARGETS)
