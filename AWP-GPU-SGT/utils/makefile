TARGETS = reformat_gpu filter_sgt interpolate_momrate

include ../../Compilers.mk

CC = $(MY_CC)
FC = $(MY_FC)

BBP_DIR = /gpfs/alpine/proj-shared/geo112/CyberShake/software/bbp/bbp-19.8.0-python3/bbp
GP_DIR = ${BBP_DIR}/src/gp
WCC_DIR = ${GP_DIR}/WccFormat/Progs
BBP_OBJS = ${WCC_DIR}/wcc_tfilter_sub.o ${WCC_DIR}/iofunc.o

GETPAR_DIR = /gpfs/alpine/proj-shared/geo112/CyberShake/software/Getpar
GETPAR_LIB = -L${GETPAR_DIR}/getpar/lib -lget
GETPAR_INC = -I${GETPAR_DIR}/getpar/include
BBP_INC = -I${WCC_DIR}


reformat_gpu : reformat_gpu.c
	${CC} -o reformat_gpu reformat_gpu.c

filter_sgt : filter_sgt.c ${BBP_OBJS}
	$(CC) -o filter_sgt filter_sgt.c ${BBP_OBJS} ${BBP_INC} ${GETPAR_INC} ${GETPAR_LIB}

filter_vel : filter_vel.c ${BBP_OBJS}
	$(CC) -o filter_vel filter_vel.c ${BBP_OBJS} ${BBP_INC} ${GETPAR_INC} ${GETPAR_LIB}

interpolate_momrate : interpolate_momrate.c 
	$(CC) -o interpolate_momrate interpolate_momrate.c -lm
