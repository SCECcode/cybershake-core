CUDA_HOME = $(subst bin/nvcc,,$(shell which nvcc))
MPI_HOME=$(dir $(shell which mpicc))../

CC      = mpicc
CFLAGS = -O3 -g
#CFLAGS  = -g -O0
#GFLAGS  = $(CUDA_HOME)bin/nvcc -use_fast_math -arch=sm_70 -O3 -lineinfo -Xptxas -v
GFLAGS = $(CUDA_HOME)bin/nvcc -O0

INCDIR  = -I$(CUDA_HOME)include -I$(MPI_HOME)/include
OBJECTS = command.o pmcl3d.o grid.o source.o mesh.o cerjan.o swap.o kernel.o io.o sgt.o
OLD_OBJECTS = command.o pmcl3d.o grid.o source.o mesh.o cerjan.o swap.o kernel.old.o io.o sgt.o
LIB     = -L$(CUDA_HOME)/lib64 -L$(MPI_HOME)/lib -ldl -lcudart -lmpi_ibm -lm -lpthread -lnvToolsExt -lstdc++
#LIB = -L$(CUDA_HOME)/lib64 -lcudart -lm

pmcl3d:	$(OBJECTS)
	$(CC) $(CFLAGS) $(INCDIR) -o	pmcl3d	$(OBJECTS)	$(LIB)
	cp pmcl3d ../bin/

pmcl3d_old: $(OLD_OBJECTS)
	$(CC) $(CFLAGS) $(INCDIR) -o    pmcl3d_old $(OLD_OBJECTS) $(LIB)
	cp pmcl3d_old ../bin/

pmcl3d.o:	pmcl3d.c
	$(CC) $(CFLAGS) $(INCDIR) -c -o pmcl3d.o	pmcl3d.c		

command.o:	command.c
	$(CC) $(CFLAGS) $(INCDIR) -c -o	command.o	command.c	

io.o:	  io.c
	$(CC) $(CFLAGS) $(INCDIR) -c -o	io.o	  io.c	

sgt.o:		sgt.c
	$(CC) $(CFLAGS) $(INCDIR) -c -o sgt.o		sgt.c		

grid.o:		grid.c
	$(CC) $(CFLAGS) $(INCDIR) -c -o grid.o		grid.c		

source.o:	source.c
	$(CC) $(CFLAGS) $(INCDIR) -c -o source.o	source.c	

mesh.o:		mesh.c
	$(CC) $(CFLAGS) $(INCDIR) -c -o mesh.o		mesh.c		

cerjan.o:	cerjan.c
	$(CC) $(CFLAGS) $(INCDIR) -c -o cerjan.o	cerjan.c

swap.o:		swap.c
	$(CC) $(CFLAGS) $(INCDIR) -c -o swap.o		swap.c

kernel.o:	kernel.cu
	$(GFLAGS) $(INCDIR) -c -o	kernel.o	kernel.cu	

kernel.old.o: kernel.old.cu
	$(GFLAGS) $(INCDIR) -c -o   kernel.old.o  kernel.old.cu

clean:	
	rm *.o
