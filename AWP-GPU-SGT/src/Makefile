CC 	= cc
CFLAGS	= -O3 -Wall
GFLAGS	= nvcc -use_fast_math -arch=sm_35

INCDIR  = -I/opt/nvidia/cudatoolkit/5.0.35.102/include #-I/usr/include/mpich2/
OBJECTS	= command.o pmcl3d.o grid.o source.o mesh.o cerjan.o swap.o kernel.o io.o sgt.o
LIB	= -lm -ldl -L/opt/nvidia/cudatoolkit/5.0.35.102/lib64 -lcudart -lmpich

pmcl3d:	$(OBJECTS)
	$(CC) $(CFLAGS) $(INCDIR) -o	pmcl3d	$(OBJECTS)	$(LIB)

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

clean:	
	rm *.o
