##
# @section LICENSE
# Copyright (c) 2013-2016, Regents of the University of California
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are pe
# rmitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of
# conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list
# of conditions and the following disclaimer in the documentation and/or other materials pr
# ovided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXP
# RESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERC
# HANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE CO
# PYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, E
# XEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTIT
# UTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER C
# AUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INC
# LUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
# IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
MPI_HOME=$(dir $(shell which mpicc))../

CC      = mpic++
LD      = mpic++
LDFLAGS = -Wl,-Bsymbolic-functions -Wl,-z,relro -O3 -g -D__HIP -lamdhip64 --offload-arch=gfx90a
CFLAGS = -O3 -ggdb -std=c++17 -Wall -D__HIP -ffast-math -D__HIP_ROCclr__ -D__HIP_ARCH_GFX90A__=1 -x hip --offload-arch=gfx90a 
GFLAGS = hipcc -ffast-math -O3 -D__HIP  -D__HIP_ROCclr__ -D__HIP_ARCH_GFX90A__=1 --amdgpu-target=gfx90a

INCDIR  = -I$(MPI_HOME)/include 
OBJECTS = command.o pmcl3d.o grid.o source.o mesh.o cerjan.o swap.o io.o sgt.o kernel.o
LIB     = -L$(MPI_HOME)/lib -L${ROCM_PATH}/lib -lamdhip64

pmcl3d:	$(OBJECTS)
	$(LD) $(LDFLAGs) $(INCDIR) -o	pmcl3d	$(OBJECTS)	$(LIB)
	cp pmcl3d ../bin/

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
