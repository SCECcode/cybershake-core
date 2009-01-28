# This makefile fragment helps us to choose the appropriate
# compilers for the site where we are compiling.

# Get the hostname we are running on
HOSTNAME = $(shell hostname -f)

# PSC BigBen  (Cray XT3)
# Note: For this to work you need to have your environment set
#       up with the gcc compilers, not the PG compilers. On
#       BigBen you need to use these commands:
#          module switch PrgEnv-pgi PrgEnv-gnu
#          module unload acml/3.0
#          module load gcc/4.0.2
ifeq (bigben,$(findstring bigben, $(HOSTNAME)))
        MY_CC = cc
        MY_FC = ftn
        MY_MPICC = cc
        MY_MPIFC = ftn 
	MY_CFLAGS = 
	MY_FFLAGS = -ffixed-line-length-132

# Default (gcc)
# Note: For this to work you need to make sure that your
#       environment is set up to use the version of mpicc
#       and mpif77 that was configured to use the GNU 
#       compilers gcc and g77. You can accomplish this on
#       most of the teragrid sites using SoftEnv. Try typing
#       'softenv' at the command prompt to get a list of 
#       available modules. Look for one that says '+mpich-*-gcc' 
#       or something similar to that. Of course, it is going
#       to work best if the mpi you choose uses the interconnect
#       supported by the site (infiniband, myrinet, etc.).
else
        MY_CC = gcc
        MY_FC = g77
#        MY_FC = gfortran
        MY_MPICC = mpicc
        MY_MPIFC = mpif77
	MY_CFLAGS = 
	MY_FFLAGS = -ffixed-line-length-132
endif

