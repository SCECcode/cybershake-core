# This makefile fragment helps us to choose the appropriate
# compilers for the site where we are compiling.

# Notes On Specific Modules:
#
# SpectralAcceleration/p2utils:
#       - attempts to use ifort/icc is available but will
#         fall back to the configured MY_FC declared here
# V4-WrapC/src:
#       - ensure it uses a fortran-77 compliant compiler (MY_FC77)
#       - will generate errors with gfortran 4.2 or lower
#       - may work with gfortran 4.3 and above since that supports
#         -finit-local-zero option.
#

# Get the hostname we are running on
HOSTNAME = $(shell hostname -f)

ifeq (summit,$(findstring summit, $(HOSTNAME)))
	MY_CC = xlc
	MY_FC = xlf
	MY_MPICC = mpicc
	MY_FC77 = xlf
	MY_MPIFC = mpif77
	FLAG = 1
endif


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
        MY_MPICC = mpicc
        MY_FC77 = ftn
        MY_MPIFC = ftn 
        MY_CFLAGS = 
        MY_FFLAGS = -ffixed-line-length-132
	FLAG = 1
endif

# NICS kraken  (Cray XT5)
# Note: For this to work you need to have your environment set
#       up with the gcc compilers, not the PG compilers. On
#       Kraken you need to use these commands:
#          module purge
#	   module load Base-opts
#	   module load PrgEnv-gnu
ifeq (kraken,$(findstring kraken, $(HOSTNAME)))
        MY_CC = cc
        MY_FC = ftn
        MY_FC77 = ftn
        MY_MPICC = cc
        MY_MPIFC = ftn
        MY_CFLAGS =
        MY_FFLAGS = -ffixed-line-length-132
	FLAG = 1
endif

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
ifeq ($(FLAG), 0)
        MY_CC = gcc
        MY_FC = g77
        MY_FC77 = g77
        MY_MPICC = mpicc
        MY_MPIFC = mpif77
        MY_CFLAGS = 
        MY_FFLAGS = -ffixed-line-length-132
	#MY_FFLAGS = -finit-local-zero -ffixed-line-length-132
endif

