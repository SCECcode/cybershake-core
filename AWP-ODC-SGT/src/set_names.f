c---------------------------------------------
c set_names():  Retrieve command line arguments (CLA) if provided.
c	        Otherwise default to original names with path in IN3D.
c Command line names in this order:
c	#1  IN3D (contains all other names)
c	#2  source
c	#3  media
c	#4  CHKP
c	#5  CHKJOB
c	#6  SSX3D
c	#7  SSY3D
c	#8  SSZ3D
c	#9  SSX3D2
c	#10 SSY3D2
c	#11 SSZ3D2
c---------------------------------------------
	subroutine set_names()

	include 'set_names.h'

        integer ::      nargs

	num_commandlineargs = 0

c retrieve command line user options.  

	c_IN3D    = 	'IN3D'

	nargs = iargc()
	if(nargs .EQ. 0)then
		c_IN3D = 	'IN3D'
	endif

		
	if(nargs .GE.  1)call getarg( 1,c_IN3D)
	if(nargs .GE.  2)call getarg( 2,c_SOURCE)
	if(nargs .GE.  3)call getarg( 3,c_MEDIA)
	if(nargs .GE.  4)call getarg( 4,c_CHKMOD)
	if(nargs .GE.  5)call getarg( 5,c_CHKJOB)
	if(nargs .GE.  6)call getarg( 6,c_SSX3D)
	if(nargs .GE.  7)call getarg( 7,c_SSY3D)
	if(nargs .GE.  8)call getarg( 8,c_SSZ3D)
	if(nargs .GE.  9)call getarg( 9,c_SSX3D2)
	if(nargs .GE. 10)call getarg(10,c_SSY3D2)
	if(nargs .GE. 11)call getarg(11,c_SSZ3D2)

	num_commandlineargs = nargs
   
	return
	end
