c Common block for command line definition of Input/Output names
c
c
c jan04	  D. Okaya   Added command line argument of directory where
c			in/out files are to be found.  This allows
c			code to be run from one place (e.g., home
c			directory) whereas inputs are in a scratch
c			directory elsewhere.
c 06aug04 D. Okaya   Adjusted command line arguments for fink.  Six
c			separate CLA's, one for each input and
c			output file names (per Vipin Gupta's request):
c			IN3D,SOURCE,MEDIA,SSX3D,SSY3D,SSZ3D.
c			Paths and names can actually be different.
c			If no CLA's, default to KBO original names in
c			local directory (no path).
c			CHKMOD,CHKJOB,VX,VY,VZ remain in local dir.
c 28jun06 D. Okaya   Added second set of output seismogram files for
c			Terashake-2 (4D volumes); first set is for
c			MapView timeslices.  Removed VX,VY,VZ.
c			Add SETTING file.
c 11mar09 Kwangyoon Lee    Remove Setting file 
c
c---------------------------------------------

	common /io_filenames/	
     +			num_commandlineargs,
     +			c_IN3D,c_SOURCE,c_MEDIA,c_CHKMOD,c_CHKJOB,
     +			c_SSX3D,c_SSY3D,c_SSZ3D,c_SSX3D2,c_SSY3D2,c_SSZ3D2

	integer		num_commandlineargs
	character*180	c_IN3D,c_SOURCE,c_MEDIA,c_CHKMOD,c_CHKJOB,
     +			c_SSX3D,c_SSY3D,c_SSZ3D,c_SSX3D2,c_SSY3D2,c_SSZ3D2
