c Fortran command line/parameter file parameter retrieval system using
c	name=value pairs.
c 
c Adapted from C-based "getpar()" routines by Caltech/Stanford.
c
c
c Nomenclature for command-line argument retrieval:
c "Token"	is "name=value" string (must have '=' delimiter) [185 char] 245
c "name"	is string which will be searched by getpar()	 [ 64 char]  64
c "value"	is string containing value of 'name'		 [120 char] 180
c
c Nomenclature for  "getpar("name","type",variable):
c "name" 	is character string for searching		 [ 64 char]  64
c "type"	is defined by C-based getpar:  d,f,s, dv, fv	 [  2 char]   2
c  variable	is recipient of value (number or string)
c----------------------
c setparm_num 	= number of name=value tokens.
c		= 0 set in setpar(); turns off all subsequent getpar() routines 
c-----------------------------------------------------------------------
c ORIGINAL PURE-FORTRAN GETPAR()
c 27jul00 DAO	Initial implementation.
c 30jul00 DAO	Further implementation and restructuring.  Concatenate
c		getpar_copt.f routines.
c 31jul00 DAO   Fix the termination of char string.
c
c 19aug03 DAO   g77 version: remove dynamic memory allocation and pointers.
c			Hardwire dimensions.
c 09oct03 DAO   g77 version: fix string termination from NULL to space.
c		    Solaris allows text to be terminated with null
c			within middle of declared char size (trailing chars
c			either don't matter or are spaces).
c		    Linux seems to not deal with the null, but wants the
c			full pad using spaces. 
c		    In this version all returned strings are padded to
c			full incoming length. 
c 02nov03 DAO	g77 version: replace char*1 definition of cname(64,500) and
c			cvalue(120,500) with explicit character lengths of
c			char*64 cname(500), char*120 cvalue(500).  The 
c			string terminations using the char*1 approach does not
c			work under Linux.  This use of hardcoded string lengths
c			trickles into other subroutines.  
c			   However, this approach may remove the NULL/SPACE issue
c			which arises (see 09Oct03).  
c		Hardcoded number of pairs allows for mem_address'ing
c			of cname/cvalue arrays to be removed.
c		    Returned value(s) still need a memory pointer because we 
c			don't know its length ahead of time.
c----------------------
c GETPARMS()
c 27nov03 DAO	modify generalized getpar() to be explicit getparmX()
c			where X indicates type of requested parameter.
c			E.g., getparmi(), getparmf(), getparms().
c			Internal value passing is easier than for getpar().
c----------------------
c SETPARMS()
c 04dec03 DAO	rename functions in order to separate from getpar().
c			loadparm(),endparm(),
c			setparmi(),  setparmf(),  setparmd(), setparms(),
c			setparmiv(), setparmfv(), setparmdv().
c		Make parameter file option functional.
c 16dec03 DAO	Rename subsidiary subroutines in order to avoid any future
c			naming conflicts by other programmers.
c		Install front-to-back or back-to-front parameter searching
c			(first mention or last mention).
c 13jan04 DAO	Fix setparmd() retrieval whose integral and remainder
c			components were getting truncated to 10 digits
c			via overuse of setparm_ConstructNumber() (which
c			was meant to retrieve 32-bit integers.  This
c			required a parallel setparm_ConstructNumber8() which
c			internally uses R*8.
c
c------------------------------------------------------------------------------
c 28feb04 DAO	Make long_line version : 64+120 = 185 -> 64+180=245
c		Increase number args from 500 to 1000
c=============================================================================

c------------------------------------------------------------------------------
c loadparm() retrieves all command line arguments, parses each into name and 
c	value strings, and stores for subsequent interpretation by setparm().
c
	subroutine loadparm(cdirection)
	implicit none
	character*8	cdirection

c----------------------------
c set GLOBALS down below right away.
c	parameter (MAXCHAR_CTOKEN=245,MAXCHAR_CNAME=64,MAXCHAR_CVALUE=180)
c	parameter (MAXARGS = 1000, M_CNAME=64, M_CVALUE=180)

c global/common variables
c  MAXCHAR_CTOKEN=245, MAXCHAR_CNAME=64, MAXCHAR_CVALUE=180
	common /globals/setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE
	integer*4	setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE

	common /cnamevalue/ MAXARGS,M_CNAME,M_CVALUE,cname,cvalue
	integer		MAXARGS,M_CNAME,M_CVALUE
	character*64 	cname(1000)
	character*180	cvalue(1000)

	common /loadparm_dir/loadparm_direction
	integer*4	loadparm_direction

c internal working variables
	
	integer numargs,iargc,itoken
	character*245	ctoken245

c------------------------------------------------------
c set GLOBALS immediately
c	parameter (MAXCHAR_CTOKEN=245,MAXCHAR_CNAME=64,MAXCHAR_CVALUE=180)
c	parameter (MAXARGS = 1000, M_CNAME=64, M_CVALUE=180)
	MAXCHAR_CTOKEN=245
	MAXCHAR_CNAME=64
	MAXCHAR_CVALUE=180
	MAXARGS = 1000
	M_CNAME=MAXCHAR_CNAME
	M_CVALUE=MAXCHAR_CVALUE
	
c------------------------------------------------------
c get number of command line arguments (tokens)
	numargs = iargc()

	if(numargs .EQ. 0)then
		setparm_num    = 0
		return
	endif

c------------------------------------------------------
c cname, cvalue working space: initialize by filling with spaces = ascii 32
	call setparm_clear()

c------------------------------------------------------
c obtain setparm direction from loadparm passed argument
	call loadparm_setdir(cdirection)

c------------------------------------------------------
c obtain arguments in ascending order

	setparm_num = 0
	do itoken = 1,numargs
		call getarg(itoken,ctoken245)
		call loadparm_parse(ctoken245)
		call loadparm_parfile()
	enddo

c when here, now have parsed all valid tokens and stored in cname(),cvalue().


	RETURN
	end

c------------------------------------------------------------------------------
c endparm() terminates getparm().
c   when dynamic memory allocation was used, endpar() released the memory.
c   here, we do nothing since all dynamic mem alloc has been replaced with
c   fixed-sized arrays.
	subroutine endparm()

	implicit none

	return
	end

c---------------------------------------------------------------------
	subroutine setparm_clear()
	implicit none

	common /cnamevalue/ MAXARGS,M_CNAME,M_CVALUE,cname,cvalue
	integer		MAXARGS,M_CNAME,M_CVALUE
	character*64 	cname(1000)
	character*180	cvalue(1000)

	character*64	a64
	character*180	b180
	character*1	a1(64),b1(180)
	equivalence	(a1,a64),(b1,b180)
	integer	j
	
	call setparm_clearstring(a1,M_CNAME)
	call setparm_clearstring(b1,M_CVALUE)
	
	do j=1,MAXARGS
		cname(j) = a64
		cvalue(j) = b180
	enddo
		
	return
	end

c-----------------------
	subroutine setparm_clearstring(cstring,nlength)
	implicit none
	integer	    nlength
	character*1 cstring(nlength)


	character*1 cspace
	integer	    i

	cspace = char(32)

	do i=1,nlength
		cstring(i) = cspace
	enddo

	return
	end


c---------------------------------------------------------------------
c Parse a command line token "name=value" by using the '=' as a delimiter
c If '=' found, increment setparm_num and store name, value strings

	subroutine loadparm_parse(ctoken245)
	
	implicit none

	character*245 ctoken245

c  MAXCHAR_CTOKEN=245, MAXCHAR_CNAME=64, MAXCHAR_CVALUE=180
	common /globals/setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE
	integer*4	setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE

	common /cnamevalue/ MAXARGS,M_CNAME,M_CVALUE,cname,cvalue
	integer		MAXARGS,M_CNAME,M_CVALUE
	character*64 	cname(1000)
	character*180	cvalue(1000)
	
	integer	length_token,idelim,i,j
	character*1	ctoken(245)
	character*245	c245
	equivalence	(c245,ctoken)

	character*64	a64
	character*180	b180
	character*1	a1(64),b1(180)
	equivalence	(a64,a1),(b180,b1)

c transfer incoming token
	call setparm_clearstring(ctoken,MAXCHAR_CTOKEN)
	c245 = ctoken245

c get size of this token
	call setparm_denull(ctoken,MAXCHAR_CTOKEN)
	call setparm_strlen(ctoken,MAXCHAR_CTOKEN,length_token)
	if(length_token .EQ. 0)return
	
c search for delimiter
	do idelim=1,length_token
		if(ctoken(idelim) .EQ. '=')goto100
	enddo

c if here, did not find '=' and so is not a valid getparm() token.
c don't store anything; just return
	return

c-----------------
c if here, found valid '=' delimiter
100	continue

c first check if x=y structure exists.
	if(idelim .EQ. 1	   ) return
	if(idelim .EQ. length_token) return

c valid structure exists, now save
	call setparm_clearstring(a1,MAXCHAR_CNAME)
	call setparm_clearstring(b1,MAXCHAR_CVALUE)
	setparm_num = setparm_num+1

	j=idelim-1
	if(j .GT. MAXCHAR_CNAME)j=MAXCHAR_CNAME
	do i=1,j
		a1(i) = ctoken(i)
	enddo
	
	j = length_token - idelim
	if(j .GT. MAXCHAR_CVALUE)j=MAXCHAR_CVALUE
	do i=1,j
		b1(i) = ctoken(idelim + i)
	enddo

	cname(setparm_num)  = a64
	cvalue(setparm_num) = b180

	RETURN
	end


	subroutine setparm_strlen(ctoken,MAXCHAR_CTOKEN,length_token)
	
	implicit none

	integer	MAXCHAR_CTOKEN,length_token
	character*1 ctoken(MAXCHAR_CTOKEN)
	
	character*1 ctrail
	integer i

c	ctrail = char(0)
	ctrail = char(32)

	do i=MAXCHAR_CTOKEN,1,-1
		if(ctoken(i) .NE. ctrail)then
			length_token = i
			return
		endif
	enddo

c if here, a blank string
	length_token = 0

	return
	end


c--------------------------------------------
c identify first space or null, and clear remainder of array with NULLs.
	subroutine setparm_despace(carray,nchars)

	implicit none

	integer	i,j,nchars
	character*1 carray(nchars),cspace,cnull

	cspace = char(32)
	cnull  = char(0)

	do i=1,nchars
		if(carray(i) .EQ. cspace .OR. carray(i) .EQ. cnull)goto100
	enddo
	return

100	do j=i,nchars
		carray(j) = cnull
	enddo
	return

	end

	subroutine setparm_denull(carray,nchars)

	implicit none

	integer	i,j,nchars
	character*1 carray(nchars),cspace,cnull

	cspace = char(32)
	cnull  = char(0)

	do i=1,nchars
		if(carray(i) .EQ. cspace .OR. carray(i) .EQ. cnull)goto100
	enddo
	return

100	do j=i,nchars
		carray(j) = cspace
	enddo
	return

	end

c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
c setparm() scans tokens for valid match, then interprets token value into
c	integer, float, or char string.
c	setparm() returns #items in value string (0, 1, >1).
c
c allow up to 200 values per name (overkill as value length is 120 chars).
c 28feb04     200 values per name (overkill as value length is 180 chars).
c------------------------------------------------------------------------------
c********************************************
	integer function setparmi(c64,kreturn)
	
	implicit none

	character*(*)	c64
	integer		kreturn

c global/common variables
	common /globals/setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE
	integer*4	setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE

	common /cnamevalue/ MAXARGS,M_CNAME,M_CVALUE,cname,cvalue
	integer		MAXARGS,M_CNAME,M_CVALUE
	character*64 	cname(1000)
	character*180	cvalue(1000)

c for internal manipulation of incoming c64,ctype arguments
	character*64 cstring64
	character*1 cstring(64)
	equivalence (cstring64,cstring)

c internal variables for extracting value to return.
c following has size of MAXCHAR_CVALUE
	character*1 kvalue(180)
	integer itoken,nvalues,jstatus
	integer	jvalue

c-------------------------
c if no command line arguments to scan, immediately return	
	if(setparm_num .EQ. 0)then
		setparmi=0
		return
	endif

c-------------------------
c search command line tokens for this setparm("cstring", value) cstring
	cstring64 = c64
	call setparm_denull(cstring,MAXCHAR_CNAME)
	
	call setparm_whichComLineToken(cname,MAXCHAR_CNAME,setparm_num,
     +							cstring,itoken)

c if no token found, return without modifying returned value.  
	if(itoken .EQ. -1)then
		setparmi=0
		return
	endif

c transfer ComLine Token text from cvalue() to kvalue()to work with 1-D string
	call setparm_getComLineToken(cvalue,MAXCHAR_CVALUE,setparm_num,
     +								itoken,kvalue)

c-------------------------
c setparm(integer)
	call setparm_valueInteger(kvalue,MAXCHAR_CVALUE,jvalue,jstatus)
     	if(jstatus .EQ. 1)then
     		kreturn = jvalue
     		nvalues=1
     	else
     		nvalues=0
     	endif

c return
1000	continue
	setparmi=nvalues

	return
	end


c********************************************
	integer function setparmf(c64,xreturn)
	
	implicit none

	character*(*)	c64
	real		xreturn

c global/common variables
	common /globals/setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE
	integer*4	setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE

	common /cnamevalue/ MAXARGS,M_CNAME,M_CVALUE,cname,cvalue
	integer		MAXARGS,M_CNAME,M_CVALUE
	character*64 	cname(1000)
	character*180	cvalue(1000)

c for internal manipulation of incoming c64,ctype arguments
	character*64 cstring64
	character*1 cstring(64)
	equivalence (cstring64,cstring)

c internal variables for extracting value to return.
c following has size of MAXCHAR_CVALUE
	character*1 kvalue(180)
	integer itoken,nvalues,jstatus
	real	xvalue

c-------------------------
c if no command line arguments to scan, immediately return	
	if(setparm_num .EQ. 0)then
		setparmf=0
		return
	endif

c-------------------------
c search command line tokens for this setparm("cstring", value) cstring
	cstring64 = c64
	call setparm_denull(cstring,MAXCHAR_CNAME)
	
	call setparm_whichComLineToken(cname,MAXCHAR_CNAME,setparm_num,
     +							cstring,itoken)

c if no token found, return without modifying returned value.  
	if(itoken .EQ. -1)then
		setparmf=0
		return
	endif

c transfer ComLine Token text from cvalue() to kvalue()to work with 1-D string
	call setparm_getComLineToken(cvalue,MAXCHAR_CVALUE,setparm_num,
     +								itoken,kvalue)

c-------------------------
c setparm(floating point / real)
	call setparm_valueFloat(kvalue,MAXCHAR_CVALUE,xvalue,jstatus)
	if(jstatus .EQ. 1)then
		xreturn = xvalue
		nvalues=1
	else
		nvalues=0
	endif

c return
1000	continue
	setparmf=nvalues

	return
	end


c********************************************
	integer function setparmd(c64,dreturn)
	
	implicit none

	character*(*)	c64
	real*8		dreturn

c global/common variables
	common /globals/setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE
	integer*4	setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE

	common /cnamevalue/ MAXARGS,M_CNAME,M_CVALUE,cname,cvalue
	integer		MAXARGS,M_CNAME,M_CVALUE
	character*64 	cname(1000)
	character*180	cvalue(1000)

c for internal manipulation of incoming c64,ctype arguments
	character*64 cstring64
	character*1 cstring(64)
	equivalence (cstring64,cstring)

c internal variables for extracting value to return.
c following has size of MAXCHAR_CVALUE
	character*1 kvalue(180)
	integer itoken,nvalues,jstatus
	real*8	dvalue

c-------------------------
c if no command line arguments to scan, immediately return	
	if(setparm_num .EQ. 0)then
		setparmd=0
		return
	endif

c-------------------------
c search command line tokens for this setparm("cstring", value) cstring
	cstring64 = c64
	call setparm_denull(cstring,MAXCHAR_CNAME)
	
	call setparm_whichComLineToken(cname,MAXCHAR_CNAME,setparm_num,
     +							cstring,itoken)

c if no token found, return without modifying returned value.  
	if(itoken .EQ. -1)then
		setparmd=0
		return
	endif

c transfer ComLine Token text from cvalue() to kvalue()to work with 1-D string
	call setparm_getComLineToken(cvalue,MAXCHAR_CVALUE,setparm_num,
     +								itoken,kvalue)

c-------------------------
c setparm(double float)
	call setparm_valueDoubleFloat(kvalue,MAXCHAR_CVALUE,dvalue,jstatus)
	if(jstatus .EQ. 1)then
		dreturn = dvalue
		nvalues=1
	else
		nvalues=0
	endif

c return
1000	continue
	setparmd=nvalues

	return
	end



c********************************************
	integer function setparmiv(c64,kreturn,narray)
	
	implicit none

	integer		narray, kreturn(narray)
	character*(*)	c64

c global/common variables
	common /globals/setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE
	integer*4	setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE

	common /cnamevalue/ MAXARGS,M_CNAME,M_CVALUE,cname,cvalue
	integer		MAXARGS,M_CNAME,M_CVALUE
	character*64 	cname(1000)
	character*180	cvalue(1000)

c for internal manipulation of incoming c64,ctype arguments
	character*64 cstring64
	character*1 cstring(64)
	equivalence (cstring64,cstring)

c internal variables for extracting value to return.
c following has size of MAXCHAR_CVALUE
	character*1 kvalue(180)
	integer itoken,nvalues,jstatus
	integer	MAX_ARRAY
	PARAMETER (MAX_ARRAY = 200)
	integer	numDelimits,nBytesPerWord,	i,ntransfer
	integer*4 ivector4(200)

c-------------------------
c if no command line arguments to scan, immediately return	
	if(setparm_num .EQ. 0)then
		setparmiv=0
		return
	endif

c-------------------------
c search command line tokens for this setparm("cstring", value) cstring
	cstring64 = c64
	call setparm_denull(cstring,MAXCHAR_CNAME)
	
	call setparm_whichComLineToken(cname,MAXCHAR_CNAME,setparm_num,
     +							cstring,itoken)

c if no token found, return without modifying returned value.  
	if(itoken .EQ. -1)then
		setparmiv=0
		return
	endif

c transfer ComLine Token text from cvalue() to kvalue()to work with 1-D string
	call setparm_getComLineToken(cvalue,MAXCHAR_CVALUE,setparm_num,
     +								itoken,kvalue)

c-------------------------
c setparm(vector of integers)
        call setparm_NumDelimits(kvalue,MAXCHAR_CVALUE,numDelimits)
	if(numDelimits .EQ. 0)then
		nvalues=0
		goto1000
	endif
	nBytesPerWord = 4
	call setparm_valueVinteger(kvalue,MAXCHAR_CVALUE,ivector4,numDelimits,
     +								jstatus)
	if(jstatus .EQ. 1)then
		ntransfer=numDelimits
		if(ntransfer .GT. narray)ntransfer=narray
		do i=1,ntransfer
			kreturn(i) = ivector4(i)
		enddo
		nvalues = ntransfer
	else
		nvalues=0
	endif

c return
1000	continue
	setparmiv=nvalues

	return
	end



c********************************************
	integer function setparmfv(c64,xreturn,narray)
	
	implicit none

	character*(*)	c64
	integer		narray
	real		xreturn(narray)

c global/common variables
	common /globals/setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE
	integer*4	setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE

	common /cnamevalue/ MAXARGS,M_CNAME,M_CVALUE,cname,cvalue
	integer		MAXARGS,M_CNAME,M_CVALUE
	character*64 	cname(1000)
	character*180	cvalue(1000)

c for internal manipulation of incoming c64,ctype arguments
	character*64 cstring64
	character*1 cstring(64)
	equivalence (cstring64,cstring)

c internal variables for extracting value to return.
c following has size of MAXCHAR_CVALUE
	character*1 kvalue(180)
	integer itoken,nvalues,jstatus
	integer	MAX_ARRAY
	PARAMETER (MAX_ARRAY = 200)
	integer	numDelimits,nBytesPerWord,	i,ntransfer
	real	vector4(200)

c-------------------------
c if no command line arguments to scan, immediately return	
	if(setparm_num .EQ. 0)then
		setparmfv=0
		return
	endif

c-------------------------
c search command line tokens for this setparm("cstring", value) cstring
	cstring64 = c64
	call setparm_denull(cstring,MAXCHAR_CNAME)
	
	call setparm_whichComLineToken(cname,MAXCHAR_CNAME,setparm_num,
     +							cstring,itoken)

c if no token found, return without modifying returned value.  
	if(itoken .EQ. -1)then
		setparmfv=0
		return
	endif

c transfer ComLine Token text from cvalue() to kvalue()to work with 1-D string
	call setparm_getComLineToken(cvalue,MAXCHAR_CVALUE,setparm_num,
     +								itoken,kvalue)

c-------------------------
c setparm(vector of floating point/reals)
        call setparm_NumDelimits(kvalue,MAXCHAR_CVALUE,numDelimits)
	if(numDelimits .EQ. 0)then
		nvalues=0
		goto1000
	endif
	nBytesPerWord = 4
	call setparm_valueVfloat(kvalue,MAXCHAR_CVALUE,vector4,numDelimits,
     +								jstatus)
	if(jstatus .EQ. 1)then
		ntransfer=numDelimits
		if(ntransfer .GT. narray)ntransfer=narray
		do i=1,ntransfer
			xreturn(i) = vector4(i)
		enddo
		nvalues = ntransfer
	else
		nvalues=0
	endif

c return
1000	continue
	setparmfv=nvalues

	return
	end



c********************************************
	integer function setparmdv(c64,xreturn,narray)
	
	implicit none

	character*(*)	c64
	integer		narray
	real*8		xreturn(narray)

c global/common variables
	common /globals/setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE
	integer*4	setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE

	common /cnamevalue/ MAXARGS,M_CNAME,M_CVALUE,cname,cvalue
	integer		MAXARGS,M_CNAME,M_CVALUE
	character*64 	cname(1000)
	character*180	cvalue(1000)

c for internal manipulation of incoming c64,ctype arguments
	character*64 cstring64
	character*1 cstring(64)
	equivalence (cstring64,cstring)

c internal variables for extracting value to return.
c following has size of MAXCHAR_CVALUE
	character*1 kvalue(180)
	integer itoken,nvalues,jstatus
	integer	MAX_ARRAY
	PARAMETER (MAX_ARRAY = 200)
	integer	numDelimits,nBytesPerWord,	i,ntransfer
	real*8	vector8(200)

c-------------------------
c if no command line arguments to scan, immediately return	
	if(setparm_num .EQ. 0)then
		setparmdv=0
		return
	endif

c-------------------------
c search command line tokens for this setparm("cstring", value) cstring
	cstring64 = c64
	call setparm_denull(cstring,MAXCHAR_CNAME)
	
	call setparm_whichComLineToken(cname,MAXCHAR_CNAME,setparm_num,
     +							cstring,itoken)

c if no token found, return without modifying returned value.  
	if(itoken .EQ. -1)then
		setparmdv=0
		return
	endif

c transfer ComLine Token text from cvalue() to kvalue()to work with 1-D string
	call setparm_getComLineToken(cvalue,MAXCHAR_CVALUE,setparm_num,
     +								itoken,kvalue)

c-------------------------
c setparm(vector of double floats)
        call setparm_NumDelimits(kvalue,MAXCHAR_CVALUE,numDelimits)
	if(numDelimits .EQ. 0)then
		nvalues=0
		goto1000
	endif
	nBytesPerWord = 8
	call setparm_valueVdoublefloat(kvalue,MAXCHAR_CVALUE,vector8,
     +							numDelimits,jstatus)
	if(jstatus .EQ. 1)then
		ntransfer=numDelimits
		if(ntransfer .GT. narray)ntransfer=narray
		do i=1,ntransfer
			xreturn(i) = vector8(i)
		enddo
		nvalues = ntransfer
	else
		nvalues=0
	endif

c return
1000	continue
	setparmdv=nvalues

	return
	end



c********************************************
	integer function setparms(c64,creturn_ptr)
	
	implicit none

c	character*64	c64
	character*(*)	c64
	character*(*) creturn_ptr

c global/common variables
	common /globals/setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE
	integer*4	setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE

	common /cnamevalue/ MAXARGS,M_CNAME,M_CVALUE,cname,cvalue
	integer		MAXARGS,M_CNAME,M_CVALUE
	character*64 	cname(1000)
	character*180	cvalue(1000)

c for internal manipulation of incoming c64,ctype arguments
	character*64 cstring64
	character*1 cstring(64)
	equivalence (cstring64,cstring)

c internal variables for extracting value to return.
c following has size of MAXCHAR_CVALUE
	character*1 kvalue(180)
	character*180 k180
	integer itoken,nvalues,jstatus, MAX_LEN

c-------------------------
c if no command line arguments to scan, immediately return	
	if(setparm_num .EQ. 0)then
		setparms=0
		return
	endif

c-------------------------
c search command line tokens for this setparm("cstring", value) cstring
	cstring64 = c64
	call setparm_denull(cstring,MAXCHAR_CNAME)
	
	call setparm_whichComLineToken(cname,MAXCHAR_CNAME,setparm_num,
     +							cstring,itoken)

c if no token found, return without modifying returned value.  
	if(itoken .EQ. -1)then
		setparms=0
		return
	endif

c transfer ComLine Token text from cvalue() to kvalue()to work with 1-D string
c	call setparm_getComLineToken(cvalue,MAXCHAR_CVALUE,setparm_num,
c     +								itoken,kvalue)
	k180 = cvalue(itoken)

c-------------------------
c setparm(character string)
	call setparm_valueString(kvalue,MAXCHAR_CVALUE,jstatus)
	if(jstatus .EQ. 1)then
c		call setreturn_string(kvalue,MAXCHAR_CVALUE,k180)
		MAX_LEN = LEN(creturn_ptr)
		if(MAX_LEN .LE. MAXCHAR_CVALUE)then
			creturn_ptr=k180(1:MAX_LEN)
		else
			creturn_ptr=k180(1:MAXCHAR_CVALUE)
		endif
		nvalues=1
	else
		nvalues=0
	endif

c return
1000	continue
	setparms=nvalues

	return
	end

c	subroutine setreturn_string(kvalue,MAXCHAR_CVALUE,k180)
c	implicit none
c	integer	    MAXCHAR_CVALUE
c	character*1 kvalue(MAXCHAR_CVALUE),  	c1(180)
c	character*180 k180,			c180
c	integer 	i
c	equivalence (c180,c1)
c
c	do i=1,180
c		c1(i) = kvalue(i)
c	enddo
c
c	k180 = c180
c
c	return
c	end



c------------------------------------------------------------------------------
c identify which command line token belongs to this getpar(cstring) call.
c
	subroutine setparm_whichComLineToken(cname,MAXCHAR_CNAME,numargs,
     +							cstring,itoken)
	
	implicit none

	integer	MAXCHAR_CNAME,numargs,itoken
	character*64 cname(numargs), c64
	character*1 cstring(64), c1(64)
	equivalence (c64,c1)

	integer  i,j
	integer	 istart,iend,idir

	common /loadparm_dir/loadparm_direction
	integer*4	loadparm_direction

	istart = 1
	iend   = numargs
	idir   = 1
	if(loadparm_direction .EQ. -1)then
		istart = numargs
		iend   = 1
		idir   = -1
	endif

	itoken = -1
	
	do j=istart,iend,idir
		c64 = cname(j)
		do i=1,MAXCHAR_CNAME
			if(cstring(i) .NE. c1(i))goto100
		enddo
c if here, found a match
		goto 200

c if here, keep looking
100		continue
	enddo

c if here, no match; return with itoken = -1
	return

c if here, found match and return valid itoken
200	continue
	itoken = j
	
	return
	end


c transfer ComLine Token text from char*180 cvalue(setparm_num) to 
c    kvalue() in order to work within 1D character array.
	subroutine setparm_getComLineToken(cvalue,MAXCHAR_CVALUE,numargs,
     +								itoken,kvalue)

	implicit none

	integer	MAXCHAR_CVALUE,numargs,itoken
	character*180 cvalue(numargs), k180
	character*1 kvalue(180),       k1(180)
	equivalence (k180,k1)
	character*1 cspace,cnull
	integer i

	cspace = char(32)
	cnull  = char(0)

c clear returned string array
	do i=1,MAXCHAR_CVALUE
		kvalue(i) = cspace
	enddo

c transfer this token's =value string
	k180 = cvalue(itoken)
	do i=1,MAXCHAR_CVALUE
		kvalue(i) = k1(i)
	enddo

	return
	end

c----------------
c Get number of delimited values that should exist within character
c   string VALUE() based on number of delimiting characters ','.
c Incoming:
c     kvalue(MAXCHAR_CVALUE)	input text string
c Outgoing:
c     numvalues 		number of values (#delimiters + 1)
c
	subroutine setparm_NumDelimits(kvalue,MAXCHAR_CVALUE,numvalues)

	implicit none

	integer	MAXCHAR_CVALUE,numvalues
	character*1 kvalue(MAXCHAR_CVALUE)

	integer i

	numvalues = 0

	do i=1,MAXCHAR_CVALUE
		if(kvalue(i) .EQ. ',')numvalues = numvalues + 1
	enddo

	numvalues = numvalues + 1

	return
	end

c----------------
c Get location of delimiter for ivalueTH value within cvalue() string.
c Incoming:
c     kvalue(MAXCHAR_CVALUE)	input text string
c     length_token		end of text contents within kvalue()
c     ivalue			ivalueTH value to find trailign delimiter
c Outgoing:
c     ldemiit 			location in kvalue() of desired delimiter.
c
	subroutine setparm_FindDelimiter(kvalue,MAXCHAR_CVALUE,length_token,
     +							ivalue,ldelimit)

	implicit none

	integer	MAXCHAR_CVALUE,length_token,ivalue,ldelimit
	character*1 kvalue(MAXCHAR_CVALUE)

	integer i,item

	ldelimit = 0
	item = 0

	do i=1,length_token
		if(kvalue(i) .EQ. ',')then
			item = item + 1
			if(item .EQ. ivalue)then
				ldelimit = i
				return
			endif
		endif
	enddo

c last valued item will not have a trailing ',' when i reaches length_token
	if(ivalue .EQ. 1)     ldelimit = length_token + 1
	if(ivalue .EQ. item+1)ldelimit = length_token + 1
	if(ivalue .NE. item+1)ldelimit = 0

	return
	end

c------------------------------------------------------------------------------
c construct integer value going left to right. Quit at first non-number.
c incoming:   kvalue(istart:iend)= section of char text to extract from.
c outgoing:   jvalue		= extracted integer number
c	      jsign 		= sign of extracted integer number
c	      jstatus		= 4 : successful over whole range.
c				= 3 : encountered exponent character in string
c				= 2 : encountered decimal point within string
c				= 1 : encountered non-valid character in string
c				= 0 : no usable value
c	      jterminate	= when jstatus = 1, 2, or 3: location
c				    of encountered character in kvalue()
c
	subroutine setparm_ConstructNumber(kvalue,nvalue,istart,iend,
     +					jsign,jvalue,jstatus,jterminate)

	implicit none

	integer	    nvalue,istart,iend,jsign,jvalue,jstatus,jterminate
	character*1 kvalue(nvalue)

	integer	i,ksign,kstart,kend

	jsign   = 1
	jvalue  = 0
	jstatus = 0
	jterminate = 0
	kstart  = istart
	kend    = iend

c check if negative number or positive sign given.
	ksign = +1
	if(kvalue(kstart) .EQ. '-')then
		if(kstart .EQ. kend)then
			jsign = 0
			jstatus = 0
			jvalue = 0
			return
		endif
		ksign = -1
		kstart = kstart + 1
	elseif(kvalue(kstart) .EQ. '+')then
		if(kstart .EQ. kend)then
			jsign = 0
			jstatus = 0
			jvalue = 0
			return
		endif
		ksign = +1
		kstart = kstart + 1
	endif

c extract integer number
	do i=kstart,kend
c check for period '.'
		if(kvalue(i) .EQ. '.')then
			jsign  = ksign
			jvalue = jvalue
			jstatus = 2
			jterminate = i
			return
c check for exponent indicator 'e' 'E'  and 'd' 'D'
		elseif(kvalue(i) .EQ. 'e' .OR. kvalue(i) .EQ. 'E' .OR.
     +				kvalue(i) .EQ. 'd' .OR. kvalue(i) .EQ. 'D')then
			jsign  = ksign
			jvalue = jvalue
			jstatus = 3
			jterminate = i
			return
c check for non-valid character
		elseif(ichar(kvalue(i)).LT.48 .OR. ichar(kvalue(i)).GT.57)then
			jsign  = ksign 
			jvalue = jvalue
			jstatus = 1
			jterminate = i
			return
		else
c accumulate this digit
c first check if potential to overflow:
			if(jvalue .GT.200000000)then
				jsign = 0
				jvalue = 0
				jstatus = 0
				return
			endif
			jvalue = jvalue*10 + (ichar(kvalue(i)) -48)
		endif
	enddo

c if here, ran length of istart:iend range and have valid number
	jsign  = ksign 
	jvalue = jvalue
	jstatus = 4
	jterminate=kend+1

	return
	end

c..............................................................................
c construct Double value going left to right. Quit at first non-number.
c incoming:   kvalue(istart:iend)= section of char text to extract from.
c outgoing:   xvalue8		= extracted double number
c	      xsign8 		= sign of extracted double number
c	      jstatus		= 4 : successful over whole range.
c				= 3 : encountered exponent character in string
c				= 2 : encountered decimal point within string
c				= 1 : encountered non-valid character in string
c				= 0 : no usable value
c	      jterminate	= when jstatus = 1, 2, or 3: location
c				    of encountered character in kvalue()
c
	subroutine setparm_ConstructNumber8(kvalue,nvalue,istart,iend,
     +					xsign8,xvalue8,jstatus,jterminate)

	implicit none

	integer	    nvalue,istart,iend, jstatus,jterminate
	real*8	    xsign8,xvalue8
	character*1 kvalue(nvalue)

	integer	i,kstart,kend,	kwork
	real*8	ysign8,setparm_INT2DOUBLE


	xsign8  = 1.0D0
	xvalue8 = 0.0D0
	jstatus = 0
	jterminate = 0
	kstart  = istart
	kend    = iend

c check if negative number or positive sign given.
	ysign8 = +1.0D0
	if(kvalue(kstart) .EQ. '-')then
		if(kstart .EQ. kend)then
			xsign8  = 0.0D0
			xvalue8 = 0.0D0
			jstatus = 0
			return
		endif
		ysign8 = -1.0D0
		kstart = kstart + 1
	elseif(kvalue(kstart) .EQ. '+')then
		if(kstart .EQ. kend)then
			xsign8  = 0.0D0
			xvalue8 = 0.0D0
			jstatus = 0
			return
		endif
		ysign8 = +1.0D0
		kstart = kstart + 1
	endif

c extract integer number
	do i=kstart,kend
c check for period '.'
		if(kvalue(i) .EQ. '.')then
			xsign8  = ysign8
			xvalue8 = xvalue8
			jstatus = 2
			jterminate = i
			return
c check for exponent indicator 'e' 'E'  and 'd' 'D'
		elseif(kvalue(i) .EQ. 'e' .OR. kvalue(i) .EQ. 'E' .OR.
     +				kvalue(i) .EQ. 'd' .OR. kvalue(i) .EQ. 'D')then
			xsign8  = ysign8
			xvalue8 = xvalue8
			jstatus = 3
			jterminate = i
			return
c check for non-valid character
		elseif(ichar(kvalue(i)).LT.48 .OR. ichar(kvalue(i)).GT.57)then
			xsign8  = ysign8
			xvalue8 = xvalue8
			jstatus = 1
			jterminate = i
			return
		else
c accumulate this digit
			kwork = ichar(kvalue(i)) -48
			xvalue8 = xvalue8*10.D0 + setparm_INT2DOUBLE(kwork)
		endif
	enddo

c if here, ran length of istart:iend range and have valid number
	xsign8  = ysign8
	xvalue8 = xvalue8
	jstatus = 4
	jterminate=kend+1

	return
	end


c------------------------------------------------------------------------------
c setparm_valueInteger: retrieve one valid integer and determine if keep or not.
c Incoming:  
c     kvalue(MAXCHAR_CVALUE)	input text string
c Outgoing:
c     jvalue 			(signed) integer value
c     jstatus 			= 1 : valid number
c 	  			= 0 : invalid number - nothing will be returned
c
	subroutine setparm_valueInteger(kvalue,MAXCHAR_CVALUE,jvalue,jstatus)
	
	implicit none

	integer	MAXCHAR_CVALUE,jvalue,jstatus
	character*1 kvalue(MAXCHAR_CVALUE)
	
	integer length_token,istart,iend,lvalue,lstatus

	jstatus=0
	jvalue =0

c get string length
	call setparm_strlen(kvalue,MAXCHAR_CVALUE,length_token)

	if(length_token .EQ. 0)then
		jstatus = 0
		jvalue  = 0
		return
	endif

	istart = 1
	iend   = length_token
	call setparm_get1integer(kvalue,MAXCHAR_CVALUE,istart,iend,
     +							lvalue,lstatus)

c keep lvalue if lstatus >= 1
	if(lstatus .GE. 1)then 
		jvalue = lvalue
		jstatus= 1
	endif

	return
	end

c--------------------
c setparm_get1integer:  get one valid integer from a specified text range
c Incoming:
c     kvalue(MAXCHAR_CVALUE)	input text string
c     istart,iend		retrieve from within this text range
c Outgoing:
c     lvalue 			(signed) integer value
c     lstatus			= 0 don't use
c				= 1 keep but terminate successive searches
c				= 2 keep but found delimiter
c				= 3 keep - used full text range
c
	subroutine setparm_get1integer(kvalue,MAXCHAR_CVALUE,istart,iend,
     +							lvalue,lstatus)
	
	implicit none

	integer	MAXCHAR_CVALUE,istart,iend,lvalue,lstatus
	character*1 kvalue(MAXCHAR_CVALUE)
	
	integer kstart,kend,kstatus,kterminate,mvalue,msign

	lstatus=0
	lvalue =0

	kstart = istart
	kend   = iend
	call setparm_ConstructNumber(kvalue,MAXCHAR_CVALUE,kstart,kend,
     +					msign,mvalue,kstatus,kterminate)

c keep mvalue if kstatus >= 1
	if(kstatus .GE. 1)then 
		lvalue = msign * mvalue
		lstatus= 1
	endif

	return
	end

c------------------------------------------------------------------------------
c setparm_valueVinteger:  retrieve a series of valid integers 
c Incoming:  
c     kvalue(MAXCHAR_CVALUE)	input text string
c     NumDelimits		number of delimited integers to process
c Outgoing:
c     jvalue 			(signed) integer values (MALLOC'ED in getpar().
c     jstatus 			= 1 : valid number
c 	  			= 0 : invalid number - nothing will be returned
c
	subroutine setparm_valueVinteger(kvalue,MAXCHAR_CVALUE,jvalue,
     +							NumDelimits,jstatus)
	
	implicit none

	integer	MAXCHAR_CVALUE,NumDelimits,jvalue(NumDelimits),jstatus
	character*1 kvalue(MAXCHAR_CVALUE)
	
	integer length_token,istart,iend,lvalue,lstatus,ldelimit
	integer ivalue,priorDelimit

	do ivalue=1,NumDelimits
		jvalue(ivalue) =0
	enddo
	jstatus=0

c get string length
	call setparm_strlen(kvalue,MAXCHAR_CVALUE,length_token)

	if(length_token .EQ. 0)then
		jstatus = 0
		return
	endif

c loop over known delimiting positions
	priorDelimit=0
	do ivalue=1,NumDelimits
		call setparm_FindDelimiter(kvalue,MAXCHAR_CVALUE,length_token,
     +							ivalue,ldelimit)
		if(ldelimit .EQ. priorDelimit+1)then
			lvalue=0
			lstatus=2
		else
			istart = priorDelimit + 1
			iend   = ldelimit - 1
			call setparm_get1integer(kvalue,MAXCHAR_CVALUE,istart,
     +							iend,lvalue,lstatus)
		endif

c		keep lvalue in all cases.
		if(lstatus .GE. 0)then 
			jvalue(ivalue) = lvalue
			priorDelimit = ldelimit
		endif
	enddo

	jstatus = 1

	return
	end
c------------------------------------------------------------------------------
c setparm_valueFloat: retrieve one valid real # and determine if to keep or not.
c Incoming:
c   kvalue(MAXCHAR_CVALUE)	input text string
c Outgoing;
c   yvalue			(signed) real value
c   jstatus			= 1 : valid number
c				= 0 : invalid number - nothing is returned

	subroutine setparm_valueFloat(kvalue,MAXCHAR_CVALUE,yvalue,jstatus)
	
	implicit none

	integer	MAXCHAR_CVALUE,jstatus
	character*1 kvalue(MAXCHAR_CVALUE)
	real	yvalue

	integer length_token,istart,iend,lstatus
	real*8	zvalue8

	jstatus = 0
	yvalue  = 0.0000000000

c get string length
	call setparm_strlen(kvalue,MAXCHAR_CVALUE,length_token)

	if(length_token .EQ. 0)then
		jstatus = 0
		yvalue  = 0.0000000000
		return
	endif

	istart = 1
	iend   = length_token
	call setparm_get1double(kvalue,MAXCHAR_CVALUE,istart,iend,
     +							zvalue8,lstatus)

c keep zvalue if lstatus >1 1
	if(lstatus .GE. 1)then
		jstatus = 1
		yvalue = SNGL(zvalue8)
	endif

	return
	end


c------------------------------------------------------------------------------
c setparm_get1double: get one valid double number from a specified text range.
c   The returned value is Double-Precision (R*8); the conversion back to
c   single-precision (R*4) must happen withing the calling routine.
c Incoming:
c   kvalue(MAXCHAR_CVALUE)	input text string
c   istart,iend			retrieve from within this text range
c Outgoing;
c   zvalue8			(signed) double precision (R*8) value
c   lstatus			= 0 : don't use
c				= 1 : keep but termiante successive searchers
c				= 2 : keep but found delimiter
c				= 3 : keep - used full text range

	subroutine setparm_get1double(kvalue,MAXCHAR_CVALUE,istart,iend,
     +							zvalue8,lstatus)
	
	implicit none

	integer	MAXCHAR_CVALUE,istart,iend,lstatus
	character*1 kvalue(MAXCHAR_CVALUE)
	real*8	zvalue8
	
	real*8	xvalue,xinteger,xremainder,xsign,xexponent
	real*8	ysign,yinteger
	integer kstart,kend,kstatus,kterminate

	lstatus=0
	zvalue8 =0.0000000000D0

c----------
	xsign      = 0.0D0
	xinteger   = 0.0D0
	xremainder = 0.0D0
	xexponent  = 0.0D0
	ysign	   = 0.0D0
	yinteger   = 0.0D0

c Start with integral portion of decimal number
	kstart = istart
	kend   = iend
	call setparm_ConstructNumber8(kvalue,MAXCHAR_CVALUE,kstart,kend,
     +					xsign,xinteger,kstatus,kterminate)
c	xsign    = setparm_INT2DOUBLE(jsign)
c	xinteger = setparm_INT2DOUBLE(jvalue)

c interpret validity of number based on kstatus
c   completely bad argument 
	if(kstatus .EQ. 0)then
		lstatus = 0
		zvalue8 = 0.0D0
		return

c   partial number then bad argument
	elseif(kstatus .EQ. 1)then
		lstatus = 1
		zvalue8 = xsign * xinteger
		return

c   integral portion composes entire string
	elseif(kstatus .EQ. 4)then
		lstatus = 1
		zvalue8 = xsign * xinteger
		return
	endif

c---------
c if here, either decimal point or exponent encountered.
	goto(1000,2000,3000,1000)kstatus
1000	return


c kstatus .EQ. 2, found decimal point: 
2000	continue

c first check if no remainder exists
	if(kterminate .EQ. iend)then
		lstatus = 1
		zvalue8 = xsign * xinteger
		return
	endif

c If here, some text exists denoting possible remainder
	kstart = kterminate+1
	ysign=1.0D0
	yinteger=0.0D0
	call setparm_ConstructNumber8(kvalue,MAXCHAR_CVALUE,kstart,kend,
     +					ysign,yinteger,kstatus,kterminate)

	if(kstatus .EQ. 0)then
		lstatus = 1
c		zvalue8 = xsign * xinteger
		zvalue8 = 0.0D0
		return
	elseif(ysign .LT. 0.0D0)then
		lstatus = 1
		zvalue8 = xsign * xinteger
		return	
	elseif(kstatus .EQ. 1 .OR. kstatus .EQ. 2 .OR. kstatus .EQ. 4)then
		xremainder = yinteger
		xremainder = xremainder/DBLE(10.**(kterminate-1 - kstart + 1))

		xvalue = xsign * (xinteger + xremainder)

		lstatus = 1
		zvalue8 = xvalue
		return
	elseif(kstatus .EQ. 3)then
		xremainder = yinteger
		xremainder = xremainder/DBLE(10.**(kterminate-1 - kstart + 1))

		goto 3000
	endif


c kstatus .EQ. 3: exponent encountered
3000	continue

c first check if no exponent exists
	if(kterminate .EQ. iend)then
		xvalue = xsign * (xinteger + xremainder)
		lstatus = 1
		zvalue8 = xvalue
		return
	endif

c If here, some text exists denoting possible exponent
	kstart = kterminate+1
	ysign = 1.0D0
	yinteger=0.0D0
	call setparm_ConstructNumber8(kvalue,MAXCHAR_CVALUE,kstart,kend,
     +					ysign,yinteger,kstatus,kterminate)

	if(kstatus .EQ. 0)then
		xvalue = xsign * (xinteger + xremainder)
		lstatus=1
		zvalue8 = xvalue
		return
	elseif(kstatus .GE. 1 .AND. kstatus .LE. 4)then
		xexponent = ysign*yinteger 
		xexponent = (10.000D0)**xexponent

		xvalue = xsign * (xinteger + xremainder) * xexponent

		lstatus = 1
		zvalue8 = xvalue
		return
	endif
	
	return
	end


c------------------------------------------------------------------------------
c setparm_valueVfloat:  retrieve a series of valid floating pt numbers 
c Incoming:  
c     kvalue(MAXCHAR_CVALUE)	input text string
c     NumDelimits		number of delimited integers to process
c Outgoing:
c     zvalue4 			(signed) real*4 values (MALLOC'ED in getpar().
c     jstatus 			= 1 : valid number
c 	  			= 0 : invalid number - nothing will be returned
c
	subroutine setparm_valueVfloat(kvalue,MAXCHAR_CVALUE,zvalue4,
     +							NumDelimits,jstatus)
	
	implicit none

	integer	MAXCHAR_CVALUE,NumDelimits,jstatus
	character*1 kvalue(MAXCHAR_CVALUE)
	real  zvalue4(NumDelimits)
	
	integer length_token,istart,iend,lstatus,ldelimit
	integer ivalue,priorDelimit
	real*8	xvalue8

	do ivalue=1,NumDelimits
		zvalue4(ivalue) =0.00000000000000
	enddo
	jstatus=0

c get string length
	call setparm_strlen(kvalue,MAXCHAR_CVALUE,length_token)

	if(length_token .EQ. 0)then
		jstatus = 0
		return
	endif

c loop over known delimiting positions
	priorDelimit=0
	do ivalue=1,NumDelimits
		call setparm_FindDelimiter(kvalue,MAXCHAR_CVALUE,length_token,
     +							ivalue,ldelimit)
		if(ldelimit .EQ. priorDelimit+1)then
			xvalue8=0.0D0
			lstatus=2
		else
			istart = priorDelimit + 1
			iend   = ldelimit - 1
			call setparm_get1double(kvalue,MAXCHAR_CVALUE,istart,
     +							iend,xvalue8,lstatus)
		endif

c		keep xvalue8 in all cases.
		if(lstatus .GE. 0)then 
			zvalue4(ivalue) = SNGL(xvalue8)
			priorDelimit = ldelimit
		endif
	enddo

	jstatus = 1

	return
	end



c====================================================================
c------------------------------------------------------------------------------
c setparm_valueDoubleFloat: retrieve one valid real # and determine if to keep or not.
c Incoming:
c   kvalue(MAXCHAR_CVALUE)	input text string
c Outgoing;
c   yvalue			(signed) real value
c   jstatus			= 1 : valid number
c				= 0 : invalid number - nothing is returned

	subroutine setparm_valueDoubleFloat(kvalue,MAXCHAR_CVALUE,yvalue,
     +								jstatus)
	
	implicit none

	integer	MAXCHAR_CVALUE,jstatus
	character*1 kvalue(MAXCHAR_CVALUE)
	real*8	yvalue

	integer length_token,istart,iend,lstatus
	real*8	zvalue8

	jstatus = 0
	yvalue  = 0.000D0

c get string length
	call setparm_strlen(kvalue,MAXCHAR_CVALUE,length_token)

	if(length_token .EQ. 0)then
		jstatus = 0
		yvalue  = 0.000D0
		return
	endif

	istart = 1
	iend   = length_token
	call setparm_get1double(kvalue,MAXCHAR_CVALUE,istart,iend,
     +							zvalue8,lstatus)

c keep zvalue if lstatus >1 1
	if(lstatus .GE. 1)then
		jstatus = 1
		yvalue = zvalue8
	endif

	return
	end



c------------------------------------------------------------------------------
c setparm_valueVdoublefloat:  retrieve a series of valid floating pt numbers 
c Incoming:  
c     kvalue(MAXCHAR_CVALUE)	input text string
c     NumDelimits		number of delimited integers to process
c Outgoing:
c     zvalue8 			(signed) real*8 values (MALLOC'ED in getpar().
c     jstatus 			= 1 : valid number
c 	  			= 0 : invalid number - nothing will be returned
c
	subroutine setparm_valueVdoublefloat(kvalue,MAXCHAR_CVALUE,zvalue8,
     +							NumDelimits,jstatus)
	
	implicit none

	integer	MAXCHAR_CVALUE,NumDelimits,jstatus
	character*1 kvalue(MAXCHAR_CVALUE)
	real*8  zvalue8(NumDelimits)
	
	integer length_token,istart,iend,lstatus,ldelimit
	integer ivalue,priorDelimit
	real*8	xvalue8

	do ivalue=1,NumDelimits
		zvalue8(ivalue) =0.0D0
	enddo
	jstatus=0

c get string length
	call setparm_strlen(kvalue,MAXCHAR_CVALUE,length_token)

	if(length_token .EQ. 0)then
		jstatus = 0
		return
	endif

c loop over known delimiting positions
	priorDelimit=0
	do ivalue=1,NumDelimits
		call setparm_FindDelimiter(kvalue,MAXCHAR_CVALUE,length_token,
     +							ivalue,ldelimit)
		if(ldelimit .EQ. priorDelimit+1)then
			xvalue8=0.0D0
			lstatus=2
		else
			istart = priorDelimit + 1
			iend   = ldelimit - 1
			call setparm_get1double(kvalue,MAXCHAR_CVALUE,istart,
     +							iend,xvalue8,lstatus)
		endif

c		keep xvalue8 in all cases.
		if(lstatus .GE. 0)then 
			zvalue8(ivalue) = xvalue8
			priorDelimit = ldelimit
		endif
	enddo

	jstatus = 1

	return
	end
c------------------------------------------------------------------------------
c for this token interpret value string as a character string
c jstatus = 1 : valid string
c 	  = 0 : invalid string - nothing will be returned
c
c NOTE:  There is no way to determine if this incoming string will
c 	 overflow the CHAR*xx declaration for the character variable
c	 as defined in the getpar("","s",char_variable) call.
c	 All we can do is pass back the resulting string as provided
c	 and hope that the recipient char_variable is long enough.
c
	subroutine setparm_valueString(kvalue,MAXCHAR_CVALUE,jstatus)
	
	implicit none

	integer	MAXCHAR_CVALUE,jstatus
	character*1 kvalue(MAXCHAR_CVALUE)
	
	integer length_token

	jstatus=0

c get string length
	call setparm_strlen(kvalue,MAXCHAR_CVALUE,length_token)

	if(length_token .EQ. 0)then
		jstatus = 0
		return
	endif

	jstatus = 1

	return
	end

c process a par filename if par=name is provided
c 	par=name activates this routine.
c	setparm_num, cname(),cvalue() is updated via common block.
c
c here, the acceptable par_file size is 1 MB characters.  This is a
c	hardware which can be increased, then recompiled.
c
c 04nov03 installation.
c 06nov03 modify fortran OPEN/READ of text lines.
c 28nov03 udpate for getparm() package.

	subroutine loadparm_parfile()
	implicit none

	common /globals/setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE
	integer*4	setparm_num,MAXCHAR_CTOKEN,MAXCHAR_CNAME,MAXCHAR_CVALUE

	common /cnamevalue/ MAXARGS,M_CNAME,M_CVALUE,cname,cvalue
	integer		MAXARGS,M_CNAME,M_CVALUE
	character*64	cname(1000)
	character*180	cvalue(1000)

c NOTE maximum par_file size defined here at 1 MB.
	character*64	cpar
	character*180	cfilename
	character*1	cfile(1000000)
	integer		PARFILESIZE,nchar_file, jstatus


	PARFILESIZE = 1000000
	cpar = "parmfile"

	if(cname(setparm_num) .NE. cpar) RETURN

c if here, found a par=name name/value entry.
	cfilename  = cvalue(setparm_num)
	nchar_file = 0
	jstatus    = 0
	call loadparm_GetParfile(cfilename,cfile,PARFILESIZE,nchar_file,jstatus)
	if(jstatus .LT. 0)return

	call loadparm_Preprocess(cfile,nchar_file,jstatus)

	call loadparm_Tokenize(cfile,nchar_file)

	return
	end




c jstatus = 0   parfile is OK
c         = -1  not able to open.
c         = -2  contents is too small (too few characters).
	subroutine loadparm_GetParfile(cparfile,cfile,PARFILESIZE,
     +							nchar_file,jstatus)
	implicit none
	character*180	cparfile,   c180
	integer		PARFILESIZE,nchar_file, jstatus, n180,nlength, i
	character*1	cfile(PARFILESIZE),c1(180), creturnJ
	equivalence	(c180,c1)

	n180	   = 180
	nchar_file = 0
	jstatus    = 0

	creturnJ = char(10)

	open(21,file=cparfile,err=900)
c       open(21,file=cinx,access='direct',recl=4*ny,form='unformatted')


100	continue
c	read(21,err=300)c1
	read(21,200,end=300,err=300)c180
200	format(a180)
		call setparm_strlen(c1,n180,nlength)

		do i=1,nlength
			nchar_file = nchar_file + 1
			cfile(nchar_file) = c1(i)
			if(nchar_file .EQ. PARFILESIZE)goto300
		enddo
		
c  terminate current line with a return. We'll remove this later,
c  but we need it to know where comment lines terminate.
		nchar_file = nchar_file + 1
		cfile(nchar_file) = creturnJ
		if(nchar_file .EQ. PARFILESIZE)goto300
		
		goto 100

c when here, read full contents of par file; now determine if OK.
c minimum size must be "x=y" which is 3 characters.
300	continue
	close(21)
	if(nchar_file .LT. 3)goto910

	return

c could not open file name
900	jstatus = -1
	return

c par file is too small ( .LT. 3 characters)
910	jstatus = -2
	return

	end


c this routine converts all characters between tokens to be spaces.
c  tab				     -> space
c  return (^J or ^M)		     -> space
c  back-slash & return combinations  -> space-space
c  return & # indicates comment line -> # through next return set to spaces.

	subroutine loadparm_Preprocess(cfile,nchar_file,jstatus)
	implicit none
	integer		nchar_file, jstatus, i,j,k
	character*1	cfile(nchar_file)

	character*1	cnull,ctab,cspace,cback,creturnJ,creturnM,ccomment

	cnull	= char(0)
	ctab 	= char(9)
	cback	= char(92)
	creturnJ = char(10)
	creturnM = char(13)
	ccomment = char(35)
	cspace	 = char(32)

	jstatus = 0

c (1) remove tabs
	do i = 1,nchar_file
		if(cfile(i) .EQ. ctab) cfile(i) = cspace
	enddo

c (2) remove nulls
	do i = 1,nchar_file
		if(cfile(i) .EQ. cnull) cfile(i) = cspace
	enddo

c (3) replace any ^M returns with ^J
	do i = 1,nchar_file
		if(cfile(i) .EQ. creturnM)cfile(i) = creturnJ
	enddo

c (4) remove comments: "# .... return"  denotes comment.

c Search for "#" symbol.  When found, comment symbol through trailing
c 	return become spaces.
c Then outer loop will probably re-examine some characters which were
c	part of the comment now converted to spaces.  This is a benign
c	search because they will be simply spaces.

	do i=1,nchar_file
		if(cfile(i) .EQ. ccomment)then
		    do j=i+1,nchar_file
		    	if(cfile(j) .EQ. creturnJ .OR. j .EQ. nchar_file)then
		    	    do k=i,j
		    	    	cfile(k) = cspace
		    	    enddo
			    goto300
		    	endif
		    enddo
		endif
300		continue
	enddo


c (5) replace backslash-return with two spaces
c	NOTE:  to be correct, this should compress by 2 characters, not be
c	replaced with two spaces.
	do i=2,nchar_file
		if(cfile(i-1) .EQ. cback .AND. cfile(i) .EQ. creturnJ)then
			cfile(i-1) = cspace
			cfile(i)   = cspace
		endif
	enddo

c (6) replace all returns with cspace
	do i=1,nchar_file
		if(cfile(i) .EQ. creturnJ)cfile(i) = cspace
	enddo

c at this point, the par file contents are either valid tokens or 
c sets of delimiting spaces.  No returns, tabs, or other delimiting characters.

	RETURN
	end





	subroutine loadparm_Tokenize(cfile,nchar_file)
	implicit none
	integer nchar_file
	character*1 cfile(nchar_file)

	character*1 ctoken(245),	cspace
	character*245 ctoken245
	equivalence (ctoken,ctoken245)
	integer	 iptr,	i

	cspace = char(32)

	iptr=0
100	continue
	iptr = iptr + 1
	if(iptr .GT. nchar_file)goto900
	if(cfile(iptr) .EQ. cspace) goto 100

c if here, found start of a token. clear ctoken array, then transfer.
	do i=1,245
		ctoken(i) = cspace
	enddo
	i = 0

c first save non-space, then check next char for termination space 
200	continue
	if(cfile(iptr) .EQ. cspace)goto 300
	i = i+1
	ctoken(i) = cfile(iptr)
	iptr = iptr + 1
	if(iptr .GT. nchar_file .AND. i .GT. 0)goto300
	if(iptr .GT. nchar_file .AND. i .EQ. 0)goto900
	
	goto200

c when here, we have extracted a valid token.  Now parse it and fill into
c cname() and cvalue(), incrementing setparm_num via common blocks.
300	continue
	call loadparm_parse(ctoken245)

c completed processing of this token, go look for next one
	goto100

900	continue
	return
	end

c----------------------------------------------------------------------------
	subroutine loadparm_setdir(cdirection)
	implicit none
	character*8 cdirection

	character*8 c8
	character*1 c1(8)
	equivalence (c8,c1)
	integer	    n8

	common /loadparm_dir/loadparm_direction
	integer*4	loadparm_direction

	c8 = cdirection

	n8 = 8
	call setparm_lowercase(c1,n8)

	loadparm_direction = 1

	if    (c8 .EQ. 'leading')then
		loadparm_direction = 1
	elseif(c8 .EQ. 'trailing')then
		loadparm_direction = -1
	else
		loadparm_direction = 1
	endif
	
	return
	end


c----------------------------------------------------------------------------
c=============================================================================
c=============================================================================
c=============================================================================
c UTILITIES TO INTERPRET GETPAR CHARACTER RESPONSES
c	These routines allow code to ask for options using
c	YES/NO, TRUE/FALSE, or CHOICE_A/CHOICE_B instead
c	of =0, =1 flags
c examples:
c	main option={YES/NO} switch={ON/OFF} logic={TRUE/FALSE} choice={A/B}
c
c	main option=yes switch=off logic=false choice=ASTRING
c	main option=N switch=ON logic=T choice=bstring
c
c jflag_return values:
c	 0 = no/off/false/A
c	 1 = yes/on/true/B
c	-1 = copt value not recognized as valid choice
c
c	when jflag_return = -1, keep whatever default may have been set.
c
c subroutine functions:
c	subroutine setparm_yesno4(copt,jflag_return)
c	subroutine setparm_onoff4(copt,jflag_return)
c	subroutine setparm_truefalse8(copt,jflag_return)
c	subroutine setparm_choiceAorB(copt,choiceA,choiceB,nchar,jflag_return)
c	subroutine setparm_lowercase(copt,nchar)
c	subroutine setparm_fixstring(cstring, ndeclaration)
c
c---------------------
c 08Nov98 DAO	Initial implementation
c 30jul00 DAO   getparm_lowercase: fix if() comparison of character to ASCII 
c			decimal value by using ichar().
c---------------------


c---------------------
c Yes/No or True/False
c   No = FALSE = OFF = 0
c   Yes = TRUE = ON  = 1
c default must be set outside prior to call; no default assumed here.
	subroutine setparm_yesno4(copt,jflag_return)
	character*4 copt,cwork
	character*1 c1(4)
	equivalence (cwork,c1)

	cwork = copt
	call setparm_fixstring(cwork,4)
	call setparm_lowercase(cwork,4)
	
	jset=0
	if(cwork .EQ. "yes " .OR. cwork .EQ. "y   ")then
		jflag_return=1
		jset=1
	endif
c	if(c1(1) .EQ. "y")then
c		jflag_return=1
c		jset=1
c	endif

	if(cwork .EQ. "no  " .OR. cwork .EQ. "n   ")then
		jflag_return=0
		jset=1
	endif
c	if(c1(1) .EQ. "n")then
c		jflag_return=0
c		jset=1
c	endif

	if(jset .EQ. 0)then
		jflag_return = -1
	endif
	
	RETURN
	end


c---------------------
c Yes/No or True/False
c   No = FALSE = OFF = 0
c   Yes = TRUE = ON  = 1
c default must be set outside prior to call; no default assumed here.
	subroutine setparm_onoff4(copt,jflag_return)
	character*4 copt,cwork
	character*1 c1(4)
	equivalence (cwork,c1)

	cwork = copt
	call setparm_fixstring(cwork,4)
	call setparm_lowercase(cwork,4)
	
	jset=0
	if(cwork .EQ. "on  ")then
		jflag_return=1
		jset=1
	endif

	if(cwork .EQ. "off ")then
		jflag_return=0
		jset=1
	endif

	if(jset .EQ. 0)then
		jflag_return = -1
	endif
	
	RETURN
	end

c---------------------
c Yes/No or True/False
c   No = FALSE = 0
c   Yes = TRUE = 1
c default must be set outside prior to call; no default assumed here.
	subroutine setparm_truefalse8(copt,jflag_return)
	character*8 copt,cwork
	character*1 c1(8)
	equivalence (cwork,c1)

	cwork = copt
	call setparm_fixstring(cwork,8)
	call setparm_lowercase(cwork,8)
	
	jset=0
	if(cwork .EQ. "true    " .OR. cwork .EQ. "t       ")then
		jflag_return=1
		jset=1
	endif
c	if(c1(1) .EQ. "t")then
c		jflag_return=1
c		jset=1
c	endif

	if(cwork .EQ. "false   " .OR. cwork .EQ. "f       ")then
		jflag_return=0
		jset=1
	endif
c	if(c1(1) .EQ. "f" )then
c		jflag_return=0
c		jset=1
c	endif

	if(jset .EQ. 0)then
		jflag_return = -1
	endif
	
	RETURN
	end


c---------------------
c Character choices: 
c   0 = Choice_A
c   1 = Choice_B
c default must be set outside prior to call; no default assumed here.
	subroutine setparm_choiceAorB(copt,choiceA,choiceB,nchar,jflag_return)

	character*1 copt(nchar),choiceA(nchar),choiceB(nchar)

	call setparm_fixstring(copt,nchar)
	call setparm_fixstring(choiceA,nchar)
	call setparm_fixstring(choiceB,nchar)
	
	call setparm_lowercase(copt,nchar)
	call setparm_lowercase(choiceA,nchar)
	call setparm_lowercase(choiceB,nchar)

	jmatch= -1
	do i=1,nchar
		if(copt(i) .NE. choiceA(i))goto100
	enddo
	jmatch=0
	goto500

100	continue
	do i=1,nchar
		if(copt(i) .NE. choiceB(i))goto500
	enddo
	jmatch=1
	goto500

500	jflag_return = jmatch
	RETURN
	end



c---------------------
c convert any upper case to lower case
c Upper case ASCII:   A to Z is 65 to 90
c Lower case ASCII:   a to z is 97 to 122
c shift is value of +32
	subroutine setparm_lowercase(copt,nchar)
	integer nchar
	character*1 copt(nchar)
	
	do i=1,nchar
		if(ichar(copt(i)) .GE. 65 .AND. ichar(copt(i)) .LE. 90)then
			copt(i) = char(ichar(copt(i)) +32)
		endif
	enddo
	RETURN
	end

c---------------------
c subroutine fixstring()
c
c for character strings brought in by getpar(), the strings are
c	terminated by a NULL if shorter than the declaration length.
c	Contents after the NULL are junk leftover from a prior use.
c
c	In FORTRAN, when a character string is shorter than its declaration
c	length, the string is padded to its full length with spaces (ASCII=32).
c
c	In order to compare or equate the getpar() string with a FORTRAN
c	character variable, the getpar() string must have its NULL and
c	trailing junk reset to spaces (ASCII=32).


	subroutine setparm_fixstring(cstring, ndeclaration)
	integer	    ndeclaration
	character*1 cstring(ndeclaration)

	integer	    i,j,NULL

	NULL = 0
	NSPACE=32

	do i=1,ndeclaration
		if(ichar(cstring(i)) .EQ. NULL)then
			do j= i,ndeclaration
				cstring(j) = char(NSPACE)
			enddo
			goto100
		endif
	enddo

100	continue
	RETURN
	end
c=============================================================================
c=============================================================================



	real*8 function setparm_INT2DOUBLE(intvalue)
	integer	intvalue
	real	xvalue
	real*8	dvalue

	xvalue = FLOAT(intvalue)
	dvalue = DBLE(xvalue)

	setparm_int2double = dvalue

	return

	end
