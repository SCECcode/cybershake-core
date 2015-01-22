c=============================================================
c NOTE: BEFORE COMPILING, SEE LINES 88-96 FOR SYSTEM DEPENDENCIES  02Sep05
c=============================================================
c This version of rspectra.f is modified to work within the SCEC/ITR
c	Community Modeling Environment (CME) framework.
c The original code operates on a single seismogram.
c This version operates on #X number of seismogams where #X can be
c	expressed as an N x M block of seismograms (such as 
c	NX by NY surface locations).
c
c Seismogram format is raw binary ("Surface Seismogram" or "KBO" format).
c Incoming seismograms are assumed to be velocity (but can be acceleration).
c Pre-process options for demean, vel->acceleration, and filter are installed.
c Dreger code originally assumes acceleration in CGS units.  So installed
c 	is a preprocess option to convert velocity seismograms from
c	MKS to CGS in preparation for vel->acceleration conversion.
c Input communication is via command line attributes in place
c	of the text parameter files.
c
c Blocks of original code not needed for this modification are
c	commented out using the string CME.
c
c Original output values are:
C   > 70     write(9,5) n,1/t(n),rd(n),rv(n),prv(n),aa(n),pra(n),b(n),t(n)
c real component of displacement  spectrum:    > rd(it)=z(1)
c real component of velocity      spectrum:    > rv(it)=z(2)
c real component of acceleration  spectrum:    > aa(it)=z(3)/981.0
c imaginary component of velocity spectrum:    > prv(it)=w*z(1)
c imaginary component of velocity spectrum:    > pra(it)=w*prv(it)/981.0
c relative real comp. acceleration:            > b(it)=aa(it)/amax
c
c imaginary vel = i*w* real displacement   where w = period = 2pi/T
c imaginary acc = i*w* real velocity    
c
c and   accel = xx / 981.  converts cm/s^2 to units of g.
c
c--------------------
c 17 may 05	D. Okaya	initial modification of rspectra.f
c 28 may 05	D. Okaya	replace C I/O with f90 I/O.
c 01 jun 05	D. Okaya	Butterworth filter (Stanford code).
c 02 jun 05	D. Okaya	options regarding input units. 
c 02 sep 05	D. Okaya	(1) Intel/Linux vs. other f90 portability:
c				 	file I/O record length sizing.
c 				(2) options regarding output units.
c 				(3) output decimation option.
c 04 sep 05	D. Okaya	Output options: text, all periods.
c 09 sep 05	D. Okaya	add a byteswap option right after read.
c 31 aug 07 G. Juve     Overwrite output; non-zero exit codes
c 13 jan 11 S. Callaghan Changed to subroutine called externally.
c 09 oct 12 S. Callaghan Added header to output
c=============================================================
c=============================================================

c      program SPECTRAD
c Otput format modified to WESSA's format
c
c
      subroutine spectrad(s_hdr, nx,ny,npts,dt,cin_units,cout_units,
     + 		cout_choice, cout_numperiods,filter_highHZ,
     +		capply_byteswap,ifile,ofile,seis,output_opt,psa_data)
	implicit none
	integer		maxperiod,maxdamp,maxnpt
      parameter(maxperiod=305, maxdamp=10, maxnpt=60000)
      real t(maxperiod),rd(maxperiod),rv(maxperiod),prv(maxperiod),
     , aa(maxperiod), pra(maxperiod), damp(maxdamp),
     , a(maxnpt),time(maxnpt),b(maxperiod),z(3), peat(112), scec(44)
	real outarray(maxperiod)
	real psa_data(maxperiod)
c
CME   character title*75,head*75,fdesc*20,ifile*40,ofile*40,rfile*40,
CME  , word*75,fmt*40
	character title*75,head*75,fdesc*20,ifile*256,ofile*256,rfile*120,
     ,	   word*75,fmt*40

	integer		ntpea,nscec,kg
	integer		i,k,kug,id,it,nt,n,	nt2
	real		pga,amax,v1,d1,vm,dm,vt,td,tmv,tmd,tma,w,dur1
c
      logical stp,free

	integer		strlen
	character*120	cperiod,cout_choice,cin_units,cout_units
	character*120 	cin_units_param,cout_units_param
	character*120	cout_choice_param,cout_numperiods_param
	character*120	capply_byteswap_param
	character*256	ifile_param,ofile_param
	character*120	capply_demean,capply_filter,capply_vel2accel
	character*120	capply_byteswap,cout_numperiods
	integer		nx,ny,npts,npts2,ndamp,output_opt
	integer		ierr,ichk, krecin, krecout
	integer		ix,iy,		jperiod_out,jperiod0,jperiod1
	integer		jout_format,jout_numperiods
	real		dt,fac,filter_lowHZ,filter_highHZ
	real		work(maxnpt), 	period_out,	scalar
	real		seis(*)
	type header
	sequence
	character*8 version
	character*8 site_name
	character*8 padding
	integer source_id
	integer rupture_id
	integer rup_var_id
	real dt
	integer nt
	integer comps
	real det_max_freq
	real stoch_max_freq
	end type

	type(header) s_hdr
	real, dimension (:,:),	allocatable :: outvalue

	common	/knodes/koutnx,koutx0,koutx1,koutdx, koutx,
     +			koutny,kouty0,kouty1,koutdy, kouty
	integer		koutnx,koutx0,koutx1,koutdx, koutx,
     +			koutny,kouty0,kouty1,koutdy, kouty

c--------------------------------------------------------------
c OS/FORTRAN PORTABILITY for binary file I/O (Linux/Intel words; others bytes)
	integer		INTELf90,OTHERf90
	parameter	(INTELf90=1, OTHERf90=4)

	common		/F90_REC_COUNTING/REC_COUNT
	integer				  REC_COUNT

c set here:
c	REC_COUNT = INTELf90
	REC_COUNT = OTHERf90
c--------------------------------------------------------------

c
1     format(/,
     , 10x,' *************************************************',/,
     , 10x,' *          Program SPECTRA version 1.0          *',/,
     , 10x,' * Copyright GEOMATRIX Consultants, January 1986 *',/,
     , 10x,' *            written by Robert Youngs           *',/,
     , 10x,' *************************************************',/)
        data  ntpea/112/, peat/
     +.2000000E+02,.1499925E+02,.1399972E+02,.1300052E+02,.1200048E+02,
     +.1099989E+02,.1000000E+02,.9500286E+01,.9000090E+01,.8499788E+01,
     +.8000000E+01,.7500187E+01,.6999860E+01,.6499837E+01,.5999880E+01,
     +.5499945E+01,.5000000E+01,.4800076E+01,.4600028E+01,.4400053E+01,
     +.4199916E+01,.4000000E+01,.3799970E+01,.3599971E+01,.3399973E+01,
     +.3200000E+01,.3000030E+01,.2800022E+01,.2599969E+01,.2399981E+01,
     +.2199978E+01,.2000000E+01,.1899985E+01,.1800018E+01,.1699986E+01,
     +.1600000E+01,.1499992E+01,.1399991E+01,.1300001E+01,.1200005E+01,
     +.1100001E+01,.1000000E+01,.9500016E+00,.9000010E+00,.8500004E+00,
     +.8000000E+00,.7500019E+00,.7000007E+00,.6666667E+00,.6500006E+00,
     +.5999989E+00,.5500006E+00,.5000000E+00,.4800008E+00,.4600007E+00,
     +.4399995E+00,.4200005E+00,.4000000E+00,.3799998E+00,.3600010E+00,
     +.3399996E+00,.3200000E+00,.3000003E+00,.2898550E+00,.2799999E+00,
     +.2600003E+00,.2399998E+00,.2200002E+00,.2000000E+00,.1899999E+00,
     +.1799998E+00,.1700001E+00,.1600000E+00,.1499999E+00,.1399999E+00,
     +.1333333E+00,.1300000E+00,.1200001E+00,.1100000E+00,.1000000E+00,
     +.9500015E-01,.9000009E-01,.8500005E-01,.8000000E-01,.7500019E-01,
     +.7000007E-01,.6666667E-01,.6500007E-01,.5999988E-01,.5500005E-01,
     +.5000000E-01,.4800008E-01,.4600007E-01,.4399995E-01,.4200005E-01,
     +.4000000E-01,.3571429E-01,.3225806E-01,.2941176E-01,.2500000E-01,
     +.2222222E-01,.2000000E-01,.1818182E-01,.1666667E-01,.1538462E-01,
     +.1428571E-01,.1333333E-01,.1250000E-01,.1176471E-01,.1111111E-01,
     +.1052632E-01,.1000000E-01/
        data  nscec/44/, scec/
     +.1000000E+02,.9500286E+01,.9000090E+01,.8499788E+01,
     +.8000000E+01,.7500187E+01,.6999860E+01,.6499837E+01,.5999880E+01,
     +.5499945E+01,.5000000E+01,.4800076E+01,.4600028E+01,.4400053E+01,
     +.4199916E+01,.4000000E+01,.3799970E+01,.3599971E+01,.3399973E+01,
     +.3200000E+01,.3000030E+01,.2800022E+01,.2599969E+01,.2399981E+01,
     +.2199978E+01,.2000000E+01,.1666667E+01,.1428571E+01,.1250000E+01,
     +.1111111E+01,.1000000E+01,.6666667E+00,.5000000E+00,.4000000E+00,
     +.3333333E+00,.2857143E+00,.2500000E+00,.2222222E+00,.2000000E+00,
     +.1666667E+00,.1428571E+00,.1250000E+00,.1111111E+00,.1000000E+00/
c
CMEc  print 1
CME   fdesc='input '
CME   close(7)
CME   call flchk('i',fdesc,ifile,stp)
CME   if(stp) STOP
CME   open(7,file=ifile)
CME99  read (7,2,end=2999) title
CME   read (7,*) nt,(t(i),i=1,nt)
CMEc
CMEc  if nt = 0 load TI periods
CMEc     nt =-1 load SCEC periods
CMEc
CME      if(nt.eq.0) then
CME        nt=ntpea
CME        do 10 i=1,nt
CME10     t(i)=peat(i)
CME   elseif(nt.eq.-1) then
CME     nt=nscec
CME     do 11 i=1,nt
CME11       t(i)=scec(i)
CME      ENDif
CMEc
CMEc...read input data
CMEc
CME   read(7,3) rfile
CME   read(7,3) ofile
CME   length=index(ofile,'  ')
CME   read(7,3) fmt
CME   if(fmt.eq.'*'.or.fmt.eq.'standard'.or.fmt.eq.'STANDARD') then
CME     free=.true.
CME    else
CME     free=.false.
CME   ENDif
CME   read(7,*) nhead,npts,fac,dt,ndamp,(damp(i),i=1,ndamp)
CMEc   added close doug
CME   close (7)
CME   if(fmt.eq.'standard'.or.fmt.eq.'STANDARD') then
CME     nhead = 2
CME   endif

CME   open(8,file=rfile)
CME   if(nhead.eq.0) then
CME     head=title
CME   else
CME     read(8,2) head
CME     if(nhead.gt.1) then
CME       do 20 i = 2,nhead
CME20     read(8,2) word
CME     ENDif
CME     if(fmt.eq.'standard'.or.fmt.eq.'STANDARD') then
CME       read(8,*) npts, dt
CME     endif
CME   ENDif
CME   ims=0
CME   if(free) then
CME     if(dt.eq.0) then
CME       read(8,*) (time(l),a(l),l=1,npts)
CME      else
CME       read(8,*) (a(l),l=1,npts)
CME     ENDif
CME    else
CME     if(dt.eq.0) then
CME       read(8,fmt) (time(l),a(l),l=1,npts)
CME      else
CME       read(8,fmt) (a(l),l=1,npts)
CME     ENDif
CME   ENDif
CME   close (8)

c	copy parameter arguments (from C) into fortran strings

       if(cout_numperiods .EQ. 'all' .OR. cout_numperiods .EQ. 'ALL')then
               jout_numperiods = 0
       else
c               ichk = setparmf("surfseis_rspectra_period",period_out)
               jout_numperiods = 1
       endif

c     print locations
c      write(*,*) loc(nx),loc(ny),loc(npts),loc(dt),
c     +          loc(cin_units) ,loc(cout_units),
c     +          loc(cout_choice), loc(cout_numperiods),loc(filter_highHZ),
c     +          loc(capply_byteswap),loc(ifile),loc(ofile),loc(seis)

	call getparms(ifile,ofile,nx,ny,npts,dt, fac,ndamp,damp,maxdamp,
     +		cperiod, capply_demean,capply_filter,capply_vel2accel,
     +		filter_lowHZ,filter_highHZ,period_out,capply_byteswap,
     +		cout_choice,cin_units,cout_units,jout_format,jout_numperiods)

	call runtime_doc(nx,ny,npts,dt, fac,ndamp,damp,maxdamp,
     +		cperiod, capply_demean,capply_filter,capply_vel2accel,
     +		filter_lowHZ,filter_highHZ,period_out,capply_byteswap,
     +		cout_choice,cin_units,cout_units,jout_format,jout_numperiods)

	call check_units_vel2accel(cin_units,capply_vel2accel)


c     if nt = 0 load TI periods
c        nt =-1 load SCEC periods
c
CME   if(nt.eq.0) then
CME        nt=ntpea
CME        do 10 i=1,nt
CME10      t(i)=peat(i)
CME   elseif(nt.eq.-1) then
CME        nt=nscec
CME        do 11 i=1,nt
CME11      t(i)=scec(i)
CME   ENDif
      if(cperiod .EQ. 'TI')then
           nt=ntpea
           do 10 i=1,nt
10         t(i)=peat(i)
      elseif(cperiod .EQ. 'SCEC') then
           nt=nscec
           do 11 i=1,nt
11         t(i)=scec(i)
      ENDif

	call get_periodout(t,nt,jout_numperiods,period_out,jperiod_out,
     +							jperiod0,jperiod1)

c define array of sample times
	do i = 1,npts
		time(i) = float(i-1) * dt
	enddo

	ierr = 0
	call open_files(ofile,npts,koutnx,jout_format,
     +				jout_numperiods,nt,ierr,output_opt)
	if(ierr .EQ. 1) write(6,21)
	if(ierr .EQ. 2) write(6,22)
	if(ierr .GT. 0) call exit(1)
21	format('Surfseis_rspectra> ERROR OPENING INPUT DATA FILE')
22	format('Surfseis_rspectra> ERROR CREATING OUTPUT DATA FILE')

c get work memory
c	allocate (outvalue(nx,ny))
c	do 1999 iy = 1,ny
c	do 1999 ix = 1,nx
c1999	outvalue(ix,iy) = 0.
 	allocate (outvalue(koutnx,koutny))
 	do 1999 iy = 1,koutny
 	do 1999 ix = 1,koutnx
 1999	outvalue(ix,iy) = 0.

c At beginning of a SurfSeis file - loop over NY x NX seismograms
c We are reading velocity seismograms, not acceleration.
	koutx = 0
	kouty = 0
	krecin = 0
	krecout = 0

c If all periods and binary, write header before loop
	if(jout_numperiods .EQ. 0)then
                if(jout_format .EQ. 0) then
			if(output_opt.EQ.1)then
				write(22) s_hdr%version
				write(22) s_hdr%site_name
				write(22) s_hdr%padding
				write(22) s_hdr%source_id, s_hdr%rupture_id
				write(22) s_hdr%rup_var_id, s_hdr%dt
				write(22) s_hdr%nt, s_hdr%comps
				write(22) s_hdr%det_max_freq, s_hdr%stoch_max_freq
			elseif (output_opt.EQ.0) then
                                write(22,rec=1) s_hdr%version
                                write(22,rec=2) s_hdr%site_name
                                write(22,rec=3) s_hdr%padding
                                write(22,rec=4) s_hdr%source_id, s_hdr%rupture_id
                                write(22,rec=5) s_hdr%rup_var_id, s_hdr%dt
                                write(22,rec=6) s_hdr%nt, s_hdr%comps
                                write(22,rec=7) s_hdr%det_max_freq, s_hdr%stoch_max_freq
				krecout=7
			endif
		endif
	endif
			


c	do 3000 iy = 1,ny
c	do 2000 ix = 1,nx
	do 3000 iy = kouty0,kouty1,koutdy
		kouty = kouty + 1
		koutx = 0
	do 2000 ix = koutx0,koutx1,koutdx
c	krecin = krecin + 1
c	write(*,*) iy, nx, ix
	krecin = (iy-1)*nx + ix
	do i=1,npts
		a(i) = seis((krecin-1)*npts+i)
	enddo
c	read (21,rec=krecin) (a(i),i=1,npts)

	if(capply_byteswap .EQ. 'yes')call seismogram_byteswap(a,npts)

	if(cin_units .EQ. 'meterpersec' .OR. cin_units .EQ. 'meterpersec2')
     +				call seismogram_mks2cgs(a,npts)

	if(capply_demean .EQ. 'yes')call seismogram_demean(a,npts,work)
	if(capply_vel2accel .EQ. 'yes')call seismogram_vel2accel(a,npts,work,dt)
	if(capply_filter .EQ. 'yes')then
c	    call seismogram_filter(a,npts,work,dt,filter_lowHZ,filter_highHZ)
	    npts2 = npts+2
	    call seismogram_stanford(a,npts,npts2,work,dt,
     +					filter_lowHZ,filter_highHz)
	endif

      pga = 0.0
      kg=npts
      if(dt.ne.0) then
        do 30 i=1,kg
        pga = max(pga, abs(a(i)))
30      time(i)=real(i-1)*dt
      ENDif
      if(fac.ne.1.0) then
        do 35 i=1,kg
35        a(i)=a(i)*fac
      ENDif
c
c       write (*,*) ofile
CME     open(9,file=ofile)
CME     write (9,*)' RESPONSE SPECTRA'
CME     write (9,*) title
CME     write (9,'(a17,a11,a25,i6,a9,i6,a7,f5.2)')
CME  1  'INPUT FILE NAME: ', rfile,
CME  1  ', WINDOW UTILIZED: POINT ', 1, ' TO POINT', npts,
CME  1  ' PGA = ',pga


CME NOTE: for below block of code, the write statement is commented
CME	out by the author.  The only parameter required beyond this
CME	block is 'amax'.  As a result we comment out this whole block
CME	and just preserve the calculation of amax.
c
c...find max velocity, displacement, and acceleration
c
CME   amax=0.0
CME   v1=0.
CME   d1=0.
CME   vm=0.
CME   dm=0.
CME   do 40 k=1,kg-1
CME     vt=v1+0.5*(time(k+1)-time(k))*(a(k)+a(k+1))
CME     td=d1+(time(k+1)-time(k))*(v1+(time(k+1)-time(k))*a(k)/3.+
CME  ,   (time(k+1)-time(k))*a(k+1)/6.)
CME     if(abs(vt).gt.vm) then
CME       vm=abs(vt)
CME       tmv=time(k+1)
CME     ENDif
CME     if(abs(td).gt.dm) then
CME       dm=abs(td)
CME       tmd=time(k+1)
CME     ENDif
CME     if(abs(a(k)).gt.amax) then
CME       amax=abs(a(k))
CME       tma=time(k+1)
CME     ENDif
CME     v1=vt
CME     d1=td
CME40   CONTINUE
CME   vm=981.*vm
CME   dm=981.*dm
CMEc  write (6,4) dm,tmd,vm,tmv,amax,tma

	amax = 0.0
	do 40 k = 1,kg-1
	    if(abs(a(k)) .GT. amax)then
		amax = abs(a(k))
	    endif
40	continue
	IF(AMAX .EQ. 0.)AMAX = 1.
c
c... convert acc's to cm/sec**2
c
      do 50 i=1,kg
50      a(i)=a(i)*981.0

c
c   compute response
c
      kug=kg-1
      do 100 id=1,ndamp
CME     do 60 it=1,nt
        do 60 it=jperiod0,jperiod1
          w=4.*asin(1.0)/t(it)
c
c...compute response
c
          if(dt.eq.0.0 .or. t(it).lt.10.*dt) then
            call ucmpmx(dur1,kug,a,time,t(it),w,damp(id),z)
           else
            call cmpmax(dur1,kug,a,t(it),w,damp(id),dt,z)
          ENDif
          rd(it)=z(1)
          rv(it)=z(2)
          aa(it)=z(3)/981.0
          prv(it)=w*z(1)
          pra(it)=w*prv(it)/981.0
          b(it)=aa(it)/amax
60        CONTINUE
c
c   output values
c
c       write(9,*)' dampimg = ',damp(id)
CME     write (9,'(i5,1x,f4.3,a63)') nt, damp(id),
CME  1 'DATA- MTIT,EQN,NO PTS,DAMP,NO,FRQ,RD,RV,PRV,AA,PAA,MAG RAT,PER'
CME
CME     do 70 n=1,nt
CME70     write(9,5) n,1/t(n),rd(n),rv(n),prv(n),aa(n),pra(n),b(n),t(n)

100     CONTINUE

CME199    format(27e12.4)
CME     if(nt.eq.nscec) then
CME       write(*,199) (aa(i)*980,i=1,nt)
CME     ENDIF
CME     close(9)
CME   go to 99


c scalar for output units (at this point r,v in CGS; a in g)
	scalar = 1.
	if(cout_units .EQ. 'meters')		scalar = .01
	if(cout_units .EQ. 'cm')		scalar = 1.
	if(cout_units .EQ. 'meterpersec')	scalar = .01
	if(cout_units .EQ. 'cmpersec')		scalar = 1.
c accel in g
	if(cout_units .EQ. 'meterpersec2')	scalar = 9.80
	if(cout_units .EQ. 'cmpersec2')		scalar = 980.
	if(cout_units .EQ. 'unitsofg')		scalar = 1.
	if(cout_units .EQ. 'percentg')		scalar = 100.
	if(scalar .EQ. 0.)scalar = 1.

c select output choice and scale per user request.
	if(cout_choice .EQ. 'rd')then
		do 200 it=jperiod0,jperiod1
200		outarray(it) = scalar * rd(it)
	elseif(cout_choice .EQ. 'rv')then
		do 210 it=jperiod0,jperiod1
210		outarray(it) = scalar * rv(it)
	elseif(cout_choice .EQ. 'prv')then
		do 220 it=jperiod0,jperiod1
220		outarray(it) = scalar * prv(it)
	elseif(cout_choice .EQ. 'aa')then
		do 230 it=jperiod0,jperiod1
230		outarray(it) = scalar * aa(it)
	elseif(cout_choice .EQ. 'paa')then
		do 240 it=jperiod0,jperiod1
240		outarray(it) = scalar * pra(it)
	elseif(cout_choice .EQ. 'magrat')then
		do 250 it=jperiod0,jperiod1
250		outarray(it) = scalar * b(it)
	endif


c SAVE this seismogram's value or OUTPUT

c One period: save desired parameter/period for later output
	if(jout_numperiods .EQ. 1)then
		koutx = koutx + 1
		outvalue(koutx,kouty) = outarray(jperiod_out)

c All periods: output outarray() now (jout_format=0 is binary, =1 is text).
	elseif(jout_numperiods .EQ. 0)then
		if(jout_format .EQ. 0)then
c			write(22,rec=krecout) (outarray(i),i=1,nt)
c Use funny-sized records (8-byte) so we can get header info in there
			if(output_opt.EQ.1) then
				write(22) (outarray(i),i=1,nt)
			elseif(output_opt.EQ.0) then
				do 299 i=1, nt, 2
c					write(*,*) "Writing record ", i
					krecout = krecout + 1
					write(22, rec=krecout) outarray(i),outarray(i+1)
299				enddo
			elseif(output_opt.EQ.2) then
				do 301 i=1, nt
					psa_data(krecin*nt+i) = outarray(i)
301				enddo
			endif
		elseif(jout_format .EQ. 1)then
			if(cperiod .EQ. 'TI')then
				write(22,300) ix,iy,(outarray(i),i=1,nt)
			elseif(cperiod .EQ. 'SCEC')then
				write(22,310) ix,iy,(outarray(i),i=1,nt)
			endif
300			format(2i7,112e12.5)
310			format(2i7,27e12.5)
		endif
	endif

c end of individual seismogram work: continue loops over seismograms
2000	continue
3000	continue


c If ONE period, write results at this time
	if(jout_numperiods .EQ. 1)then
		if(jout_format .EQ. 0)then
c 		now save results - grid of NX x NY elements:
		    krecout=0
		    do iy = 1,koutny
		    	krecout = krecout + 1
		    	write(22,rec=krecout) (outvalue(i,iy),i=1,koutnx)
		    enddo

		elseif(jout_format .EQ. 1)then
		    do iy = 1,koutny
		    do ix = 1,koutnx
			write(22,400)koutx0+koutdx*(ix-1),kouty0+koutdy*(iy-1),
     +							outvalue(ix,iy)
400			format(2i7,e12.5)
		    enddo
		    enddo

		endif
	endif

c shutdown
	deallocate(outvalue)
c	close(21)
	write(*,*) "Closing output file."
	close(22)
	write(*,*) "Output file closed."

c
2     format(a75)
3     format(a40)
4     format(/,' max disp = ',f10.5,' cm at time = ', f10.5, ' secs',/
     , ' max vel = ',f10.5,' cm/sec at time = ', f10.5, ' secs',/
     , ' max accel = ',f10.5,' g at time = ', f10.5, ' secs')
5     format(i3,8e15.7)
6     format(i7,' points uneq dt  parameters are t,sd,sv,psv,sa,mr',
     , /,i5)
7     format(i7,' points',f6.3,' dt  parameters are t,sd,sv,psv,sa,mr',
     , /,i5)
c    2999  stop
2999  write(*,*) "Finished PSA calculation."
      END

      subroutine ucmpmx(dur1,kug,ug,time,pr,w,d,z)
      real ug(*),time(*),z(*),t(3),c(3),x(2,3)
c
      wd=sqrt(1.-d*d)*w
      w2=w*w
      w3=w2*w
      do 10 i=1,3
        x(1,i)=0.
 10     z(i)=0.
      f2=1./w2
      f3=d*w
      f4=1./wd
      f5=f3*f4
      f6=2.*f3
      do 100 k=1,kug
        dt=time(k+1)-time(k)
        ns=nint(10.*dt/pr)+1
        dt=dt/real(ns)
        f1=2.*d/w3/dt
        e=exp(-f3*dt)
        g1=e*sin(wd*dt)
        g2=e*cos(wd*dt)
        h1=wd*g2-f3*g1
        h2=wd*g1+f3*g2
        dug=(ug(k+1)-ug(k))/real(ns)
        g=ug(k)
        z1=f2*dug
        z3=f1*dug
        z4=z1/dt
        do 100 is=1,ns
          z2=f2*g
          b=x(1,1)+z2-z3
          a=f4*x(1,2)+f5*b+f4*z4
          x(2,1)=a*g1+b*g2+z3-z2-z1
          x(2,2)=a*h1-b*h2-z4
          x(2,3)=-f6*x(2,2)-w2*x(2,1)
          do 80 l=1,3
            c(l)=abs(x(2,l))
            if(c(l).gt.z(l)) then
              z(l)=c(l)
              t(l)=time(k)+is*dt+dur1
            ENDif
80          x(1,l)=x(2,l)
          g=g+dug
100       CONTINUE
c     write(6,1) pr,(t(l),l=1,3)
      return
1     format(' ucmpmx t=',f6.3,' td = ',f8.4,' tv = ',f8.4,' ta = ',
     , f8.4)
      end

      subroutine cmpmax(dur1,kug,ug,pr,w,d,dt,z)
      real ug(*),x(2,3),t(3),z(*),c(3)
c
      wd=sqrt(1.-d*d)*w
      w2=w*w
      w3=w2*w
      do 10 i=1,3
        x(1,i)=0.
10       z(i)=0.
      f1=2.*d/(w3*dt)
      f2=1./w2
      f3=d*w
      f4=1./wd
      f5=f3*f4
      f6=2.*f3
      e=exp(-f3*dt)
      g1=e*sin(wd*dt)
      g2=e*cos(wd*dt)
      h1=wd*g2-f3*g1
      h2=wd*g1+f3*g2
      do 100 k=1,kug
        dug=ug(k+1)-ug(k)
        z1=f2*dug
        z2=f2*ug(k)
        z3=f1*dug
        z4=z1/dt
        b=x(1,1)+z2-z3
        a=f4*x(1,2)+f5*b+f4*z4
        x(2,1)=a*g1+b*g2+z3-z2-z1
        x(2,2)=a*h1-b*h2-z4
        x(2,3)=-f6*x(2,2)-w2*x(2,1)
        do 80 l=1,3
          c(l)=abs(x(2,l))
          if(c(l).gt.z(l)) then
            z(l)=c(l)
            t(l)=dt*real(k)+dur1
          ENDif
80        x(1,l)=x(2,l)
100     CONTINUE
c     write(6,1) pr,(t(l),l=1,3)
      return
1     format(' cmpmax t=',f6.3,' td = ',f8.4,' tv = ',f8.4,' ta = ',
     , f8.4)
      end

      subroutine flchk(ftype,fdesc,filen,stp)
      character type,ftype,filen*(*),answer*40,fdesc*20
      logical fexst,stp
c
      if(ftype.eq.'d') then
        type='i'
        go to 20
       elseif(ftype.eq.'w') then
        type='o'
        go to 20
       else
        type=ftype
      ENDif
c     print'('' enter '',a20,''file name(q to quit): '')',fdesc
10    read(5,1) answer
      if(answer.eq.'q'.or.answer.eq.'Q') then
        stp=.true.
        return
       else
        filen=answer
        stp=.false.
      ENDif
20      inquire(file=filen,exist=fexst)
        if(type.eq.'i') then
          if(.not.fexst) then
            print'(1x,a20,''file '',a30,'' does not exist'',/,
     ,       '' enter new name or q to quit: '')',fdesc,filen
            go to 10
          END if
         else if(type.eq.'o') then
          if(fexst) then
            print'(1x,a20,''file '',a30,'' exists'',/,'' enter y to '',
     ,       ''overwrite, q to quit or new name:  '')',fdesc,filen
            read(5,1) answer
            if (answer.eq.'y' .or. answer.eq.'Y') then
              stp=.false.
             else if(answer.eq.'q'.or.answer.eq.'Q') then
              stp=.true.
             else
              filen=answer
              go to 20
            END if

          END if
        END if
c     END if
      return
1     format(a40)
      END



c==================================================================
	subroutine getparms(ifile,ofile,nx,ny,npts,dt, fac,ndamp,damp,maxdamp,
     +		cperiod, capply_demean,capply_filter,capply_vel2accel,
     +		filter_lowHZ,filter_highHZ,period_out,capply_byteswap,
     +		cout_choice,cin_units,cout_units,jout_format,jout_numperiods)
	implicit none

	character*256	ifile,ofile
	character*120   cperiod, cout_choice,cin_units,cout_units
	character*120	capply_demean,capply_filter,capply_vel2accel
	character*120	cout_format,cout_numperiods,	capply_byteswap
	integer		nx,ny,npts,ndamp,maxdamp
	integer		jout_format,jout_numperiods
	real		dt,fac,damp(maxdamp),filter_lowHZ,filter_highHZ
	real		period_out

	common	/knodes/koutnx,koutx0,koutx1,koutdx, koutx,
     +			koutny,kouty0,kouty1,koutdy, kouty
	integer		koutnx,koutx0,koutx1,koutdx, koutx,
     +			koutny,kouty0,kouty1,koutdy, kouty

	integer setparms,setparmi,setparmf
	integer	ichk,numfiles,ierr,	i,	kout_units

	write(*,*) nx

	fac = .0010204
	do i=1,maxdamp
		damp(i) = 0.
	enddo	
	ndamp = 1
	damp(ndamp) = .05

	cperiod = 'SCEC'
c	cin_units   = 'meterpersec'
c	cout_choice = 'aa'
c	cout_units  = 'unitsofg'

c	capply_byteswap = 'no'
	capply_demean = 'yes'
	capply_filter = 'yes'
	capply_vel2accel = 'yes'
	filter_lowHZ = 0.
c	filter_highHZ = 0.
c	period_out = 2.

	jout_format = 0
c	jout_numperiods = 1

	ichk = 0
	numfiles = 0
	kout_units = 1

c	call loadparm('trailing')

c	numfiles = numfiles + setparms("in",ifile)
c	numfiles = numfiles + setparms("out",ofile)
	numfiles = 2

c get dimensions from input mesh size
C commented out since we are passing things in instead
c	ichk = setparmi("mesh_nx",nx)
c	ichk = setparmi("mesh_ny",ny)
c	ichk = setparmi("simulation_timesamples",npts)
c	ichk = setparmf("simulation_dt",dt)

c override with output (possibly decimated) mesh
c	ichk = setparmi("simulation_out_pointsX",nx)
c	ichk = setparmi("simulation_out_pointsY",ny)
c	ichk = setparmi("simulation_out_timesamples",npts)
c	ichk = setparmf("simulation_out_timeskip",dt)

c rspectra_specific parameters:
c	ichk = setparms("surfseis_rspectra_output_type",cout_choice)
c	ichk = setparmf("surfseis_rspectra_galstog",fac)

c	ichk = setparmi("surfseis_rspectra_numdampcoeffs",ndamp)
c	ichk = setparmf("surfseis_rspectra_dampcoeff1",damp(1))
c	ichk = setparmf("surfseis_rspectra_dampcoeff2",damp(2))
c	ichk = setparmf("surfseis_rspectra_dampcoeff3",damp(3))
c	ichk = setparmf("surfseis_rspectra_dampcoeff4",damp(4))
c	ichk = setparmf("surfseis_rspectra_dampcoeff5",damp(5))
c	ichk = setparmf("surfseis_rspectra_dampcoeff6",damp(6))
c	ichk = setparmf("surfseis_rspectra_dampcoeff7",damp(7))
c	ichk = setparmf("surfseis_rspectra_dampcoeff8",damp(8))
c	ichk = setparmf("surfseis_rspectra_dampcoeff9",damp(9))
c	ichk = setparmf("surfseis_rspectra_dampcoeff10",damp(10))

c	ichk = setparms("surfseis_rspectra_periodtype",cperiod)
c	ichk = setparms("surfseis_rspectra_period",cout_numperiods)
c	if(cout_numperiods .EQ. 'all' .OR. cout_numperiods .EQ. 'ALL')then
c		jout_numperiods = 0
c	else 
c		ichk = setparmf("surfseis_rspectra_period",period_out)
c		jout_numperiods = 1
c	endif

c	ichk = setparms("surfseis_rspectra_seismogram_units",cin_units)
c	kout_units = setparms("surfseis_rspectra_output_units",cout_units)
c	ichk = setparms("surfseis_rspectra_output_format",cout_format)

c these are yes/no answers
c	ichk = setparms("surfseis_rspectra_apply_byteswap",capply_byteswap)
c	ichk = setparms("surfseis_rspectra_apply_demean",capply_demean)
c	ichk = setparms("surfseis_rspectra_apply_filter",capply_filter)
c	ichk = setparms("surfseis_rspectra_apply_vel2accel",capply_vel2accel)

		if(dt .GT. 0)filter_highHZ = 1/(2.*dt)
c	ichk = setparmf("surfseis_rspectra_apply_filter_lowHZ",filter_lowHZ)
c	ichk = setparmf("surfseis_rspectra_apply_filter_highHZ",
c     +								filter_highHZ)

		koutnx = nx
		koutx0 = 1
		koutx1 = nx
		koutdx = 1
			koutny = ny
			kouty0 = 1
			kouty1 = ny
			koutdy = 1
c	ichk = setparmi("surfseis_rspectra_out_pointsXstart",koutx0)
c	ichk = setparmi("surfseis_rspectra_out_pointsXend"  ,koutx1)
c	ichk = setparmi("surfseis_rspectra_out_pointsXdel"  ,koutdx)
c	ichk = setparmi("surfseis_rspectra_out_pointsYstart",kouty0)
c	ichk = setparmi("surfseis_rspectra_out_pointsYend"  ,kouty1)
c	ichk = setparmi("surfseis_rspectra_out_pointsYdel"  ,koutdy)

c	call endparm()

	ierr = 0
	if(numfiles .NE. 2)	ierr = 1
	if(nx*ny*npts .EQ. 0)	ierr = 1
	if(dt .EQ. 0.)		ierr = 1
	if(fac .EQ. 0.)		ierr = 1
	if(ndamp .EQ. 0)	ierr = 1

	if(cperiod .EQ. 'SCEC' .OR. cperiod .EQ. 'scec')then
			cperiod='SCEC'
	else
			cperiod='TI'
	endif

	if(capply_byteswap .EQ. 'yes' .OR. capply_byteswap .EQ. 'YES')then
			capply_byteswap = 'yes'
	else
			capply_byteswap = 'no'
	endif

	if(capply_demean .EQ. 'yes' .OR. capply_demean .EQ. 'YES')then
			capply_demean = 'yes'
	else
			capply_demean = 'no'
	endif
	if(capply_filter .EQ. 'yes' .OR. capply_filter .EQ. 'YES')then
			capply_filter = 'yes'
	else
			capply_filter = 'no'
	endif
	if(capply_vel2accel .EQ. 'yes' .OR. capply_vel2accel .EQ. 'YES')then
			capply_vel2accel = 'yes'
	else
			capply_vel2accel = 'no'
	endif

	if(cin_units .EQ. 'meterpersec' .OR. cin_units .EQ. 'METERPERSEC')then
			cin_units = 'meterpersec'
	elseif(cin_units .EQ. 'cmpersec' .OR. cin_units .EQ. 'CMPERSEC')then
			cin_units = 'cmpersec'
	elseif(cin_units .EQ. 'meterpersec2' .OR. 
     +					cin_units .EQ. 'METERPERSEC2')then
			cin_units = 'meterpersec2'
	elseif(cin_units .EQ. 'cmpersec2' .OR. cin_units .EQ. 'CMPERSEC2')then
			cin_units = 'cmpersec2'
	else
		write(*,*) "cin_units",cin_units // ":"
		ierr = 1
	endif

	if(cout_choice .EQ. 'rd' .OR. cout_choice .EQ. 'RD')then
			cout_choice = 'rd'
	elseif(cout_choice .EQ. 'rv' .OR. cout_choice .EQ. 'RV')then
			cout_choice = 'rv'
	elseif(cout_choice .EQ. 'prv' .OR. cout_choice .EQ. 'PRV')then
			cout_choice = 'prv'
	elseif(cout_choice .EQ. 'aa' .OR. cout_choice .EQ. 'AA')then
			cout_choice = 'aa'
	elseif(cout_choice .EQ. 'paa' .OR. cout_choice .EQ. 'PAA')then
			cout_choice = 'paa'
	elseif(cout_choice .EQ. 'magrat' .OR. cout_choice .EQ. 'MAGRAT')then
			cout_choice = 'magrat'
	else
		write(*,*) "cout_choice", cout_choice
		ierr = 1
	endif

	if(cout_units .EQ. 'meter' .OR. cout_units .EQ. 'METER')then
			cout_units = 'meter'
	elseif(cout_units .EQ. 'cm' .OR. cout_units .EQ. 'CM')then
			cout_units = 'cm'
	elseif(cout_units .EQ. 'meterpersec' .OR. 
     +					cout_units .EQ. 'METERPERSEC')then
			cout_units = 'meterpersec'
	elseif(cout_units .EQ. 'cmpersec' .OR. cout_units .EQ. 'CMPERSEC')then
			cout_units = 'cmpersec'
	elseif(cout_units .EQ. 'meterpersec2' .OR. 
     +					cout_units .EQ. 'METERPERSEC2')then
			cout_units = 'meterpersec2'
	elseif(cout_units .EQ. 'cmpersec2' .OR. cout_units .EQ. 'CMPERSEC2')then
			cout_units = 'cmpersec2'
	elseif(cout_units .EQ. 'unitsofg' .OR. cout_units .EQ. 'UNITSOFG')then
			cout_units = 'unitsofg'
	elseif(cout_units .EQ. 'percentg' .OR. cout_units .EQ. 'PERCENTG')then
			cout_units = 'percentg'
	else
		write(*,*) "cout_units", cout_units
		ierr = 1
	endif


	if(cout_format .EQ. 'text' .OR. cout_format .EQ. 'TEXT')then
			jout_format = 1
	else
			jout_format = 0
	endif


	if(ierr .EQ. 1)then
		write(6,10)
		write(6,20)
		write(6,11)
		write(6,12)
		write(6,13)
		write(6,14)
		call exit(1)
	endif
10	format('#surfseis_rspectra in= out=',/,
     +	'#           simulation_out_pointsX= simulation_out_pointsY=',/,
     +	'#           simulation_out_timesamples= simulation_out_timeskip=',/,
     +	'#           surfseis_rspectra_output_type= ',/,
     +				'surfseis_rspectra_output_format=',/,
     +	'#           surfseis_rspectra_galstog=',/,
     +	'#           surfseis_rspectra_periodtype= ',
     +				'surfseis_rspectra_period=',/,
     +	'#           surfseis_rspectra_seismogram_units= ',/,
     +	'#           surfseis_rspectra_output_units= ',/,
     +	'#           surfseis_rspectra_apply_byteswap= ',/,
     +	'#           surfseis_rspectra_apply_vel2accel= ',/,
     +	'#           surfseis_rspectra_apply_demean=',/,
     +	'#           surfseis_rspectra_apply_filter=',/,
     +	'#           surfseis_rspectra_apply_filter_lowHZ= ',/,
     +	'#           surfseis_rspectra_apply_filter_highHZ= ')
20	format(
     +	'#           surfseis_rspectra_numdampcoeffs=',/,
     +	'#           surfseis_rspectra_dampcoeffs1=',/,
     +	'#           surfseis_rspectra_dampcoeffs2=',/,
     +	'#               :        :        :      :',/,
     +	'#           surfseis_rspectra_dampcoeffs10= ',/,
     +	'# surfseis_rspectra_out_pointsXstart= ',
     +					'surfseis_rspectra_out_pointsXend=',/,
     +	'# surfseis_rspectra_out_pointsYstart= ',
     +	 				'surfseis_rspectra_out_pointsYend=',/,
     +	'# surfseis_rspectra_out_pointsXdel= ',
     +					'surfseis_rspectra_out_pointsYdel=')

11	format('#',/,
     +	'# in=                                       ',
     +				'input Surf-seis file      [none]',/,
     +	'# out=                                      ',
     +				'output spectral file      [none]',/,
     +	'# simulation_out_pointsX=                   ',
     +				'# grid points in X        [none]',/,
     +	'# simulation_out_pointsY=                   ',
     +				'# grid points in Y        [none]',/,
     +	'# simulation_out_timesamples=               ',
     +				'# time points             [none]',/,
     +	'# simulation_out_timeskip=                  ',
     +				'sample rate (sec)         [none]')
12	format(
     +	'# surfseis_rspectra_output_type=     output:',
     +				' rd,rv,prv,aa,paa,magrat    [aa]',/,
     +	'# surfseis_rspectra_output_format=          ',
     +				'binary or text          [binary]',/,
     +	'# surfseis_rspectra_galstog=                ',
     +				'gals to g (1/980)     [.0010204]',/,
     +	'# surfseis_rspectra_periodtype=             ',
     +				'SCEC or TI                [SCEC]',/,
     +	'# surfseis_rspectra_period=                 ',
     +				'output period (sec)or "all" [2.]',/,
     +	'# surfseis_rspectra_seismogram_units=       ',
     +				'                   [meterpersec]',/,
     +	'#                                           ',
     +				'vel:   "meterpersec"  "cmpersec"',/,
     +	'#                                           ',
     +				'accel: "meterpersec2" "cmpersec2"',/,
     +	'# surfseis_rspectra_output_units=           ',
     +				'                      [unitsofg]',/,
     +	'#',43x, 		'rd:    "meter"        "cm"',/,
     +	'#',43x, 		'rv,prv:"meterpersec"  "cmpersec"',/,
     +	'#',43x, 		'aa,paa:"meterpersec2" "cmpersec2"',/,
     +	'#',43x, 		'aa,paa:"unitsofg"     "percentg"',/,
     +	'# surfseis_rspectra_apply_byteswap=         ',
     +				'apply: "yes" or "no"        [no]',/,
     +	'# surfseis_rspectra_apply_vel2accel=        ',
     +				'apply: "yes" or "no"       [yes]',/,
     +	'# surfseis_rspectra_apply_demean=           ',
     +				'apply: "yes" or "no"       [yes]')
13	format(
     +	'# surfseis_rspectra_apply_filter=           ',
     +				'apply: "yes" or "no"       [yes]',/,
     +	'# surfseis_rspectra_apply_filter_lowHZ=     ',
     +				'low HZ cutoff               [0.]',/,
     +	'# surfseis_rspectra_apply_filter_highHZ=    ',
     +				'high HZ cutoff         [Nyquist]',/,
     +	'# surfseis_rspectra_numdampcoeffs=          ',
     +				'#damping coeffs              [1]',/,
     +	'# surfseis_rspectra_dampcoeffs1=            ',
     +				'damping value #1          [0.05]',/,
     +	'# surfseis_rspectra_dampcoeffs2=            ',
     +				'damping value #2          [none]',/,
     +	'#                                           ',
     +				'               :',/,
     +	'# surfseis_rspectra_dampcoeffs10=           ',
     +				'damping value #10         [none]')

14	format('#      More optional: output decimation in surface nodes:',/,
     +	'# surfseis_rspectra_out_pointsXstart=       ',
     +				'first seismogram in X        [1]',/,
     +	'# surfseis_rspectra_out_pointsXend=         ',
     +				'last  seismogram in X    [all X]',/,
     +	'# surfseis_rspectra_out_pointsYstart=       ',
     +				'first seismogram in Y        [1]',/,
     +	'# surfseis_rspectra_out_pointsYend=         ',
     +				'last  seismogram in Y    [all Y]',/,
     +	'# surfseis_rspectra_out_pointsXdel=         ',
     +				'increment in X               [1]',/,
     +	'# surfseis_rspectra_out_pointsYdel=         ',
     +				'increment in Y               [1]')


	if(kout_units .LE. 0)then
		if(cout_choice .EQ. 'rd')then
			cout_units = 'meter'
		elseif(cout_choice .EQ. 'rv')then
			cout_units = 'meterpersec'
		elseif(cout_choice .EQ. 'prv')then
			cout_units = 'meterpersec'
		elseif(cout_choice .EQ. 'aa')then
			cout_units = 'unitsofg'
		elseif(cout_choice .EQ. 'paa')then
			cout_units = 'unitsofg'
		elseif(cout_choice .EQ. 'magrat')then
			cout_units = 'unitsofg'
		else
			cout_units = 'unitsofg'
		endif
	endif



c determine output node totals
	koutnx = 0
	do i = koutx0,koutx1,koutdx
		koutnx = koutnx + 1
	enddo

	koutny = 0
	do i = kouty0,kouty1,koutdy
		koutny = koutny + 1
	enddo


	return
	end


c==================================================================
	subroutine runtime_doc(nx,ny,npts,dt, fac,ndamp,damp,maxdamp,
     +		cperiod, capply_demean,capply_filter,capply_vel2accel,
     +		filter_lowHZ,filter_highHZ,period_out,capply_byteswap,
     +		cout_choice,cin_units,cout_units,jout_format,jout_numperiods)
	implicit none
	character*120	cperiod, cout_choice,cin_units,cout_units
	character*120	capply_demean,capply_filter,capply_vel2accel
	character*120	capply_byteswap
	integer		nx,ny,npts,ndamp,maxdamp
	integer		jout_format,jout_numperiods
	real		dt,fac,damp(maxdamp),filter_lowHZ,filter_highHZ
	real		period_out

	integer		i
	character*20	cline
	character*48	cout_option

	common	/knodes/koutnx,koutx0,koutx1,koutdx, koutx,
     +			koutny,kouty0,kouty1,koutdy, kouty
	integer		koutnx,koutx0,koutx1,koutdx, koutx,
     +			koutny,kouty0,kouty1,koutdy, kouty

	cline = '#surfseis_rspectra> '

	write(6,100)cline,ny
	write(6,110)cline,nx
	write(6,120)cline,npts
	write(6,130)cline,dt

	write(6,150)cline,koutx0,koutx1,koutdx
	write(6,160)cline,kouty0,kouty1,koutdy
	write(6,170)cline,koutnx,koutny

	if(cout_choice .EQ. 'rd')then
		cout_option = 'real component of displacement spectrum.'
	elseif(cout_choice .EQ. 'rv')then
		cout_option = 'real component of velocity spectrum.'
	elseif(cout_choice .EQ. 'prv')then
		cout_option = 'imaginary component of velocity spectrum.'
	elseif(cout_choice .EQ. 'aa')then
		cout_option = 'real component of acceleration spectrum.'
	elseif(cout_choice .EQ. 'paa')then
		cout_option = 'imaginary component of acceleration spectrum.'
	elseif(cout_choice .EQ. 'magrat')then
		cout_option = 'real comp. of accel. normalized to max accel.'
	endif
	write(6,200)cline,cout_option
	write(6,210)cline,cout_units
	if(jout_format .EQ. 0 .AND. jout_numperiods .EQ. 1)write(6,220)cline
	if(jout_format .EQ. 0 .AND. jout_numperiods .EQ. 0)write(6,230)cline
	if(jout_format .EQ. 1 .AND. jout_numperiods .EQ. 1)write(6,240)cline
	if(jout_format .EQ. 1 .AND. jout_numperiods .EQ. 0)write(6,250)cline

	write(6,300)cline,cperiod
	if(jout_numperiods .EQ. 1)write(6,310)cline,period_out
	if(jout_numperiods .EQ. 0)write(6,320)cline

	write(6,390)cline,cin_units
	if(capply_byteswap .EQ. 'yes')write(6,395)cline
	if(capply_byteswap .NE. 'yes')write(6,396)cline
	if(capply_demean .EQ. 'yes')write(6,400)cline
	if(capply_demean .NE. 'yes')write(6,410)cline
	if(capply_vel2accel .EQ. 'yes')write(6,420)cline
	if(capply_vel2accel .NE. 'yes')write(6,430)cline
	if(capply_filter .NE. 'yes')write(6,440)cline
	if(capply_filter .EQ. 'yes')then
		write(6,450)cline
		write(6,460)cline,filter_lowHZ
		write(6,470)cline,filter_highHZ
	endif

	write(6,500)cline,fac
	write(6,600)cline,ndamp
	do i=1,ndamp
		write(6,610)cline,i,damp(i)
	enddo


100	format(a20,' outer dimension of #seismograms is',i7)
110	format(a20,' inner dimension of #seismograms is',i7)
120	format(a20,' #time samples per seismogram    is',i7)
130	format(a20,' time sampling interval (sec)    is',f7.4)

150	format(a20,' Seismograms to be used in X: ',i6,' to',i6,'; del:',i6)
160	format(a20,' Seismograms to be used in Y: ',i6,' to',i6,'; del:',i6)
170	format(a20,' Output #Seismograms in X, Y: ',i6,',  ',i6)

200	format(a20,' Output will be ',a48)
210	format(a20,' Output units:  ',a48)
220	format(a20,' Output format: ','binary file of one period')
230	format(a20,' Output format: ','binary file of all periods')
240	format(a20,' Output format: ','text file of one period')
250	format(a20,' Output format: ','text file of all periods')

300	format(a20,' Period lookup table (TI or SCEC) is ',a8)
310	format(a20,' Requested period is ',f10.5,' sec.')
320	format(a20,' Requested period is table of all periods.')

390	format(a20,' Incoming seismograms have units specified as ',a16)
395	format(a20,' Preprocessing 0) apply seismogram byteswap: YES') 
396	format(a20,' Preprocessing 0) apply seismogram byteswap: NO') 
400	format(a20,' Preprocessing 1) apply seismogram DEMEAN: YES') 
410	format(a20,' Preprocessing 1) apply seismogram DEMEAN: NO') 
420	format(a20,' Preprocessing 2) apply seismogram vel->accel: YES') 
430	format(a20,' Preprocessing 2) apply seismogram vel->accel: NO') 
440	format(a20,' Preprocessing 3) apply seismogram filter: NO') 
450	format(a20,' Preprocessing 3) apply seismogram filter: YES') 
460	format(a20,' Preprocessing 3) filter low  frequency (Hz):',f8.4)
470	format(a20,' Preprocessing 3) filter high frequency (Hz):',f8.4)

500	format(a20,' Factor to convert gals to g:',f10.7)
600	format(a20,' Number of user-requested damping coefficients:',i5)
610	format(a20,'        Damping value #',i2,' is',f10.7)


	return
	end

c==================================================================
	subroutine check_units_vel2accel(cunits,cvel2accel)
	implicit none
	character*120	cunits,cvel2accel
	integer		ierr

	ierr = 0

	if    (cunits .EQ. 'meterpersec'  .AND. cvel2accel .EQ. 'yes')then
			ierr = 0
	elseif(cunits .EQ. 'cmpersec'     .AND. cvel2accel .EQ. 'yes')then
			ierr = 0
	elseif(cunits .EQ. 'meterpersec2' .AND. cvel2accel .EQ. 'no')then
			ierr = 0
	elseif(cunits .EQ. 'cmpersec2'    .AND. cvel2accel .EQ. 'no')then
			ierr = 0

	elseif(cunits .EQ. 'meterpersec'  .AND. cvel2accel .EQ. 'no')then
			ierr = 1
	elseif(cunits .EQ. 'cmpersec'     .AND. cvel2accel .EQ. 'no')then
			ierr = 1
	elseif(cunits .EQ. 'meterpersec2' .AND. cvel2accel .EQ. 'yes')then
			ierr = 2
	elseif(cunits .EQ. 'cmpersec2'    .AND. cvel2accel .EQ. 'yes')then
			ierr = 2
	endif

	if(ierr .EQ. 1)then
		write(6,100)
		cvel2accel = 'yes'
	elseif(ierr .EQ. 2)then
		write(6,100)
		cvel2accel = 'no'
	endif

100	format('#surfseis_rspectra> Specified units are VELOCITY but ',
     +						'vel->accel NOT requested;',/,
     +	'#surfseis_rspectra>    vel->accel WILL BE performed (override).')

200	format('#surfseis_rspectra> Specified units are ACCELERATION but ',
     +						'vel->accel REQUESTED;',/,
     +	'#surfseis_rspectra>    vel->accel WILL NOT be performed (override).')

	return
	end
c==================================================================
	subroutine open_files(ofile,npts,koutnx,jout_format,
     +				jout_numperiods,nt,ierr,output_opt)
	implicit none
	character*256	ofile
	integer npts,koutnx,jout_format,jout_numperiods,ierr,nt,output_opt
	logical exist

	common		/F90_REC_COUNTING/REC_COUNT
	integer				  REC_COUNT
c value of REC_COUNT defined at top of main.

c output file
c   binary - one period (koutnx X koutny)
	if(jout_format .EQ. 0 .AND. jout_numperiods .EQ. 1)then
		open(22,file=ofile,err=910,access='direct',form='unformatted',
     +			status='replace', recl=REC_COUNT*koutnx)

c   binary all periods (nt X koutnx X koutny)
	elseif(jout_format .EQ. 0 .AND. jout_numperiods .EQ. 0)then
c		open(22,file=ofile,err=910,access='direct',form='unformatted',
c     +			status='replace', recl=REC_COUNT*nt)
c  If we're not using pipe forwarding or traditional, don't do any output
		if(output_opt.EQ.1) then
c  Using pipe forwarding, so use stream access
			write(*,*) "Using pipe forwarding."
			open(22,file=ofile,err=910,access='stream',
     +			form='unformatted',position='append')
		elseif(output_opt.EQ.0) then
			inquire(file=ofile, exist=exist)
                        if (exist) then
				open(22,file=ofile,err=910,access='direct',
     +				status='old', form='unformatted',recl=REC_COUNT*2)
			else
				open(22,file=ofile,err=910,access='direct',
     +				form='unformatted',recl=REC_COUNT*2)
			endif
		endif

	elseif(jout_format .EQ. 1)then
		open(22,file=ofile,err=910)

	endif


	return

900	ierr=1
	return

910	ierr=2
	return
	end

c==================================================================
	subroutine get_periodout(t,nt,jout_numperiods,period_out,jperiod_out,
     +							jperiod0,jperiod1)
	implicit none
	integer	nt,jout_numperiods,jperiod_out, jperiod0,jperiod1
	real	t(nt),period_out,	del_t,del_j,ratio

	integer	i

	jperiod_out = 0

	if(jout_numperiods .EQ. 0)then
		jperiod0 = 1
		jperiod1 = nt

	elseif(jout_numperiods .EQ. 1)then
c first check if desired period is exactly on a defined value
	do i=1,nt
		if(period_out .EQ. t(i))then
			jperiod_out = i
			goto 100
		endif
	enddo

c if here, desired period is not exactly on a defined value. 
	if (period_out .GT. t(1))then
		jperiod_out = 1
		goto 100
	elseif(period_out .LT. t(nt))then
		jperiod_out = nt
		goto 100
	else
		do i=1,nt-1
		if(period_out .LT. t(i) .AND. period_out .GT. t(i+1))then
			del_t = t(i) - t(i+1)
			del_j = period_out - t(i+1)
			ratio = del_j / del_t
			if(ratio .GE. 0.50)then
				jperiod_out = i
				goto100
			elseif(ratio .LT. 0.50) then
				jperiod_out = i+1
				goto100
			endif
		endif
		enddo

	endif

100	continue
		jperiod0 = jperiod_out
		jperiod1 = jperiod_out

	write(6,200)jperiod_out,t(jperiod_out)
200	format(
     +	'#surfseis_rspectra> Actual output period is index ',i3,
     +				'; period (sec)',f10.5)

	endif
	return
	end

c==================================================================
	subroutine seismogram_byteswap(a,npts)
	implicit none

	integer	npts,	i
	real	a(npts)
	real	xxin,xout
	integer	kkin,kout
	equivalence	(xxin,kkin),(xout,kout)

	do i=1,npts
		xxin = a(i)
		call byte_swap4(kkin,kout)
		a(i) = xout
	enddo

	return
	end

        subroutine byte_swap4(cin4,cout4)
        character*1 cin4(4),cout4(4)

        cout4(1) = cin4(4)
        cout4(2) = cin4(3)
        cout4(3) = cin4(2)
        cout4(4) = cin4(1)

        return
        end

c==================================================================
	subroutine seismogram_mks2cgs(a,npts)
	implicit none

	integer	npts,	i
	real	a(npts)

	do i=1,npts
		a(i) = a(i) * 100.
	enddo

	return
	end

c==================================================================
	subroutine seismogram_demean(a,npts,work)
	implicit none

	integer	npts,	i
	real	a(npts),work(npts),xmean,xnpts

	do i=1,npts
		work(i) = a(i)
	enddo
	do i=1,npts
		a(i) = 0.
	enddo

	xmean = 0.	
	xnpts = float(npts)
	do i=1,npts
		xmean = xmean + work(i) / xnpts
	enddo

	do i=1,npts
		a(i) = work(i) - xmean
	enddo

	return
	end

c==================================================================
	subroutine seismogram_vel2accel(a,npts,work,dt)
	implicit none

	integer	npts,	i
	real	a(npts),work(npts),dt

	do i=1,npts
		work(i) = a(i)
	enddo
	do i=1,npts
		a(i) = 0.
	enddo

	a(1) = 0.
	do i=2,npts
		a(i) = (work(i) - work(i-1)) / dt
	enddo

	return
	end


c==================================================================
	subroutine seismogram_filter(a,npts,work,dt,
     +					filter_lowHZ,filter_highHZ)
	implicit none

	integer	npts,	i
	real	a(npts),work(npts),dt,filter_lowHZ,filter_highHZ
	real	flow,fhigh
	integer	nplo,nphi,phase

	do i=1,npts
		work(i) = 0.
	enddo

	flow  = filter_lowHZ
	fhigh = filter_highHZ
	nplo  = 6
	nphi  = 6

c =0 zero-phase, =1 minimum phase
	phase = 1

	call stanford_bandpass(a,work,npts,dt, flow, fhigh, nplo, nphi, phase)

	do i=1,npts
		a(i) = work(i)
	enddo

	return
	end


c==================================================================
c note: stanford bandpass code requires incoming data and work
c		arrays to have two extra samples (npts2 = npts+2)
c
	subroutine seismogram_stanford(a,npts,npts2,work,dt,
     +					filter_lowHZ,filter_highHZ)
	implicit none

	integer	npts,npts2,	i
	real	a(npts2),work(npts2),dt,filter_lowHZ,filter_highHZ
	real	flow,fhigh
	integer	nplo,nphi,phase


	do i=1,npts2
		work(i) = 0.
	enddo

	flow  = filter_lowHZ
	fhigh = filter_highHZ
	nplo  = 6
	nphi  = 6

c =0 zero-phase, =1 minimum phase
	phase = 1

	call stanford_bandpass(a,work,npts,dt, flow, fhigh, nplo, nphi, phase)

c transfer only npts, not npts2:
	do i=1,npts
		a(i) = work(i)
	enddo

	return
	end

	integer function strlen(str)
	
	integer i
	character str*(*)
	i = len(str)
	do while (str(i:i) .EQ. ' ')
		i = i - 1
	enddo
	strlen = i
	write(*,*) str, " is length ", strlen
	return
	end
