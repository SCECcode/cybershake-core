C   Efecan added ntiskp_sgt parameter to control SGT strain outputs on Oct 4, 2012
C
C   'io.f' READS AND WRITES SIMULATION PARAMETERS AND VARIABLES


      subroutine open3dp(io,comm,fhx,fhy,fhz,sxgro_timestep,sygro_timestep,szgro_timestep,
     + 				fhx_md5,fhy_md5,fhz_md5,sxgro_md5,sygro_md5,szgro_md5,imd5,ifault,is_z,filen)

CCC   ROUTINE TO OPEN AND CLOSE MPI FILES

C     IO   UNIT #       (INTEGER)(SENT)
C     COMM COMMUNICATOR (INTEGER)(SENT)
C     FH*  FILE POINTER (INTEGER)(RETURN)

      include 'mpif.h'

      integer :: io,comm,fhx,fhy,fhz,fhx_md5,fhy_md5,fhz_md5,err,ifault,is_z,filen

      character (len=filen) :: sxgro_timestep,sygro_timestep,szgro_timestep
      character (len=filen) :: sxgro_md5,sygro_md5,szgro_md5

C     OUTPUT: SEISMOGRAMS

      if (io.eq.1) then

        call MPI_FILE_OPEN(comm,sxgro_timestep,MPI_MODE_WRONLY +
     +                   MPI_MODE_CREATE,MPI_INFO_NULL,fhx,err)

        call MPI_FILE_OPEN(comm,sygro_timestep,MPI_MODE_WRONLY +
     +                   MPI_MODE_CREATE,MPI_INFO_NULL,fhy,err)

        if (is_z .eq. 1) then
        call MPI_FILE_OPEN(comm,szgro_timestep,MPI_MODE_WRONLY +
     +                   MPI_MODE_CREATE,MPI_INFO_NULL,fhz,err)
        endif
        if (imd5==1) then
          call MPI_FILE_OPEN(comm,sxgro_md5,MPI_MODE_WRONLY +
     +                   MPI_MODE_CREATE,MPI_INFO_NULL,fhx_md5,err)

          call MPI_FILE_OPEN(comm,sygro_md5,MPI_MODE_WRONLY +
     +                   MPI_MODE_CREATE,MPI_INFO_NULL,fhy_md5,err)

          call MPI_FILE_OPEN(comm,szgro_md5,MPI_MODE_WRONLY +
     +                   MPI_MODE_CREATE,MPI_INFO_NULL,fhz_md5,err)
        end if
      else

C     CLOSE ALL FILES

        call MPI_FILE_CLOSE(fhx,err)
        call MPI_FILE_CLOSE(fhy,err)
        if (is_z .eq. 1) then
        call MPI_FILE_CLOSE(fhz,err) 
        endif

        if (imd5==1) then
          call MPI_FILE_CLOSE(fhx_md5,err)
          call MPI_FILE_CLOSE(fhy_md5,err)
          call MPI_FILE_CLOSE(fhz_md5,err)
        end if 

      end if
      
      return
      end

      subroutine open3dp2(io,comm,fhx2,fhy2,fhz2,sxgro2_timestep,sygro2_timestep,
     +                    szgro2_timestep,fhx2_md5,fhy2_md5,fhz2_md5,sxgro2_md5,
     +                    sygro2_md5,szgro2_md5,imd5,ifault,is_z,filen)

CCC   ROUTINE TO OPEN AND CLOSE MPI FILES

C     IO   UNIT #       (INTEGER)(SENT)
C     COMM COMMUNICATOR (INTEGER)(SENT)
C     FH*  FILE POINTER (INTEGER)(RETURN)

      include 'mpif.h'

      integer :: io,comm,fhx2,fhy2,fhz2,fhx2_md5,fhy2_md5,fhz2_md5,err,ifault,is_z,filen

      character (len=filen) :: sxgro2_timestep,sygro2_timestep,szgro2_timestep
      character (len=filen) :: sxgro2_md5,sygro2_md5,szgro2_md5

C     OUTPUT: SEISMOGRAMS

      if (io.eq.1) then

        call MPI_FILE_OPEN(comm,sxgro2_timestep,MPI_MODE_WRONLY +
     +                   MPI_MODE_CREATE,MPI_INFO_NULL,fhx2,err)

        call MPI_FILE_OPEN(comm,sygro2_timestep,MPI_MODE_WRONLY +
     +                   MPI_MODE_CREATE,MPI_INFO_NULL,fhy2,err)

        if (is_z .eq. 1) then
        call MPI_FILE_OPEN(comm,szgro2_timestep,MPI_MODE_WRONLY +
     +                   MPI_MODE_CREATE,MPI_INFO_NULL,fhz2,err)
        endif

        if (imd5==1) then
          call MPI_FILE_OPEN(comm,sxgro2_md5,MPI_MODE_WRONLY +
     +                   MPI_MODE_CREATE,MPI_INFO_NULL,fhx2_md5,err)

          call MPI_FILE_OPEN(comm,sygro2_md5,MPI_MODE_WRONLY +
     +                   MPI_MODE_CREATE,MPI_INFO_NULL,fhy2_md5,err)

          call MPI_FILE_OPEN(comm,szgro2_md5,MPI_MODE_WRONLY +
     +                   MPI_MODE_CREATE,MPI_INFO_NULL,fhz2_md5,err)
        end if
      else

C     CLOSE ALL FILES

        call MPI_FILE_CLOSE(fhx2,err)
        call MPI_FILE_CLOSE(fhy2,err)
        if (is_z .eq. 1) then
        call MPI_FILE_CLOSE(fhz2,err)
        endif

        if (imd5==1) then
          call MPI_FILE_CLOSE(fhx2_md5,err)
          call MPI_FILE_CLOSE(fhy2_md5,err)
          call MPI_FILE_CLOSE(fhz2_md5,err)
        end if

      end if
      
      return
      end

      subroutine openpar(io)

CCC   ROUTINE TO OPEN FIXED NAME PARAMETER FILE

C     IO   UNIT #       (INTEGER)(SENT)

      include 'set_names.h'

      if (io.eq.1) then

        open(22,file=c_IN3D)
        rewind(22)

      else
   
        close(22)

      end if

      return 
      end


      subroutine open3d(io,nve,chkp,chkj,insrc,invel,filen,ifault,nxt,nst)

CCC   ROUTINE TO OPEN AND CLOSE FILES

C     IO       UNIT NUMBER              (INTEGER)  (SENT)
C     CHKP     PARAMETER CHECK FILE     (CHARACTER)(SENT)
C     CHKP     JOB VALUE CHECK FILE     (CHARACTER)(SENT)
C     INSRC    SOURCE FILE              (CHARACTER)(SENT)
C     INPAR    MEDIA FILE               (CHARACTER)(SENT)
C     FILEN    MAX LENGTH OF FILENAME   (INTEGER)  (SENT)

      integer :: filen,nve,nx,nyt,nst
      integer :: bsize,bsize2
      integer :: varnum(8)
      integer :: varnum2(6)
      character (len=filen) :: chkp,chkj,insrc,invel
      bsize=4*size(varnum)
      bsize2=4*size(varnum2)
      inquire(iolength=bsize)varnum
      inquire(iolength=bsize2)varnum2
      if (io.eq.1) then

C     INPUT: SIMULATION SPACE, SOURCE, STRUCTURE

        if (ifault == 3 .or. ifault==4) then
           open(33, file=insrc,form='formatted')
c          open(33, file=insrc,access='direct',form='unformatted',recl=bsize2)
c          print *, 'bsize2=',bsize2
          rewind(33)
        else

         if (ifault.eq.1) then      ! Binary format
           open(25,file=insrc,form='unformatted', access='direct',
     +            recl=bsize*(6*nst+3)/8)   
           rewind(25)
         else if (ifault.eq.0) then ! ASCII format
           open(25,file=insrc)
         end if
        end if

        if(nve==1) then
          open(27,file=invel,access='direct',form='unformatted',
     +         recl=bsize*nxt)
          print *, 'bsize=',bsize
        else
c          open(27,file=invel,access='direct',form='unformatted',recl=bsize2*nx*nyt)
          open(27,file=invel,access='direct',form='unformatted',recl=bsize2*nxt)
        end if

C     OUTPUT: SIMULATION PARAMETERS, WAVEFIELD

        open(21,file=chkp)
        open(24,file=chkj)

C     ENSURE SEQUENTIAL FILES START AT RECORD 1

        rewind(21)
        rewind(24)

      else


C     CLOSE ALL FILES

        close(21)
        close(27)
        close(24)
        if (ifault==3 .or. ifault==4) then
          close(33)
        else
          close(25)
        endif

      end if

      return
      end

      subroutine rddt3d(igreen,tmax,dh,dt,npc,nd,arbc,nsrc,nst,nx,ny,nz,
     +           npx,npy,npz,nbgx,nedx,nskpx,nbgy,nedy,nskpy,
     +           nbgz,nedz,nskpz,ntiskp,nbgx2,nedx2,nskpx2,nbgy2,
     +           nedy2,nskpy2,nbgz2,nedz2,nskpz2,ntiskp2,ntiskp_sgt,nve,
     +           mu_ss,mu_dd,fl,
     +           fh,fp,chkp,chkj,insrc,invel,insgt,sxgro,sygro,szgro,
     +           sxgro2,sygro2,szgro2,sgtgro,read_step,write_step,write_step2,
     +           sgsn,filen,pht1,
     +           ifault,checkpoint,isfcvlm,imd5,ivelocity,
     +           media,io_opt,iperf,idyna,nvar,iost,partdeg,SoCalQ)

CCC   ROUTINE TO READ MODEL AND SOURCE PARAMETERS

C     igreen SGT INPUT COLUMN INDEX                (INTEGER)(RETURNED)
C     TMAX   PROPAGATION TIME                      (REAL)   (RETURNED)
C     DH     SPACIAL DISCRETIZATION STEP           (REAL)   (RETURNED)
C     DT     TIME DISCRETIZATION STEP              (REAL)   (RETURNED)
C     NPC    CHOOCE PML OR CERJAN                  (INTEGER)(RETURNED)
C     ND     NUMBER OF PML NODES                   (INTEGER)(RETURNED)
C     ARBC   ARBITRARY DAMPING COEFFICIENT - PML   (REAL)   (RETURNED)
C     NSRC   NUMBER OF NODES ON FAULT              (INTEGER)(RETURNED)
C     NST    LENGTH OF SOURCE TIME FUNCTION (NODES)(INTEGER)(RETURNED)
C     NX     X MODEL DIMENSION                     (INTEGER)(RETURNED)
C     NY     Y MODEL DIMENSION                     (INTEGER)(RETURNED)
C     NZ     Z MODEL DIMENSION                     (INTEGER)(RETURNED)
C     NPX    NUMBER OF PROCS IN X DIRECTION        (INTEGER)(RETURNED)
C     NPY    NUMBER OF PROCS IN Y DIRECTION        (INTEGER)(RETURNED)
C     NPZ    NUMBER OF PROCS IN Z DIRECTION        (INTEGER)(RETURNED)
C     IFAULT     FAULT OR INITIAL STRESS MODEL     (INTEGER)(RETURNED) 
C     CHECKPOINT CHECKPOINT AT STEPS               (INTEGER)(RETURNED)
C     ISCFVLM    SURFACE OR VOLUME OUTPUT          (INTEGER)(RETURNED)
C     IMD5       MD5 OUTPUT                        (INTEGER)(RETURNED)
C     IVELOCITY  OUTPUT ACCUMULATION               (INTEGER)(RETURNED)
C     MEDIA      MEDIA SELECTION                   (INTEGER)(RETURNED)
C     IO_OPT     OUTPUT DATA ON/OFF                (INTEGER)(RETURNED)
C     INI_HOMO   FAST INIT                         (INTEGER)(RETURNED)
C     IPERF      PERF TIMING MEASUREMENT           (INTEGER)(RETURNED)
C     IDYNA      DYNAMIC RUPTURE MODEL             (INTEGER)(RETURNED)
C     SOCALQ     SO CAL VP-VS Q RELATIONSHIP       (INTEGER)(RETURNED)
C     NVE    CHOOSE VISCO OR ELASTIC FD SCHEME     (INTEGER)(RETURNED)
C     MU_SS  STATIC FRICTION COEF                  (REAL)   (RETURNED)
C     MU_DD  DYNAMIC FRICTION COEF                 (REAL)   (RETURNED)
C     FL     Q FREQUENCY LOW                       (REAL)   (RETURNED)
C     FH     Q FREQUENCY HIGH                      (REAL)   (RETURNED)
C     FP     Q FREQUENCY CENTER                    (REAL)   (RETURNED)
C     READ_STEP                                    (INTEGER)(RETURNED)
C     WRITE_STEP                                   (INTEGER)(RETURNED)
C     WRITE_STEP2                                  (INTEGER)(RETURNED)
C     NBGX   FIRST LINE IN X DIR CONTAIN RECEIVERS (INTEGER)(RETURNED)
C     NEDX   LAST LINE IN X DIR CONTAIN RECEIVERS  (INTEGER)(RETURNED)
C     NSKPX  LOOP INCREMENT IN X DIR FOR RECEIVERS (INTEGER)(RETURNED)
C     NBGY   FIRST LINE IN Y DIR CONTAIN RECEIVERS (INTEGER)(RETURNED)
C     NEDY   LAST LINE IN Y DIR CONTAIN RECEIVERS  (INTEGER)(RETURNED)
C     NSKPY  LOOP INCREMENT IN Y DIR FOR RECEIVERS (INTEGER)(RETURNED)
C     NBGZ   FIRST LINE IN Z DIR CONTAIN RECEIVERS (INTEGER)(RETURNED)
C     NEDZ   LAST LINE IN Z DIR CONTAIN RECEIVERS  (INTEGER)(RETURNED)
C     NSKPZ  LOOP INCREMENT IN Z DIR FOR RECEIVERS (INTEGER)(RETURNED)
C     NTISKP LOOP INCREMENT FOR TIME SAMPLES       (INTEGER)(RETURNED)
C     NBGX2  FIRST LINE IN X DIR CONTAIN RECEIVERS (INTEGER)(RETURNED)
C     NEDX2  LAST LINE IN X DIR CONTAIN RECEIVERS  (INTEGER)(RETURNED)
C     NSKPX2 LOOP INCREMENT IN X DIR FOR RECEIVERS (INTEGER)(RETURNED)
C     NBGY2  FIRST LINE IN Y DIR CONTAIN RECEIVERS (INTEGER)(RETURNED)
C     NEDY2  LAST LINE IN Y DIR CONTAIN RECEIVERS  (INTEGER)(RETURNED)
C     NSKPY2 LOOP INCREMENT IN Y DIR FOR RECEIVERS (INTEGER)(RETURNED)
C     NBGZ2  FIRST LINE IN Z DIR CONTAIN RECEIVERS (INTEGER)(RETURNED)
C     NEDZ2  LAST LINE IN Z DIR CONTAIN RECEIVERS  (INTEGER)(RETURNED)
C     NSKPZ2 LOOP INCREMENT IN Z DIR FOR RECEIVERS (INTEGER)(RETURNED)
C     NTISKP2 LOOP INCREMENT FOR TIME SAMPLES      (INTEGER)(RETURNED)
C     NTISKP_SGT LOOP INCREMENT FOR SGT STRAINS    (INTEGER)(RETURNED)
C     CHKP   CHECK PARAMETERS FILE                 (CHAR)   (RETURNED)
C     CHKJ   CHECK JOB VALUES FILE                 (CHAR)   (RETURNED)
C     INSRC  INPUT MT BASED FAULT FILE             (CHAR)   (RETURNED)
C     INVEL  INPUT ASCII LIST FORM VELOCITY FILE   (CHAR)   (RETURNED)
C     INSGT  INPUT ASCII STATION LIST FOR SGT      (CHAR)   (RETURNED)
C     SXGRO  REGULAR-GRID SEISMOGRAM OUTPUT FILE   (CHAR)   (RETURNED)
C     SYGRO  REGULAR-GRID SEISMOGRAM OUTPUT FILE   (CHAR)   (RETURNED)
C     SZGRO  REGULAR-GRID SEISMOGRAM OUTPUT FILE   (CHAR)   (RETURNED)
C     SXGRO2 REGULAR-GRID SEISMOGRAM OUTPUT FILE   (CHAR)   (RETURNED)
C     SYGRO2 REGULAR-GRID SEISMOGRAM OUTPUT FILE   (CHAR)   (RETURNED)
C     SZGRO2 REGULAR-GRID SEISMOGRAM OUTPUT FILE   (CHAR)   (RETURNED)
C     SGSN OUTPUT FILE OF ALL OTHER 19 VARIABLES   (CHAR)   (RETURNED)
      use parstat
      integer :: filen,read_step, write_step, write_step2,igreen, i
      real :: mu_ss, mu_dd,dt,tmax
      character (len=filen) :: chkp,chkj,insrc,invel,insgt,sgtgro
      character (len=filen) :: sxgro,sygro,szgro,sxgro2,sygro2,szgro2
      character (len=filen) :: sgsn
      integer :: checkpoint,ifault,isfcvlm,imd5,ivelocity
      integer :: media,io_opt,iperf,idyna,SoCalQ,nvar,iost,partdeg
      logical :: file_exists
      real :: lonsgt, latsgt, depsgt
      integer (kind=8) :: sgt_index
      character (len=1024) :: sgt_buf

      include 'set_names.h'

      read(22,*) igreen
      read(22,*) tmax
      read(22,*) dh
      read(22,*) dt
      read(22,*) npc

      read(22,*) nd
      read(22,*) arbc
      read(22,*) pht1

      read(22,*) nsrc
      read(22,*) nst

      read(22,*) nx
      read(22,*) ny
      read(22,*) nz

      read(22,*) npx
      read(22,*) npy
      read(22,*) npz

      read(22,*) ifault
      read(22,*) checkpoint
      read(22,*) isfcvlm
      read(22,*) imd5
      read(22,*) ivelocity
      read(22,*) media
      read(22,*) nvar
      read(22,*) iost
      read(22,*) partdeg
      read(22,*) io_opt
      read(22,*) iperf
      read(22,*) idyna
      read(22,*) SoCalQ

      read(22,*) nve

      read(22,*) mu_ss
      read(22,*) mu_dd

      read(22,*) fl
      read(22,*) fh
      read(22,*) fp


      read(22,*) read_step
      read(22,*) write_step
      read(22,*) write_step2
                
      read(22,*) nbgx
      read(22,*) nedx
      read(22,*) nskpx
      read(22,*) nbgy
      read(22,*) nedy
      read(22,*) nskpy
      read(22,*) nbgz
      read(22,*) nedz
      read(22,*) nskpz

      read(22,*) ntiskp

      read(22,*) nbgx2
      read(22,*) nedx2
      read(22,*) nskpx2
      read(22,*) nbgy2
      read(22,*) nedy2
      read(22,*) nskpy2
      read(22,*) nbgz2
      read(22,*) nedz2
      read(22,*) nskpz2

      read(22,*) ntiskp2
      read(22,*) ntiskp_sgt

      read(22,*) chkp
      read(22,*) chkj
      read(22,*) insrc
      read(22,*) invel
      read(22,*) insgt

      read(22,*) sxgro
      read(22,*) sygro
      read(22,*) szgro

      read(22,*) sxgro2
      read(22,*) sygro2
      read(22,*) szgro2

      read(22,*) sgtgro

      read(22,*) sgsn

cSGT -------start---------------------------------------
      sgt_io_out=0
      if (igreen .ne. -1) then
        INQUIRE(FILE=insgt, EXIST=file_exists)
        if (file_exists) then
            open(20,file=insgt)
            sgt_buf='#'
            do while (sgt_buf(1:1)=='#')
                read(20,*) sgt_buf
            end do
            read(sgt_buf,*) sgt_numsta
            allocate(sgt_sta_coord(sgt_numsta,3))
            do i=1,sgt_numsta
                read(20,*) sgt_sta_coord(i,1),sgt_sta_coord(i,2),
     +              sgt_sta_coord(i,3)
            end do
            close(20)
            sgt_io_out=1;
            write(*,*) 'Use SGT output mode, sgt_numsta=',sgt_numsta
        else
            write(*,*) 'Station list does not exist. Use normal output mode'
        endif
      endif
cSGT --------end----------------------------------------

      call override_names(chkp,chkj,insrc,invel,sxgro,sygro,szgro,
     +                                  sxgro2,sygro2,szgro2,filen)

!     print *, mu_s, mu_d, fl, fp,nsrc,npsrc
!     print *, ifault, nve, szgro2
  
c      write(*,*) tmax,dh,dt,npc,nd,arbc,nsrc,nst,nx,ny,nz,
c     +           npx,npy,npz,nbgx,nedx,nskpx,nbgy,nedy,nskpy,
c     +           nbgz,nedz,nskpz,ntiskp,nbgx2,nedx2,nskpx2,nbgy2,
c     +           nedy2,nskpy2,nbgz2,nedz2,nskpz2,ntiskp2,nve,fl,
c     +           fh,fp,chkp,chkj,insrc,invel,sxgro,sygro,szgro,
c     +           sxgro2,sygro2,szgro2,filen
      return
      end

      subroutine override_names(chkp,chkj,insrc,invel,sxgro,sygro,szgro,
     +                                  sxgro2,sygro2,szgro2,filen)
      integer :: filen
      character (len=filen) :: chkp,chkj,insrc,invel
      character (len=filen) :: sxgro,sygro,szgro,sxgro2,sygro2,szgro2

      include 'set_names.h'

      if(num_commandlineargs .GE.  2)insrc  = c_SOURCE
      if(num_commandlineargs .GE.  3)invel  = c_MEDIA
      if(num_commandlineargs .GE.  4)chkp   = c_CHKMOD
      if(num_commandlineargs .GE.  5)chkj   = c_CHKJOB
      if(num_commandlineargs .GE.  6)sxgro  = c_SSX3D
      if(num_commandlineargs .GE.  7)sygro  = c_SSY3D
      if(num_commandlineargs .GE.  8)szgro  = c_SSZ3D
      if(num_commandlineargs .GE.  9)sxgro2 = c_SSX3D2
      if(num_commandlineargs .GE. 10)sygro2 = c_SSY3D2
      if(num_commandlineargs .GE. 11)szgro2 = c_SSZ3D2

      return
      end

      subroutine wrpm3d(igreen,nxt,nyt,nzt,nt,dh,dt,npc,
     +  arbc,ntiskp,nbgx,nedx,nskpx,nbgy,nedy,nskpy,
     +  nbgz,nedz,nskpz,vpe,vse,den,nve,fl,fh,fp)

CCC   WRITE INPUT AND COMPUTED SIMULATION PARAMETERS TO CHECK FILE

      integer, parameter :: nl=2
      real, dimension(nl) :: vpe,vse,den
      integer igreen
      write(21,*) 'igreen = ', igreen

      write(21,70) vpe(2)*dt/dh
      write(21,55) nxt
      write(21,56) nyt
      write(21,57) nzt
      write(21,58) nt
      write(21,59) dh
      write(21,60) dt
      write(21,62) arbc
      write(21,84) vpe(2)
      write(21,85) vpe(1)
      write(21,86) vse(2)
      write(21,87) vse(1)
      write(21,88) den(2)
      write(21,89) den(1)
      write(21,28) ntiskp
      write(21,65) npc
      write(21,63) nve
      write(21,113) fl,fp,fh

   28 format('SKIP OF SEISMOGRAMS IN TIME (LOOP COUNTER)',I5)
   84 format('HIGHEST P-VELOCITY ENCOUNTERED',F10.2)
   85 format('LOWEST P-VELOCITY ENCOUNTERED',F10.2)
   86 format('HIGHEST S-VELOCITY ENCOUNTERED',F10.2)
   87 format('LOWEST S-VELOCITY ENCOUNTERED',F10.2)
   88 format('HIGHEST DENSITY ENCOUNTERED',F15.2)
   89 format('LOWEST  DENSITY ENCOUNTERED',F15.2)
   55 format('# OF X NODES PER PROC',I5)
   56 format('# OF Y NODES PER PROC',I5)
   57 format('# OF Z NODES PER PROC',I5)
   58 format('# OF TIME STEPS',I5)
   59 format('DISCRETIZATION IN SPACE',F10.4)
   60 format('DISCRETIZATION IN TIME',F10.4)
   65 format('ABC CONDITION, PML=1 OR CERJAN=0:',I5)
   63 format('FD SCHEME, VISCO=1 OR ELASTIC=0:',I5)
   62 format('PML REFLECTION COEFFICIENT',G13.4)
   70 format('STABILITY CRITERIA .5 > CMAX*DT/DX = ',E15.7)
  113 format('Q, FL,FP,FH:',3G14.5)

c     use the following line on platforms other than XT3
      call flush(21)
c     use the following line on XT3
c     call flush_(21)

      return
      end
