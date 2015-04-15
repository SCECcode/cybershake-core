
C     MPI-2 VERSION OF 3D, 4TH-ORDER, STAGGERED-GRID, ELASTIC
C     OR VISCOELASTIC FD WAVE PROPAGATION SCHEME. INCLUDES PML
C     ABSORBING BOUNDARY CONDITION, 3D STRUCTURE INPUT, EXTENDED
C     SOURCE FAULT INPUT AND REGULAR GRID 4D VOLUME OUTPUT

C
C     AWM INTEGRATIED VERSION
C
C     C. MARCINKOVICH AND K.B. OLSEN (2004)
C
C     Y.F. CUI (2004), MPI/MPIIO/INITIALIZATION CHANGES, ENHANCE
C                      CODE WITH CAPABILITY OF RUN ON LARGE MESH SIZE
C                      - MAINLY MEMORY IMPROVEMENTS, ADD MORE SETTINGS,
C                      WITH INPUTS FROM G. CHUKKAPALLI
C     G. ELY & B. SHKOLLER (2004), ADD RESTART/CHECKPOINT CAPABILITY
C     Y.F. HU (2004), ADD MD5 CAPABILITY
C     L. BRIEGER (2004), SIMPLIFIED WAY USING NEW MC1 RANK ID
C     Y.F. CUI & Y.F. HU (2004), INTEGRATION OF WAVE PROPAGATION FEATURES
C     Y.F. CUI & Y.F. HU (2005), INCORPORATION OF DYNAMIC RUPTURE FEATURES
C                                BASED ON K.B. OLSEN (2005)
C     Y.F. CUI (2005), MERGE OF DYNAMIC RUPTURE AND WAVE PROPAGATION CODE
C
C     L.A. DALGUER (2006), INCORPORATION OF SGSN DYNAMIC FAULT MODEL
C                              BASED ON DALGUER L.A. AND S.M. DAY (2006)
c     D. Okaya (2006), Command-line definition of I/O file names. In this
c                               main.f, only involves a call to set_names().
c                               Corresponding set_names.h now exists.
c                               "IN3D" itself can now have any name and
c                               is first command line argument.  If no CLA's
c                               provided, defaults are original names.
c                               Set filen=180 from 128.
c     D. Okaya (2006), Fixed typos in structure.f for tmpx usage for
c                               elastic (versus visco) case.  Elastic
c                               uses media records of 6 4-byte words,
c                               not 8.  But current version of subroutine
c                               ini3d() incorectly uses 8 word records.
c     Y.F. CUI & J. ZHU (2007), MODIFIED RUPTURE OUTPUT TO BE ON THE ENTIRE FAULT FOR DYNAMIC MODE
c                               ADD NEW I/O FEATURE FOR SGSN MODE (outupt 19 more variables)
c     Y.F. CUI & K.Y. LEE (2009), MERGED IN3D and SETTING into IN3D
c                                 PARAMERTER ORDER IN IN3D IS ALSO CHANGED
c     Y.F. CUI & K.Y. LEE (2009), ASYNCHRONOUS COMMUNICATION CODE
c     Y.F. CUI & K.Y. LEE (2009), 64BIT MPIIO CODE ADDED

      use parstat
      implicit none
      !include "mpif.h"

C     INTEGRATION/PERFORMANCE/EXTRA VARIABLES    

      integer :: MCW,MC1,MFS,MCR
      integer :: send_d2u_1,send_u2d_1,send_d2u_2,send_u2d_2
      integer :: recv_d2u_1,recv_u2d_1,recv_d2u_2,recv_u2d_2
      integer :: send_s2n_1,send_n2s_1,send_s2n_2,send_n2s_2
      integer :: recv_s2n_1,recv_n2s_1,recv_s2n_2,recv_n2s_2
      integer :: send_w2e_1,send_e2w_1,send_w2e_2,send_e2w_2
      integer :: recv_w2e_1,recv_e2w_1,recv_w2e_2,recv_e2w_2
      integer :: nmtype(3),ntype(3),send_offset(3),recv_offset(3)
      integer :: nd,nx,ny,nz,nxt,nyt,nzt,ndt,nt,npts
      integer :: nbgx,nedx,nskpx,nbgy,nedy,nskpy,nbgz,nedz,nskpz,ntiskp
      integer :: nbgx2,nedx2,nskpx2,nbgy2,nedy2,nskpy2,nbgz2,nedz2,nskpz2,ntiskp2
      integer :: ntiskp_sgt
      real    :: dh,arbc,taumax,taumin,saveallcheckpoints
      real    :: uav,vav,wav

      integer :: j,k,it,isfcvlm,imd5,ivelocity
      integer :: filetype, filetype3, filetype4
      integer :: xn,yn,zn
      integer :: ifault,io_opt,nbyte,mediarestart,iperf,ckp_count,idtout,nvar,iost,partdeg
c     integer, parameter :: filen=128
      integer, parameter :: filen=180
      integer(i8) :: nrec,nrec2,nrec3
      integer (kind=MPI_OFFSET_KIND) :: disp,disp2,disp3,disp4
      integer (kind=MPI_OFFSET_KIND) :: offset
      integer ( kind=MPI_OFFSET_KIND ) :: nbyte_l
      integer ( kind=MPI_OFFSET_KIND ) :: nrec_l,nrec2_l,nsrc_l,nrec3_l
      integer varnum(1)
      integer,dimension(:), allocatable :: block
      real :: mu_ss, mu_dd,tmax,dt,pi,qpinv,qsinv,tmp1,tmp2,w0,ww1,w2,vpvs
      real, dimension(:), allocatable :: bufx,bufy,bufz,bufxyz,bufx2,bufy2,bufz2
      real, dimension(:,:,:), allocatable :: sgt_buf
      real, dimension(:), allocatable :: sgt_buf_buf
      integer (kind=MPI_OFFSET_KIND) :: sgt_offset
      real(b8), dimension(11) :: dtime,dtout
      real(b8) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10=0.0,t11,t12=0.,t13=0.,t14=0.,t16=0.0,
     +          t18,t19,t20=0.,t21,t22=0.,t23=0.,xmb,xmb2,xmbt,xmb2t,t15

C added by yhu for media file read-in bug fixing

      integer :: tmp_count, tmp_count2
      real, dimension(:), allocatable :: tmpall
 
C     CHECKPOINT RESTART VARIABLES

      integer :: checkpoint, bsize, iwrite
      integer :: tmp, t0=0
      integer :: count,count2, sgt_count

      character (len=filen) :: ckpfile,mediafile

C     ANNOUNCE SELECT VARIABLES
 
      integer :: rank,size1,err
      integer, dimension(MPI_STATUS_SIZE) :: status

C     CHOOSE VISCO OR ELASTIC, PML OR CERJAN

      integer :: nve,npc

C     PROC INTER-RELATIONS FOR SWAPPING
 
      integer :: fromw,toe,frome,tow
      integer :: froms,ton,fromn,tos
      integer :: fromd,tou,fromu,tod

      integer :: swd2neu,neu2swd,swu2ned,ned2swu
      integer :: sed2nwu,nwu2sed,seu2nwd,nwd2seu
      integer :: sw2ne,ne2sw,se2nw,nw2se,su2nd,nd2su
      integer :: sd2nu,nu2sd,wu2ed,ed2wu,wd2eu,eu2wd

C     MASTER PROC, 3D SPACE AND # OF COMPUTATIONAL ZONES
 
      integer, parameter :: maxdim=3, maxzone=18
      integer :: master, transmaster, rankmc1

C     MAXIMUM LENGTH OF FILENAME (PATH AND NAME)

 
C     ALLOW PROC REORDER, DISALLOW GRID PERIODICITY

      logical :: reorder=.true.
      logical, dimension(maxdim) :: periodic=.false.


C     PROC DIMENSIONS, DECOMPOSITION AND COORDINATES IN PROC SPACE

      integer :: npx,npy,npz
      integer, dimension(maxdim) :: dims,coords

C     MEMORY INITIALIZATION, FREE SURFACE AND RECEIVER FLAGS

      integer :: fs,fw,fn,fe,fb, fsf
      integer :: color
      integer :: shape
      integer :: i      
      

C     EXTREME MEDIA VALUES AND Q BANDWIDTH DESCRIPTION
 
      real, dimension(2) :: vpe,vse,dde,tmpvpe,tmpvse,tmpdde
      real :: fl,fh,fp,maxvpe,minvpe,maxvse,minvse,maxdde,mindde

C     SOURCE PROC AND SOURCE-TIME FUNCTION

      integer :: srcproc, nsrc,nst,npsrc
      integer :: recproc, nprec 
      integer :: recproc2, nprec2
      integer :: recproc3, nprec3
      integer :: recproc4, nprec4
      integer :: sgt_recproc,sgt_nprec
      integer :: cur_step,read_step,write_step,write_step2

C     MPI-2 SEISMOGRAM FILE HANDLING AND OTHER IO FILES

      integer :: fhx_md5,fhy_md5,fhz_md5,fhx2_md5,fhy2_md5,fhz2_md5,
     $           fhx,fhy,fhz,fhx2,fhy2,fhz2,fhsgsn,sgt_mpi_file
      integer (kind=MPI_OFFSET_KIND) :: disp_md5
 
      character (len=filen) :: chkp,chkj,insrc,invel,insgt
      character (len=filen) :: sxgro,sygro,szgro, sgtgro
      character (len=filen) :: sxgro2,sygro2,szgro2,sgsn
      character (len=filen) :: sxgro_timestep,sygro_timestep,szgro_timestep
      character (len=filen) :: sgt_timestep,sgt_string_buf,sgt_filename
      character (len=filen) :: sxgro2_timestep,sygro2_timestep,szgro2_timestep
      character (len=filen) :: sgsn_timestep
      character (len=filen) :: sxgro_md5,sygro_md5,szgro_md5
      character (len=filen) :: sxgro2_md5,sygro2_md5,szgro2_md5
      character *32 :: resblock

C     FOR SGSN DYNAMIC FAULT MODEL

      real :: eta               !=0.32
      integer :: ix_d,fx_d,y_d,iz_d,fz_d,idyna
      integer :: fsr            !=0
      integer :: y_dt           !=-1

C     SGT Starts
      integer :: igreen, floatsize
      real, dimension(:,:,:), allocatable :: sg1, sg2, sg3
      real, dimension(:,:,:,:), allocatable :: sgt
      integer ix, iy, iz
C     SGT ENDS

C     SOUTHERN CALIFORNIA VP-VS Q RELATIONSHIP FLAG
      integer :: SoCalQ

C     =============END OF VARIABLE DECLARATION==============

C     =======================
C     VARIABLE INITIALIZATION
C     =======================

c  integer
      disp=0
      disp2=0
      disp3=0
      disp4=0
      offset=0
      t0=0
      count=0
      count2=0
      sgt_count=0
      master=0
      fs=0
      fw=0
      fn=0
      fe=0
      fb=0
      fsf=0
      color=0
      shape=0
      srcproc=-1
      recproc=-1
      recproc3=-1
      recproc4=-1
      y_dt=-1
      floatsize=4
      sgt_mpi_file=0

c  real
      vpe=0
      vse=0
      dde=0
      tmpvpe=0
      tmpvse=0
      tmpdde=0
      eta=0.32

C     ========================================
C     GET COMMAND LINE ARGUMENTS (FILE NAMES)
C     DEFAULT (NO ARGUMENTS) OPENS IN3D
C     ========================================

      call set_names()

      bsize=4*size(varnum)
      inquire(iolength=bsize)varnum

C     =========
C     START MPI 
C     =========

      call MPI_INIT(err)
      call MPI_COMM_DUP(MPI_COMM_WORLD,MCW,err)
      call MPI_COMM_RANK(MCW,rank,err)
      call MPI_COMM_SIZE(MCW,size1,err)

C     ====================================================
C     READ PARAMETERS, COMPUTE PROC SIZE AND OPEN IO FILES
C     ====================================================

      call MPI_BARRIER(MCW,err)  
      t1=mpi_wtime()
      if(rank==master) then

        call openpar(1)
C       call openlog(1)
        call rddt3d(igreen,tmax,dh,dt,npc,nd,arbc,nsrc,nst,nx,ny,nz,
     +       npx,npy,npz,nbgx,nedx,nskpx,nbgy,nedy,nskpy,
     +       nbgz,nedz,nskpz,ntiskp,nbgx2,nedx2,nskpx2,nbgy2,
     +       nedy2,nskpy2,nbgz2,nedz2,nskpz2,ntiskp2,ntiskp_sgt,nve,
     +       mu_ss,mu_dd,fl,
     +       fh,fp,chkp,chkj,insrc,invel,insgt,sxgro,sygro,szgro,
     +       sxgro2,sygro2,szgro2,sgtgro,read_step,write_step,write_step2,
     +       sgsn,filen,pht,
     +       ifault,checkpoint,isfcvlm,imd5,ivelocity,
     +       mediarestart,io_opt,iperf,idyna,nvar,iost,partdeg,SoCalQ)
        if (rank==master)  write(*,*) 'in main, after rddt3d: pht=', pht
        nxt = nx/npx
        nyt = ny/npy
        nzt = nz/npz
        ndt = nd

c     ---------------
c     define timestep
c     ---------------
        nt = nint(tmax/dt) + 1
c     print *, 'tmax, dt, nt', tmax, dt, nt
        print *,'pmcl3d mediarestart',ifault, mediarestart
      end if

C     ==========================
C     BROADCAST A FEW PARAMETERS
C     ==========================

      call MPI_BCAST(ifault,1,MPI_INTEGER, master,MCW,err)      
      call MPI_BCAST(idyna,1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(SoCalQ,1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nvar,1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(iost,1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(partdeg,1,MPI_INTEGER, master,MCW,err)
      
C     ========================
C     
C     ========================

      if(rank==master) then
        call open3d(1,nve,chkp,chkj,insrc,invel,filen,ifault,nxt,nst)
      end if

      call MPI_BARRIER(MCW,err)

C     ======================================
C     BROADCAST PARAMETERS TO ALL PROCESSORS
C     ======================================

c     ----------------------------------------------------------
c     mediarestart = 1 read in media file
c                    2 initialization with writing media file
c                    0 initialization without writing media file
c     ----------------------------------------------------------
       
      call MPI_BCAST(mediarestart,1,MPI_INTEGER, master,MCW,err)

      call MPI_BCAST(checkpoint,1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(io_opt,1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(iperf,1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(isfcvlm,1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(imd5,1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(ivelocity,1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(igreen,1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(invel,filen,MPI_CHARACTER,master,MCW,err)
      call MPI_BCAST(tmax,  1,MPI_REAL,    master,MCW,err)
      call MPI_BCAST(dh,    1,MPI_REAL,    master,MCW,err)
      call MPI_BCAST(dt,    1,MPI_REAL,    master,MCW,err)
      call MPI_BCAST(arbc,  1,MPI_REAL,    master,MCW,err)
      call MPI_BCAST(pht,   1,MPI_REAL,    master,MCW,err)
      call MPI_BCAST(nt,    1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(ndt,   1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nx,    1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(ny,    1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nz,    1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(npx,   1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(npy,   1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(npz,   1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nxt,   1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nyt,   1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nzt,   1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nsrc,  1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nst,   1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nbgx,  1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nedx,  1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nskpx, 1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nbgy,  1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nedy,  1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nskpy, 1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nbgz,  1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nedz,  1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nskpz, 1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(ntiskp,1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nve,   1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(npc,   1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(sxgro,filen,MPI_CHARACTER,master,MCW,err)
      call MPI_BCAST(sygro,filen,MPI_CHARACTER,master,MCW,err)
      call MPI_BCAST(szgro,filen,MPI_CHARACTER,master,MCW,err)
      call MPI_BCAST(nbgx2,  1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nedx2,  1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nskpx2, 1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nbgy2,  1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nedy2,  1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nskpy2, 1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nbgz2,  1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nedz2,  1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(nskpz2, 1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(ntiskp2,1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(ntiskp_sgt,1,MPI_INTEGER, master,MCW,err)
      call MPI_BCAST(sxgro2, filen,MPI_CHARACTER,master,MCW,err)
      call MPI_BCAST(sygro2, filen,MPI_CHARACTER,master,MCW,err)
      call MPI_BCAST(szgro2, filen,MPI_CHARACTER,master,MCW,err)
cSGT -------start---------------------------------------
      call MPI_BCAST(sgt_io_out,1,MPI_INTEGER,master,MCW,err)
      call MPI_BCAST(sgt_numsta,1,MPI_INTEGER,master,MCW,err)
      if (rank .ne. master) allocate(sgt_sta_coord(sgt_numsta,3))
      call MPI_BCAST(sgt_sta_coord,sgt_numsta*3,MPI_INTEGER,master,MCW,err)
      call MPI_BCAST(sgtgro, filen,MPI_CHARACTER,master,MCW,err)
cSGT --------end----------------------------------------
      call MPI_BCAST(read_step, 1,MPI_INTEGER,master,MCW,err)
      call MPI_BCAST(write_step, 1,MPI_INTEGER,master,MCW,err)
      call MPI_BCAST(write_step2, 1,MPI_INTEGER,master,MCW,err)
      call MPI_BCAST(sgsn, filen,MPI_CHARACTER,master,MCW,err)
      call MPI_BCAST(mu_dd, 1,MPI_REAL, master,MCW,err)
      call MPI_BCAST(mu_ss, 1,MPI_REAL, master,MCW,err)
      call MPI_BCAST(fl,    1,MPI_REAL,    master,MCW,err)
      call MPI_BCAST(fh,    1,MPI_REAL,    master,MCW,err)
      call MPI_BCAST(fp,    1,MPI_REAL,    master,MCW,err)
      call MPI_BARRIER(MCW,err)
    
C     ==============================================================
C     DECOMPOSE, ORGANIZE FOR MEMORY DISTRIBUTION AND FIND NEIGHBORS
C     ==============================================================

      dims(1)=npx
      dims(2)=npy
      dims(3)=npz
      call MPI_CART_CREATE(MCW,maxdim,dims,periodic,reorder,MC1,err)      
      call MPI_CART_GET(MC1,maxdim,dims,periodic,coords,err)

c     -------------------
c     MC1 has new rank id 
c     -------------------

      call MPI_CART_RANK(MC1, coords, rankmc1,err)     
      if (rank .eq. 0) transmaster = rankmc1

c     ------------------------------------------------
c     Send this new value of "master" to all processes
c     ------------------------------------------------

      call MPI_BCAST(transmaster,1,MPI_INTEGER,master,MCW,err)
      call MPI_BARRIER(MCW,err)
      master = transmaster 
      rank = rankmc1

      call MPI_CART_SHIFT(MC1,0, 1,fromw,toe,err)
      call MPI_CART_SHIFT(MC1,0,-1,frome,tow,err)
      call MPI_CART_SHIFT(MC1,1, 1,froms,ton,err)
      call MPI_CART_SHIFT(MC1,1,-1,fromn,tos,err)
      call MPI_CART_SHIFT(MC1,2, 1,fromd,tou,err)
      call MPI_CART_SHIFT(MC1,2,-1,fromu,tod,err)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank, 1, 1, 1,swd2neu)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank,-1,-1,-1,neu2swd)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank, 1, 1,-1,swu2ned)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank,-1,-1, 1,ned2swu)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank,-1, 1, 1,sed2nwu)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank, 1,-1,-1,nwu2sed)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank,-1, 1,-1,seu2nwd)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank, 1,-1, 1,nwd2seu)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank, 1, 1, 0,sw2ne)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank,-1,-1, 0,ne2sw)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank,-1, 1, 0,se2nw)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank, 1,-1, 0,nw2se)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank, 0, 1,-1,su2nd)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank, 0,-1, 1,nd2su)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank, 0, 1, 1,sd2nu)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank, 0,-1,-1,nu2sd)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank, 1, 0,-1,wu2ed)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank,-1, 0, 1,ed2wu)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank, 1, 0, 1,wd2eu)
      call CMPI_CART_JUMP(MC1,coords,dims,maxdim,rank,-1, 0,-1,eu2wd)

      call MPI_BARRIER(MC1,err)

C     ===============================================
C     DETERMINE FREE SURFACE DISTRIBUTION AMONG PROCS
C     ===============================================

      if(dims(3)==1) then
        fsf=1
      else if(mod(coords(3)+1,dims(3))==0) then
        fsf=1
      end if

C     ==============================================
C     CREATE FREE SURFACE COMM INCLUDING MASTER PROC
C     ==============================================

      if(fsf==1 .or. rank==master) color=1
      call MPI_COMM_SPLIT(MC1,color,0,MFS,err)

C     ==============================================
C     BEGINING of SOURCE INITIALIZATION
C     ==============================================

      call inisource(rank,size1,MCW,master,invel,dims,coords,maxdim,
     +               npx,npy,npz,nxt,nyt,nzt,nx,ny,nz,
     +               nsrc,read_step,nst,cur_step,srcproc,npsrc,
     +               ix_d,fx_d,y_d,iz_d,fz_d,y_dt,
     +               ifault,idyna)

C     --------END of SOURCE INITIALIZATION-----------

      npts = nxt*nyt*nzt

C     call MPI_BARRIER(MC1,err)
      t4=mpi_wtime()

C     ==============================================
C     RECEIVER PROCS, LOCATION AND INITIALIZATION
C     ==============================================
     
      if (io_opt == 1) then

      allocate(tprec(npts,maxdim))
      allocate(tpmap(npts))
      call getrec(igreen,coords,maxdim,rank,nxt,nyt,nzt,nz,
     +     recproc,nprec,nrec,sgt_recproc,sgt_nprec,
     +     nbgx,nedx,nskpx,nbgy,nedy,nskpy,nbgz,nedz,nskpz,1)
cSGT -------start---------------------------------------
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,sgt_nprec>0,rank,sgt_comm,err)
cSGT --------end----------------------------------------

      if (ifault==3 .or. ifault==4) then
         allocate(tprec(npts,maxdim))
         allocate(tpmap(npts))
         call getrec(igreen,coords,maxdim,rank,nxt,nyt,nzt,nz,
     +        recproc3,nprec3,nrec3,sgt_recproc,sgt_nprec,
     +        1,nx,1,nbgy,nedy,nskpy,1,nz,1,3)
      endif

c     ------------------------------------------------------------------------------
c     no need to call getrec -4, because nbgx...nedz are all the same with getrec -1
c     ------------------------------------------------------------------------------

      nrec_l=nrec                                                 
      nrec2_l=nrec2
      nrec3_l=nrec3
      call MPI_BARRIER(MC1,err)
      t15=mpi_wtime()

C     ===================================================
C     CREATE RECEIVER COMMUNICATOR OPEN FILES FOR WRITING
C     ===================================================

      if(rank==recproc) shape=1
      call MPI_COMM_SPLIT(MC1,shape,0,MCR,err)

C     ============================================
C     PREPARE FOR MPI-IO FILE TYPES     
C
C     nbyte as 4 bytes, nbyte_l as MPI_OFFSET_KIND
C     ============================================

      nbyte = 4
      nbyte_l = nbyte

      allocate(tpmap(nprec*write_step))
      do i = 1, write_step
        do j = 1, nprec
          tpmap((i-1)*nprec+j) = pmap(j)+int((i-1)*nrec,MPI_ADDRESS_KIND)
        end do
      end do
      deallocate(pmap)
      allocate(pmap(nprec*write_step))

c     ----------------------------------------------
c     use the following line on XT3, BG/L, XeonEM64T
c     ----------------------------------------------
      allocate(block(nprec*write_step))

      do i = 1, write_step*nprec
        pmap(i) = tpmap(i)
c     ----------------------------------------------
c     use the following line on XT3, BG/L, XeonEM64T
c     ----------------------------------------------
        block(i) = 1

      end do
      deallocate(tpmap)

c     -------------------------------------
c     use the following line on PWR4, IA-64
c     -------------------------------------
c     call MPI_TYPE_CREATE_INDEXED_BLOCK(nprec*write_step,1,pmap,MPI_REAL,
c
c     ---------------------------------------------------------
c     use the following line on XT3, BG/L, XeonEM64T - 32BIT IO
c     ---------------------------------------------------------
c     call MPI_TYPE_INDEXED(nprec*write_step,block,pmap,MPI_REAL,
c
c     ------------------------------------------------------------
c     use the following 2 lines on XT3, BG/L, XeonEM64T - 64BIT IO
c     ------------------------------------------------------------
      pmap=pmap*4               ! byte addressing instead of element addressing
      call MPI_TYPE_CREATE_HINDEXED(nprec*write_step,block,pmap,MPI_REAL,
c
c     ----------------------------------------------------------
c     each of the 3 options above will use the following line to
c     close off MPI-IO data type creation
c     ----------------------------------------------------------
     +     filetype,err)

      call MPI_TYPE_COMMIT(filetype,err)

c     ----------------------------------------------
c     use the following line on XT3, BG/L, XeonEM64T
c     ----------------------------------------------
      deallocate(block)

c     -----------------------------------------------------------------------
c     set pmap3, filetype3 for both DYN modes; pmap4, filetype4 for SGSN mode
c     pmap3 not needed to host write_step, only last timestep
c     -----------------------------------------------------------------------
      if (ifault==3 .or. ifault==4) then

c     ----------------------------------------------
c     use the following 4 lines on XT3, BG/L, XeonEM64T
c     ----------------------------------------------
         allocate(block(nprec3))
         do i=1, nprec3
            block(i) = 1
         end do

c     -------------------------------------
c     use the following line on PWR4, IA-64
c     -------------------------------------
c     call MPI_TYPE_CREATE_INDEXED_BLOCK(nprec3,1,pmap3,MPI_REAL,
c
c     --------------------------------------------------------
c     use the following line on XT3, BG/L, XeonEM64T -32BIT IO
c     --------------------------------------------------------
c     call MPI_TYPE_INDEXED(nprec3,block,pmap3,MPI_REAL,
c
c     -----------------------------------------------------------
c     use the following 2 lines on XT3, BG/L, XeonEM64T -64BIT IO
c     -----------------------------------------------------------
      pmap3=pmap3*4       ! byte addressing instead of element addressing
      call MPI_TYPE_CREATE_HINDEXED(nprec3,block,pmap3,MPI_REAL,
c
c     ----------------------------------------------------------
c     each of the 3 options above will use the following line to
c     close off MPI-IO data type creation
c     ----------------------------------------------------------
     +     filetype3,err)

      call MPI_TYPE_COMMIT(filetype3,err)
      deallocate(pmap3)

c     ----------------------------------------------
c     use the following line on XT3, BG/L, XeonEM64T
c     ----------------------------------------------
      deallocate(block)

      end if

c     -----------------------------------------------
c     output 19 more variables for SGSN+ifault=4 mode
c     -----------------------------------------------
      if (idyna==1 .and. ifault==4) then
         allocate(tpmap(nprec*write_step*19))
         do i = 1, write_step
            do k=1,19
               do j = 1, nprec
                  tpmap((i-1)*19*nprec+(k-1)*nprec+j) = pmap(j)+(k-1)*nrec+(i-1)*19*nrec
               end do
            end do
         end do

         allocate(pmap4(nprec*write_step*19))
c     -----------------------------------------------
c     use the following line on XT3, BG/L, XeonEM64T
c     -----------------------------------------------
         allocate(block(nprec*write_step*19))

         do i = 1, write_step*nprec*19
            pmap4(i) = tpmap(i)

c     ----------------------------------------------
c     use the following line on XT3, BG/L, XeonEM64T
c     ----------------------------------------------
          block(i) = 1

       end do
       deallocate(tpmap)

c     -------------------------------------
c     use the following line on PWR4, IA-64
c     -------------------------------------
c     call MPI_TYPE_CREATE_INDEXED_BLOCK(nprec*write_step*19,1,pmap4,MPI_REAL,
c
c     ---------------------------------------------------------
c     use the following line on XT3, BG/L, XeonEM64T - 32BIT IO
c     ---------------------------------------------------------
c     call MPI_TYPE_INDEXED(nprec*write_step*19,block,pmap4,MPI_REAL,
c
c     ------------------------------------------------------------
c     use the following 2 lines on XT3, BG/L, XeonEM64T - 64BIT IO
c     ------------------------------------------------------------
       pmap4=pmap4*4            ! byte addressing instead of element addressing
       call MPI_TYPE_CREATE_HINDEXED(nprec*write_step*19,block,pmap4,MPI_REAL,
c
c     ----------------------------------------------------------
c     each of the 3 options above will use the following line to
c     close off MPI-IO data type creation
c     ----------------------------------------------------------
     +      filetype4,err)

       call MPI_TYPE_COMMIT(filetype4,err)
       deallocate(pmap4)

c     ----------------------------------------------
c     use the following line on XT3, BG/L, XeonEM64T
c     ----------------------------------------------
       deallocate(block)

      end if
      end if  ! (io_opt == 1)
C     -----------END OF MPI-IO FILE TYPE CREATION---------------

C     ========================================================
C     MEMORY AND BOUNDS FOR WAVEFIELD VARIABLES, MEDIA AND ABC
C     ========================================================

c     define mpi_subarry begins

      nmtype(1) = nxt+4
      nmtype(2) = nyt+4
      nmtype(3) = nzt+4

      ntype(1)  = nxt+4
      ntype(2)  = nyt+4
      ntype(3)  = 1
      send_offset(1) = 0
      send_offset(2) = 0
      send_offset(3) = 2
      recv_offset(1) = 0
      recv_offset(2) = 0
      recv_offset(3) = nzt+2
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, send_offset, mpi_order_fortran, MPI_REAL, send_u2d_1, err )
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, recv_offset, mpi_order_fortran, MPI_REAL, recv_u2d_1, err )
      call MPI_TYPE_COMMIT( send_u2d_1, err )
      call MPI_TYPE_COMMIT( recv_u2d_1, err )

      ntype(1)  = nxt+4
      ntype(2)  = nyt+4
      ntype(3)  = 2
      send_offset(1) = 0
      send_offset(2) = 0
      send_offset(3) = 2
      recv_offset(1) = 0
      recv_offset(2) = 0
      recv_offset(3) = nzt+2
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, send_offset, mpi_order_fortran, MPI_REAL, send_u2d_2, err )
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, recv_offset, mpi_order_fortran, MPI_REAL, recv_u2d_2, err )
      call MPI_TYPE_COMMIT( send_u2d_2, err )
      call MPI_TYPE_COMMIT( recv_u2d_2, err )

      ntype(1)  = nxt+4
      ntype(2)  = nyt+4
      ntype(3)  = 1
      send_offset(1) = 0
      send_offset(2) = 0
      send_offset(3) = nzt+1
      recv_offset(1) = 0
      recv_offset(2) = 0
      recv_offset(3) = 1
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, send_offset, mpi_order_fortran, MPI_REAL, send_d2u_1, err )
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, recv_offset, mpi_order_fortran, MPI_REAL, recv_d2u_1, err )
      call MPI_TYPE_COMMIT( send_d2u_1, err )
      call MPI_TYPE_COMMIT( recv_d2u_1, err )

      ntype(1)  = nxt+4
      ntype(2)  = nyt+4
      ntype(3)  = 2
      send_offset(1) = 0
      send_offset(2) = 0
      send_offset(3) = nzt
      recv_offset(1) = 0
      recv_offset(2) = 0
      recv_offset(3) = 0
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, send_offset, mpi_order_fortran, MPI_REAL, send_d2u_2, err )
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, recv_offset, mpi_order_fortran, MPI_REAL, recv_d2u_2, err )
      call MPI_TYPE_COMMIT( send_d2u_2, err )
      call MPI_TYPE_COMMIT( recv_d2u_2, err )

      ntype(1)  = nxt+4
      ntype(2)  = 1
      ntype(3)  = nzt+4
      send_offset(1) = 0
      send_offset(2) = nyt+1
      send_offset(3) = 0
      recv_offset(1) = 0
      recv_offset(2) = 1
      recv_offset(3) = 0
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, send_offset, mpi_order_fortran, MPI_REAL, send_s2n_1, err )
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, recv_offset, mpi_order_fortran, MPI_REAL, recv_s2n_1, err )
      call MPI_TYPE_COMMIT( send_s2n_1, err )
      call MPI_TYPE_COMMIT( recv_s2n_1, err )

      ntype(1)  = nxt+4
      ntype(2)  = 2
      ntype(3)  = nzt+4
      send_offset(1) = 0
      send_offset(2) = nyt
      send_offset(3) = 0
      recv_offset(1) = 0
      recv_offset(2) = 0
      recv_offset(3) = 0
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, send_offset, mpi_order_fortran, MPI_REAL, send_s2n_2, err )
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, recv_offset, mpi_order_fortran, MPI_REAL, recv_s2n_2, err )
      call MPI_TYPE_COMMIT( send_s2n_2, err )
      call MPI_TYPE_COMMIT( recv_s2n_2, err )

      ntype(1)  = nxt+4
      ntype(2)  = 1
      ntype(3)  = nzt+4
      send_offset(1) = 0
      send_offset(2) = 2
      send_offset(3) = 0
      recv_offset(1) = 0
      recv_offset(2) = nyt+2
      recv_offset(3) = 0
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, send_offset, mpi_order_fortran, MPI_REAL, send_n2s_1, err )
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, recv_offset, mpi_order_fortran, MPI_REAL, recv_n2s_1, err )
      call MPI_TYPE_COMMIT( send_n2s_1, err )
      call MPI_TYPE_COMMIT( recv_n2s_1, err )

      ntype(1)  = nxt+4
      ntype(2)  = 2
      ntype(3)  = nzt+4
      send_offset(1) = 0
      send_offset(2) = 2
      send_offset(3) = 0
      recv_offset(1) = 0
      recv_offset(2) = nyt+2
      recv_offset(3) = 0
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, send_offset, mpi_order_fortran, MPI_REAL, send_n2s_2, err )
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, recv_offset, mpi_order_fortran, MPI_REAL, recv_n2s_2, err )
      call MPI_TYPE_COMMIT( send_n2s_2, err )
      call MPI_TYPE_COMMIT( recv_n2s_2, err )

      ntype(1)  = 1
      ntype(2)  = nyt+4
      ntype(3)  = nzt+4
      send_offset(1) = nxt+1
      send_offset(2) = 0
      send_offset(3) = 0
      recv_offset(1) = 1
      recv_offset(2) = 0
      recv_offset(3) = 0
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, send_offset, mpi_order_fortran, MPI_REAL, send_w2e_1, err )
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, recv_offset, mpi_order_fortran, MPI_REAL, recv_w2e_1, err )
      call MPI_TYPE_COMMIT( send_w2e_1, err )
      call MPI_TYPE_COMMIT( recv_w2e_1, err )

      ntype(1)  = 2
      ntype(2)  = nyt+4
      ntype(3)  = nzt+4
      send_offset(1) = nxt
      send_offset(2) = 0
      send_offset(3) = 0
      recv_offset(1) = 0
      recv_offset(2) = 0
      recv_offset(3) = 0
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, send_offset, mpi_order_fortran, MPI_REAL, send_w2e_2, err )
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, recv_offset, mpi_order_fortran, MPI_REAL, recv_w2e_2, err )
      call MPI_TYPE_COMMIT( send_w2e_2, err )
      call MPI_TYPE_COMMIT( recv_w2e_2, err )

      ntype(1)  = 1
      ntype(2)  = nyt+4
      ntype(3)  = nzt+4
      send_offset(1) = 2
      send_offset(2) = 0
      send_offset(3) = 0
      recv_offset(1) = nxt+2
      recv_offset(2) = 0
      recv_offset(3) = 0
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, send_offset, mpi_order_fortran, MPI_REAL, send_e2w_1, err )
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, recv_offset, mpi_order_fortran, MPI_REAL, recv_e2w_1, err )
      call MPI_TYPE_COMMIT( send_e2w_1, err )
      call MPI_TYPE_COMMIT( recv_e2w_1, err )

      ntype(1)  = 2
      ntype(2)  = nyt+4
      ntype(3)  = nzt+4
      send_offset(1) = 2
      send_offset(2) = 0
      send_offset(3) = 0
      recv_offset(1) = nxt+2
      recv_offset(2) = 0
      recv_offset(3) = 0
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, send_offset, mpi_order_fortran, MPI_REAL, send_e2w_2, err )
      call MPI_TYPE_CREATE_SUBARRAY( 3, nmtype, ntype, recv_offset, mpi_order_fortran, MPI_REAL, recv_e2w_2, err )
      call MPI_TYPE_COMMIT( send_e2w_2, err )
      call MPI_TYPE_COMMIT( recv_e2w_2, err )

c     define mpi_subarry ends

      call vsemem(nve,nxt,nyt,nzt) 

C     ========================================================
C     DETERMINE PROC BOUNDS FOR EACH COMPUTATIONAL ZONE IF PML
C     ========================================================

      if(npc==1) then
        call pmlmem(dims,coords,maxdim,nxt,nyt,nzt,ndt,fs,fw,fn,fe,fb) 
        call bndmem(maxzone)
        call bounds(nxt,nyt,nzt,ndt,fs,fw,fn,fe,fb)
      end if

C     ===============================
C     BEGINNIG OF MESH INITIALIZATION 
C     ===============================

      if (iperf == 1) then
        call MPI_BARRIER(MC1,err)
        t2=mpi_wtime()
      end if 

      if (rank==master) print *,'before inimesh'
      call inimesh(rank,size1,MC1,master,invel,dims,coords,maxdim,
     +            npx,npy,npz,nxt,nyt,nzt,nx,ny,nz,vpe,vse,dde,
     +            taumin,taumax,fl,fh,fp,bsize,y_dt,idyna,nve,mediarestart,nvar,iost,partdeg,SoCalQ)

      if (rank==master) print *,'after inimesh' 
      if (iperf == 1) then
        call MPI_BARRIER(MC1,err)
        t3=mpi_wtime()
      end if 

C     -------END OF MESH INITIALIZATION-----------

C     =====================
C     SWAP MESH INFORMATION
C     =====================

      call mediaswp(rank,nve,MC1,nxt,nyt,nzt,ton,froms,tos,fromn,
     +              toe,fromw,tow,frome,tou,fromd,tod,fromu,
     +              swd2neu,neu2swd,swu2ned,ned2swu,
     +              sed2nwu,nwu2sed,seu2nwd,nwd2seu,
     +              sw2ne,ne2sw,se2nw,nw2se,su2nd,nd2su,
     +              sd2nu,nu2sd,wu2ed,ed2wu,wd2eu,eu2wd)

C     ==============================================================
C     PML OR CERJAN DAMPING AND TAU INITIALIZATION IF VISCO (VSE==1)
C     ==============================================================

      if(npc==1) then
        call inipml(dh,vse,ndt,arbc)
      else
        call inicrj(coords,maxdim,nxt,nyt,nzt,ndt,nx,ny,arbc)
      end if

      if(nve==1) call tausub(taumin,taumax) 
 
C     ===================================================
C     CHECK WRITE OUT OF PARAMETERS AND SEISMOGRAM HEADER
C     ===================================================

      if(rank==master) then

        call wrpm3d(igreen,nxt,nyt,nzt,nt,dh,dt,npc,
     +       arbc,ntiskp,nbgx,nedx,nskpx,nbgy,nedy,nskpy,
     +       nbgz,nedz,nskpz,vpe,vse,dde,nve,fl,fh,fp)

      end if

C     ======================
C     FOR SGSN DYNAMIC MODEL
C     ======================

c     ------------------------------------------------------------
c     estimation of  eta factor of viscous damping for SGSN dyna
c     eta=0.3 is the eta base value for STABILITY CRITERIA = 0.396
c     ------------------------------------------------------------

      eta=eta*0.396*dh/(vpe(2)*dt)

      call MPI_BARRIER(MC1,err)

C     =======================================
C     INITIATE POINT-SOURCE AND SWAP STRESSES
C     OR CALCULATE DYNAMIC SOURCE
C     =======================================

      i=1

      if(rank==srcproc) then
         if(ifault==3 .or. ifault==4) then
            if(idyna == 0) then              ! stress-glut method
               call dyna(i,mu_dd,mu_ss,
     +              dh,dt,npsrc,nxt,nyt,nzt,
     +              nsrc,nst,rank)
            else                             ! sgsn method
               call srcprmswp(rank,MC1,nxt,nzt,fx_d,fz_d,ix_d,iz_d,
     +              toe,fromw,tow,frome,tou,fromd,tod,fromu,
     +              wu2ed,ed2wu,wd2eu,eu2wd)
               ! ---------------------------------------------------------
               ! determine if dynamic fault reach free-surface among procs
               ! ---------------------------------------------------------
               if(fsf.eq.1.and.nzt.eq.fz_d) fsr=1
               if(fsr==1) call fsrcprm(ix_d,fx_d,nzt)
          end if
        else     ! match 'if (ifault==3|4)'
           ! -------------------------
           ! add another point source
           ! ------------------------- 
           call addsrc(i,dh,dt,nst,npsrc,read_step,ifault,rank,igreen,nzt)
        endif
      endif

      call MPI_BARRIER(MC1,err)
 
      call strswp(rank,MC1,nxt,nyt,nzt,ton,froms,tos,fromn,
     +            toe,fromw,tow,frome,tou,fromd,tod,fromu)
      if (iperf == 1) then
        call MPI_BARRIER(MC1,err)
        t5=mpi_wtime()
      end if



C     ----------------------------------------------------------
C     BEGIN CHECKPOINTING
C     ----------------------------------------------------------      
      t8=0.0
      t10=0.0
      t12=0.0
      t16=0.0       
      t20=0.0
      ckp_count=0
      saveallcheckpoints = 1

c     print *, 'before write chk'
c     ----------------------------------------
c     If present, restart from checkpoint file
c     ----------------------------------------
      if (checkpoint.ne.0) then
        tmp = 0
        ! --------------------------------------
        ! Maximally allow IOST amount of readers
        ! --------------------------------------
        do iwrite=0,IOST
          if (mod(rank,IOST+1).eq.iwrite) then
            write( ckpfile, '(a,i7.7,a)' ) 'output_ckp/ckp', 
     $      rank, '.hdr'
            open( 9, file=ckpfile, status='old', iostat=err )
            if( err .eq. 0 ) then
              read( 9, * ) tmp 
              close( 9 )
            end if
          end if ! IOST
          call MPI_BARRIER(MC1,err)
        end do 

        call MPI_ALLREDUCE( tmp, t0, 1, MPI_INTEGER,
     $        MPI_MIN, MC1, err )
        if ( t0 .ne. 0 ) then
          if( rank .eq. 0 ) write( 0, * )
     $           'Checkpoint found, starting from step',t0
          ! --------------------------------------
          ! Maximally allow IOST amount of readers
          ! --------------------------------------
          do iwrite=0,IOST
            if (mod(rank,IOST+1).eq.iwrite) then
              write( ckpfile, '(a,i7.7,i7.7,a)') 'output_ckp/ckp', rank, t0, '.bin' 
              open( 9, file=ckpfile, form='unformatted', access='direct',
     $              recl=bsize*(nxt+4)*(nyt+4)*(nzt+4)*3, status='old',
     $              iostat=err)
              read( 9, rec=1 ) u1, v1, w1
              read( 9, rec=2 ) xx, yy, zz
              read( 9, rec=3 ) xy, xz, yz
              read( 9, rec=4 ) r1, r2, r3
              read( 9, rec=5 ) r4, r5, r6
              close( 9 )

              write( ckpfile, '(a,i7.7,i7.7,a)') 'output_ckp/ckpdy', rank, t0, '.bin'
              open( 9, file=ckpfile, form='unformatted', access='direct',
     $              recl=bsize*(nxt+4)*(nyt+4)*(nzt+4)*3, status='old',
     $              iostat=err)
              read( 9, rec=1 ) u2,broken,mu
              close( 9 )

              if (npc==1) then
              write( ckpfile, '(a,i7.7,i7.7,a)') 'output_ckp/ckpew', rank, t0, '.bin'
              open( 9, file=ckpfile, form='unformatted', access='direct',
     $                recl=bsize*nzt*nyt*ndt*6, status='old', iostat=err)
                read( 9, rec=1 ) ueT, ueH, veT, veH, weT, weH
                read( 9, rec=2 ) xxeT, xxeH, yyeT, yyeH, zzeT, zzeH
                read( 9, rec=3 ) xzeT, xzeH, yzeT, yzeH, xyeT, xyeH
                read( 9, rec=4 ) uwT, uwH, vwT, vwH, wwT, wwH
                read( 9, rec=5 ) xxwT, xxwH, yywT, yywH, zzwT, zzwH
                read( 9, rec=6 ) xzwT, xzwH, yzwT, yzwH, xywT, xywH
                close( 9 )
                write( ckpfile, '(a,i7.7,i7.7,a)') 'output_ckp/ckpns', rank, t0, '.bin'
              open( 9, file=ckpfile, form='unformatted', access='direct',
     $              recl=bsize*nxt*nzt*ndt*6, status='old', iostat=err) 
                read( 9, rec=1 ) unT, unH, vnT, vnH, wnT, wnH
                read( 9, rec=2 ) xxnT, xxnH, yynT, yynH, zznT, zznH
                read( 9, rec=3 ) xznT, xznH, yznT, yznH, xynT, xynH
                read( 9, rec=4 ) usT, usH, vsT, vsH, wsT, wsH
                read( 9, rec=5 ) xxsT, xxsH, yysT, yysH, zzsT, zzsH
                read( 9, rec=6 ) xzsT, xzsH, yzsT, yzsH, xysT, xysH
                close( 9 )
                write( ckpfile, '(a,i7.7,i7.7,a)') 'output_ckp/ckpb', rank, t0, '.bin'
              open( 9, file=ckpfile, form='unformatted', access='direct',
     $              recl=bsize*nxt*nyt*ndt*6, status='old', iostat=err)
                read( 9, rec=1 ) ubT, ubH, vbT, vbH, wbT, wbH
                read( 9, rec=2 ) xxbT, xxbH, yybT, yybH, zzbT, zzbH
                read( 9, rec=3 ) xzbT, xzbH, yzbT, yzbH, xybT, xybH
                close( 9 )
              end if ! NPC
            end if ! IOST
            call MPI_BARRIER(MC1,err)
          end do ! iwrite
        end if ! t0
      endif ! checkpoint
C     -------------------END OF CHECKPOINT SETUP-----------------------

C     ----------------------------------
C     OPEN FILES IF OUTPUT IS AGGREGATED
C     ----------------------------------
      t0 = t0 + 1
      if((ivelocity .ne. 1) .and. (io_opt .eq. 1)) then         ! if output aggregation....
        write(sxgro_timestep,'(a)') sxgro
        write(sygro_timestep,'(a)') sygro
        write(szgro_timestep,'(a)') szgro
        write(sxgro2_timestep,'(a)') sxgro2
        write(sygro2_timestep,'(a)') sygro2
        write(szgro2_timestep,'(a)') szgro2
        write(sgt_timestep,'(a)') sgtgro
        write(sgsn_timestep,'(a)') sgsn
        write(sxgro_md5,'(a,a)') trim(sxgro),'.md5'
        write(sygro_md5,'(a,a)') trim(sygro),'.md5'
        write(szgro_md5,'(a,a)') trim(szgro),'.md5'
        write(sxgro2_md5,'(a,a)') trim(sxgro2),'.md5'
        write(sygro2_md5,'(a,a)') trim(sygro2),'.md5'
        write(szgro2_md5,'(a,a)') trim(szgro2),'.md5'

        call open3dp(1,MCR,fhx,fhy,fhz,sxgro_timestep,sygro_timestep,
     +            szgro_timestep,fhx_md5,fhy_md5,fhz_md5,sxgro_md5,
     +            sygro_md5,szgro_md5,imd5,ifault,1,filen)

c     -----------------------------------------------
c     if running sgsn mode, more files must be opened
c     -----------------------------------------------
        if (idyna==1 .and. ifault==4) then
          call MPI_FILE_OPEN(MCR,sgsn_timestep,MPI_MODE_WRONLY +
     +                   MPI_MODE_CREATE,MPI_INFO_NULL,fhsgsn,err)

        end if

      end if    ! match 'if (ivelocity .ne. 1)'
      call MPI_BARRIER(MC1,err)

C     =========================
C     GET INITIALIZATION TIMES
C     =========================
      if (iperf==1) then
        dtime(3)=t15-t4
        dtime(4)=t5-t1
        dtime(5)=t3-t2
        dtime(7)=t4-t1
        do idtout=1,11
          dtout(idtout)=0.0d-00
        enddo
        call MPI_ALLREDUCE(dtime,dtout,11,MPI_REAL8,MPI_MAX,MC1,err)
        if (rank==master) then
          write (*,'(" inialization time =",f10.3," sec")') dtout(4)
          write (*,'("    get station time =",f10.3," sec")') dtout(3)
          write (*,'("    read source time =",f10.3," sec")') dtout(7)
          write (*,'("    read media time  =",f10.3," sec")') dtout(5)
        endif
      end if    ! match 'if (iperf==1)'


cSGT -------start---------------------------------------
        if (igreen .ne. -1) then
           allocate(sg1(1:nxt,1:nyt,1:nzt))
           allocate(sg2(1:nxt,1:nyt,1:nzt))
           allocate(sg3(1:nxt,1:nyt,1:nzt))

           do iz = 1,nzt
             do iy = 1,nyt
               do ix = 1,nxt
                  sg1(ix,iy,iz)=(1./lam(ix,iy,iz)+1./mu(ix,iy,iz))/(1./mu(ix,iy,iz)*(3.0/lam(ix,iy,iz)+2.0/mu(ix,iy,iz)))   
                  sg2(ix,iy,iz)=-1./lam(ix,iy,iz)/(2.0/mu(ix,iy,iz)*(3.0/lam(ix,iy,iz)+2.0/mu(ix,iy,iz)))
                  sg3(ix,iy,iz)=mu(ix,iy,iz)
               enddo
             enddo
           enddo
           sgt=0.0
        endif
cSGT --------end--------------------------------------

      if(rank==0) then
        write(*,*) 'npc=',npc,'nve=',nve,'ifault=',ifault
      end if

C     ===============================
C     BEGINNING OF COMPUTATIONAL LOOP
C     ===============================
      do i=t0,nt
         ! ------------------------------
         ! files if output not aggregated
         ! -------------------------------
         if (ivelocity==1) then
            write(sxgro_timestep,'(a,i7.7)') trim(sxgro),i
            write(sygro_timestep,'(a,i7.7)') trim(sygro),i
            write(szgro_timestep,'(a,i7.7)') trim(szgro),i
            write(sxgro2_timestep,'(a,i7.7)') trim(sxgro2),i
            write(sygro2_timestep,'(a,i7.7)') trim(sygro2),i
            write(szgro2_timestep,'(a,i7.7)') trim(szgro2),i
            write(sgt_string_buf,'(i7)') i
            write(sgt_string_buf,'(a)') trim(adjustl(sgt_string_buf))
!******************** non-MPI Output **************************
!            write(sgt_timestep,'(a,a,a)') trim(sgtgro),'-',trim(sgt_string_buf)
!******************** non-MPI Output **************************
!*********************** MPI Output ***************************
            write(sgt_timestep,'(a)') trim(sgtgro)
!*********************** MPI Output ***************************
            write(sgsn_timestep,'(a,i7.7)') trim(sgsn),i
            write(sxgro_md5,'(a,i7.7,a)') trim(sxgro),i,'.md5'
            write(sygro_md5,'(a,i7.7,a)') trim(sygro),i,'.md5'
            write(szgro_md5,'(a,i7.7,a)') trim(szgro),i,'.md5'
            write(sxgro2_md5,'(a,i7.7,a)') trim(sxgro2),i,'.md5'
            write(sygro2_md5,'(a,i7.7,a)') trim(sygro2),i,'.md5'
            write(szgro2_md5,'(a,i7.7,a)') trim(szgro2),i,'.md5'
         end if
         ! ------------------------------------------------
         ! single barrier to synchronize computational loop
         ! ------------------------------------------------
         call MPI_BARRIER(MC1,err)

         if (iperf == 1) then
            call MPI_BARRIER(MC1,err)
         end if

         if(rank==master) write(*,*) 'TIME STEP ', i,' OF ', nt
         t12= t12 + 1.0
         t6 = mpi_wtime()

         call velswp_r(rank,MC1,froms,fromn,fromw,frome,fromd,fromu,
     +        recv_d2u_1,recv_u2d_1,recv_d2u_2,recv_u2d_2,
     +        recv_s2n_1,recv_n2s_1,recv_s2n_2,recv_n2s_2,
     +        recv_w2e_1,recv_e2w_1,recv_w2e_2,recv_e2w_2)

         ! -----------------------------
         ! compute velocity
         ! npc = 1 ==> pml
         ! else    ==> cerjan
         ! -----------------------------

         if(npc==1) then
            call dvel(dh,dt)                 
            call pmlvel(nxt,nyt,nzt,dh,dt,rank)
         else
            call dvelc(dh,dt,nxt,nyt,nzt)
            call addcrjvel(nxt,nyt,nzt)
         end if
        
         ! -------------------------------
         ! compute (sgsn dynamic mode....
         ! 1. velocity on fault plane
         ! 2. swap velocities 
         ! -------------------------------
         if(rank==srcproc) then        
            if(ifault==3.or.ifault==4) then 
               if(idyna == 1) then        
                  call dynavel(ix_d,fx_d,y_d,iz_d,fz_d,nxt,nzt,dh,dt,fsr)
                  call splitvelswp(rank,MC1,nxt,nzt,fx_d,y_d,fz_d,ix_d,iz_d,
     +                 toe,fromw,tow,frome,tou,fromd,tod,fromu,
     +                 wu2ed,ed2wu,wd2eu,eu2wd)
               end if
            endif
         endif         

         ! ------------------------------------------------------------
         ! swap internal,pml regional, and free surface velocity planes
         ! ------------------------------------------------------------
         call velswp_s(rank,MC1,ton,tos,toe,tow,tou,tod,
     +        send_d2u_1,send_u2d_1,send_d2u_2,send_u2d_2,
     +        send_s2n_1,send_n2s_1,send_s2n_2,send_n2s_2,
     +        send_w2e_1,send_e2w_1,send_w2e_2,send_e2w_2)

         call wait_all_vel()         
         ! ------------------------------------------------
         ! compute slip and slip rate for sgsn dynamic mode
         ! ------------------------------------------------
         if(rank==srcproc)  then
            if(ifault==3.or.ifault==4) then
               if(idyna == 1)
     +              call sourceslip(ix_d,fx_d,y_d,iz_d,fz_d,nzt,fsr,dh,dt)
            end if
         endif
 
         ! ---------------------------------
         ! compute (x,y) velocity above
         ! free surface
         ! ---------------------------------
         if(fsf==1) call fvelxy(nxt,nyt,nzt)

         ! ---------------------------------
         ! compute (x,y) split velocities
         ! above free surface 
         ! ---------------------------------
         if(fsr==1) call fvelxysplit(ix_d,fx_d,nzt)        
         
         ! -----------------------------
         ! swap free surface velocities
         ! -----------------------------
         if(fsf==1) call fvelswpxy(rank,MC1,nxt,nyt,nzt,
     +        ton,froms,tow,frome)
         
         ! -------------------------------
         ! compute z velocity above free
         ! surface, swapping not necessary
         ! -------------------------------
         !     if(fsf==1) call fvelzma(nxt,nyt,nzt)

         if(fsf==1) call fvelz(nxt,nyt,nzt)


         ! -----------------------------------
         ! compute split z velocity above free
         ! surface, swapping not necessary
         ! -----------------------------------
         if(fsr==1) call fvelzsplit(ix_d,fx_d,y_d,nzt)

         ! -----------------------------------------------
         ! compute stress
         !
         ! nve==1 and npc==1 ==> visco material and pml
         ! nve==1 and npc==0 ==> visco material and cerjan
         ! nve==0 and npc==1 ==> elas. material and pml
         ! nve==0 and npc==0 ==> elas. material and cerjan
         ! -----------------------------------------------


         call strswp_r(rank,MC1,froms,fromn,fromw,frome,fromd,fromu,
     +        recv_d2u_1,recv_u2d_1,recv_d2u_2,recv_u2d_2,
     +        recv_s2n_1,recv_n2s_1,recv_s2n_2,recv_n2s_2,
     +        recv_w2e_1,recv_e2w_1,recv_w2e_2,recv_e2w_2)

         if(nve==1 .and. npc==1) then
            call dstrq(coords,nxt,nyt,nzt,dh,dt,nz)
            call pmlstr(nxt,nyt,nzt,dh,dt,rank)
         else if(nve==1 .and. npc==0) then
            call dstrqc(coords,nxt,nyt,nzt,dh,dt,nz,rank) 
            call addcrjstr(nxt,nyt,nzt)
         else if(nve==0 .and. npc==1) then 
            call dstr(dh,dt)
            call pmlstr(nxt,nyt,nzt,dh,dt,rank)
         else
            call dstrc(dh,dt,nxt,nyt,nzt)
            call addcrjstr(nxt,nyt,nzt)
         end if

         ! ---------------------------------
         ! compute stress above free surface
         ! ---------------------------------
c     if(fsf==1) call fstr(nxt,nyt,nzt)

         if(rank==srcproc) then
            if(ifault==3 .or. ifault==4) then
               ! -------------------------------------------
               ! compute (sg and sgsn dynamic methods)......
               ! 1. dynamic stress on fault
               ! 2. stress above free surface
               ! 3. swap stress
               ! -------------------------------------------
               if(idyna == 0) then
                  call dyna(i,mu_dd,mu_ss,
     +                 dh,dt,npsrc,nxt,nyt,nzt,
     +                 nsrc,nst,rank)
               else
                  call dynastress(i,ix_d,fx_d,y_d,iz_d,fz_d,nxt,nzt,
     +                 fsr,dh,dt)   
                  if(fsr==1) call fstrsplit(ix_d,fx_d,y_d,nzt)          
                  call splitstrswp(rank,MC1,nxt,nzt,fx_d,y_d,fz_d,ix_d,iz_d,
     +                 toe,fromw,tow,frome,tou,fromd,tod,fromu)   
               end if           

               if(fsf==1) call fstr(nxt,nyt,nzt)
            else                
               ! -------------------------------------------
               ! 1. determine if processor is source bearing
               ! 2. read subgrid fault files if processor
               !    is source bearing.
               ! 3. add source if necessary 
               ! -------------------------------------------
               if (ifault==2 .and. nst.gt.read_step .and. mod(i,read_step)==0 .and. i.lt.nst)  then
                  call read_src_ifault_2(coords,maxdim,MC1,rank,master,nxt,nyt,nzt,srcproc,
     +                 nsrc,nst,npsrc,nz,i/read_step +1,read_step)
               end if

c  add KBO 7/21/12 for alt surface force

               if ((igreen.ge.4) .and. (igreen.le.6)) then
                  if (fsf==1) call fstr(nxt,nyt,nzt)
                  call addsrc(i+1,dh,dt,nst,npsrc,read_step,ifault,rank,igreen,nzt)
               else
                  call addsrc(i+1,dh,dt,nst,npsrc,read_step,ifault,rank,igreen,nzt)
                  if (fsf==1) call fstr(nxt,nyt,nzt)
               end if
            endif               !end of (ifault==)
         else
            if (fsf==1) call fstr(nxt,nyt,nzt)
         endif                  !end of (rank=srcproc)

c  add KBO 7/21/12 for alt surface force


c              call addsrc(i+1,dh,dt,nst,npsrc,read_step,ifault,rank,igreen)
c           endif               !end of (ifault==)
c        endif                  !end of (rank=srcproc)
c        
c        if(fsf==1) call fstr(nxt,nyt,nzt)        
         
         ! ---------------------
         ! swap planes of stress
         ! ---------------------
         call strswp_s(rank,MC1,ton,tos,toe,tow,tou,tod,
     +        send_d2u_1,send_u2d_1,send_d2u_2,send_u2d_2,
     +        send_s2n_1,send_n2s_1,send_s2n_2,send_n2s_2,
     +        send_w2e_1,send_e2w_1,send_w2e_2,send_e2w_2)

         call wait_all_vel()

         ! -----------------------------------------
         ! compute restoring forces/traction on fault
         ! for sgsn dynamic mode
         ! 1. force/traction computation on fault
         ! 2. traction swap among neighbors
         ! 3. friction computation on fault
         ! -----------------------------------------
         if(rank==srcproc) then 
            if(ifault==3.or.ifault==4) then 
               if(idyna == 1) then 
                  call dynarestrac(i,ix_d,fx_d,y_d,iz_d,fz_d,nxt,nzt,
     +                 fsr,dh,dt,eta)        
                  call tracswp(rank,MC1,nxt,nzt,fx_d,y_d,fz_d,ix_d,iz_d,
     +                 toe,fromw,tow,frome,tou,fromd,tod,fromu,
     +                 wu2ed,ed2wu,wd2eu,eu2wd)     
                  call sourcestres(i,ix_d,fx_d,y_d,iz_d,fz_d,nzt,fsr,dh,dt)                    
               end if
            end if
         end if       

cSGT -------start---------------------------------------
!        if (igreen .ne. -1 .and. mod(i,ntiskp_sgt).eq.0) then
! off-diagonal stress on staggered-grid interpolate to diagonal stress (P. Chen ?)
c           xyn=(xy(0:nxt-1,1:nyt,1:nzt)+xy(1:nxt,1:nyt,1:nzt)+xy(0:nxt-1,0:nyt-1,1:nzt)+xy(1:nxt,0:nyt-1,1:nzt))/4.0
c           xzn=(xz(0:nxt-1,1:nyt,1:nzt)+xz(1:nxt,1:nyt,1:nzt)+xz(0:nxt-1,1:nyt,0:nzt-1)+xz(1:nxt,1:nyt,0:nzt-1))/4.0
c           yzn=(yz(1:nxt,0:nyt-1,1:nzt)+yz(1:nxt,1:nyt,1:nzt)+yz(1:nxt,0:nyt-1,0:nzt-1)+yz(1:nxt,1:nyt,0:nzt-1))/4.0

c corrected averaging of shear stresses to normal stresses, KBO 7/13/12
c           xyn=(xy(2:nxt+1,1:nyt,1:nzt)+xy(1:nxt,1:nyt,1:nzt)+xy(2:nxt+1,0:nyt-1,1:nzt)+xy(1:nxt,0:nyt-1,1:nzt))/4.0
c           xzn=(xz(2:nxt+1,1:nyt,1:nzt)+xz(1:nxt,1:nyt,1:nzt)+xz(2:nxt+1,1:nyt,0:nzt-1)+xz(1:nxt,1:nyt,0:nzt-1))/4.0
c           yzn=(yz(1:nxt,0:nyt-1,1:nzt)+yz(1:nxt,1:nyt,1:nzt)+yz(1:nxt,0:nyt-1,0:nzt-1)+yz(1:nxt,1:nyt,0:nzt-1))/4.0
!
!           IMPLEMENTED IN THE OUTPUT WRITING PART BELOW!
!
!           xyn=xy(1:nxt,1:nyt,1:nzt)
!           xzn=xz(1:nxt,1:nyt,1:nzt)
!           yzn=yz(1:nxt,1:nyt,1:nzt)
!
!           do j = 1,sgt_nprec
!               ix = sgt_prec(j,1)
!               iy = sgt_prec(j,2)
!               iz = sgt_prec(j,3)
!               sgt(1,j)=sg1(ix,iy,iz)*xx(ix,iy,iz)+sg2(ix,iy,iz)*yy(ix,iy,iz)+sg2(ix,iy,iz)*zz(ix,iy,iz)
!               sgt(2,j)=sg2(ix,iy,iz)*xx(ix,iy,iz)+sg1(ix,iy,iz)*yy(ix,iy,iz)+sg2(ix,iy,iz)*zz(ix,iy,iz)
!               sgt(3,j)=sg2(ix,iy,iz)*xx(ix,iy,iz)+sg2(ix,iy,iz)*yy(ix,iy,iz)+sg1(ix,iy,iz)*zz(ix,iy,iz)
! off-diagonal stress on staggered-grid interpolate to diagonal stress
!               sgt(4,j)=sg3(ix,iy,iz)*xyn(ix,iy,iz)
!               sgt(5,j)=sg3(ix,iy,iz)*xzn(ix,iy,iz)
!               sgt(6,j)=sg3(ix,iy,iz)*yzn(ix,iy,iz)
!           enddo
!        endif
cSGT --------end--------------------------------------
         
         ! ------------------------------------------
         ! write seismograms for dynamic mode
         ! ------------------------------------------
         if (iperf == 1) then
            call MPI_BARRIER(MC1,err)
            t7 = mpi_wtime()
            t8 = t8+ (t7 - t6)
         end if
         
         if(rank==master) then
            if (ifault==3 .or. ifault==4) then
               write(24,'(I6,3E10.3)') i,u2(ndt+2,ndt+2,ndt+2),
     +              rupt(ndt+2,ndt+2,ndt+2),u3(ndt+2,ndt+2,ndt+2)
               
            else
               write(24,'(I6,3E10.3)') i, u1(ndt+2,ndt+2,ndt+2),
     +              v1(ndt+2,ndt+2,ndt+2),w1(ndt+2,ndt+2,ndt+2)

c     ----------------------------------------------------
c     use the following line on platforms other than BG/L 
c     ----------------------------------------------------
               call flush(24)

c     ----------------------------------------------------
c     use the following line on BG/L
c     ----------------------------------------------------
c     call flush_(24)

            end if
         endif

         ! -----------------------------------------------
         ! move wrtrec here to save memory use
         ! move file_set_view, rather moving indexed_block
         ! call write_at_all rather write_all
         ! use only one bufxyz, rather 3 bufx/bufy/bufz
         ! define disp 8-bytes
         ! -----------------------------------------------
c     call wrtrec(i,ntiskp,fhx,fhy,fhz,nprec,nrec_l)
c     wrtrec added in here
         
         nsrc_l=nsrc
         
!     -----------------------------------------
!     Begin writing to files if io_opt = 1
!
!     The following cases are considered
!     1. wave propagation mode (ifault = 0,1,2)
!        * isfcvlm = 0 ==> 
!          only surface velocity output enabled
!        * isfcvlm = 1 ==> 
!          volume and surface velocity output enabled
!     2. dynamic rupture mode (ifault = 3,4)
!        * ifault  = 3 and idyna = 0 ==>
!          surface velocity output disbaled but
!          sg dynamic output enabled
!        * ifault  = 4 and idyna = 1 ==>
!          surface velocity output enabled and sgsn 
!          dynamic output enabled
!        * ifault  = 4 and idyna = 0 ==>
!          surface velocity output enabled and sg
!          dynamic output enabled
!
!     The file writing procedure consists of following steps
!     1. open necessary files
!     2. set position in opened files to commence writing
!     3. collect data in buffers
!     4. write buffered data to files 
!     ------------------------------------------------------

         if (io_opt == 1) then
cSGT -------start---------------------------------------
            !SGT's extended source output mode
            if ((sgt_io_out .eq. 1) .and. (igreen .ne. -1)) then
                it = ntiskp_sgt*int(i/ntiskp_sgt)
                if ((i==it) .and. (sgt_nprec .gt. 0)) then
                    ! buffer allocation
                    if (sgt_count == 0) then
                        allocate(sgt_buf(sgt_nprec,6,write_step))
                        sgt_buf=-1.0
                    endif
                    sgt_count=sgt_count+1
                    ! store data to buffer
                    do j=1,sgt_nprec
                        ix = sgt_prec(j,1)
                        iy = sgt_prec(j,2)
                        iz = sgt_prec(j,3)
                        sgt_buf(j,1,sgt_count)=sg1(ix,iy,iz)*xx(ix,iy,iz)+sg2(ix,iy,iz)*yy(ix,iy,iz)+sg2(ix,iy,iz)*zz(ix,iy,iz)
                        sgt_buf(j,2,sgt_count)=sg2(ix,iy,iz)*xx(ix,iy,iz)+sg1(ix,iy,iz)*yy(ix,iy,iz)+sg2(ix,iy,iz)*zz(ix,iy,iz)
                        sgt_buf(j,3,sgt_count)=sg2(ix,iy,iz)*xx(ix,iy,iz)+sg2(ix,iy,iz)*yy(ix,iy,iz)+sg1(ix,iy,iz)*zz(ix,iy,iz)
! off-diagonal stress on staggered-grid interpolate to diagonal stress
                        sgt_buf(j,4,sgt_count)=sg3(ix,iy,iz)*xy(ix,iy,iz)
                        sgt_buf(j,5,sgt_count)=sg3(ix,iy,iz)*xz(ix,iy,iz)
                        sgt_buf(j,6,sgt_count)=sg3(ix,iy,iz)*yz(ix,iy,iz)
                    end do
                    if (sgt_count == write_step) then
! non-MPI Output doesn't work and needs update.
!******************** non-MPI Output **************************
!                        do j=1,sgt_nprec
!                            write(sgt_string_buf,'(i7)') sgt_local_sta_coord(j,1)
!                            write(sgt_string_buf,'(a)') TRIM(ADJUSTL(sgt_string_buf))
!                            write(sgt_filename,'(a,a,a)') trim(sgt_timestep),'-',
!     +                            trim(sgt_string_buf)
!                            write(sgt_string_buf,'(i7)') sgt_local_sta_coord(j,2)
!                            write(sgt_string_buf,'(a)') TRIM(ADJUSTL(sgt_string_buf))
!                            write(sgt_filename,'(a,a,a)') trim(sgt_filename),'-',
!     +                            trim(sgt_string_buf)
!                            write(sgt_string_buf,'(i7)') sgt_local_sta_coord(j,3)
!                            write(sgt_string_buf,'(a)') TRIM(ADJUSTL(sgt_string_buf))
!                            write(sgt_filename,'(a,a,a)') trim(sgt_filename),'-',
!     +                            trim(sgt_string_buf)
!                            open(20,file=sgt_filename,form='unformatted',
!     +                           access='direct',recl=floatsize*6*write_step)
!                            write(20,rec=1) (sgt_buf(j,k),k=1,6*write_step)
!                            close(20)
!                        enddo
!******************** non-MPI Output **************************

!********************   MPI Output   **************************
                        if (sgt_mpi_file .eq. 0) then
			                call MPI_FILE_OPEN(sgt_comm,sgt_timestep,MPI_MODE_WRONLY + MPI_MODE_CREATE,
     +                          MPI_INFO_NULL, sgt_mpi_file,err)
                        endif
                        allocate(sgt_buf_buf(write_step))
			            do j=1,sgt_nprec
			                do k=1,6
                                sgt_offset=1
                                sgt_offset=sgt_offset*sgt_local_sta_coord(j,4)*6*floatsize*
     +                              int(int(nt/write_step)*write_step/ntiskp_sgt)
!                                write (*,*), 'sgt_local_sta_coord=', sgt_local_sta_coord(j,4),'offset=',sgt_offset
                                sgt_offset=sgt_offset+(k-1)*floatsize*
     +                              int(int(nt/write_step)*write_step/ntiskp_sgt)
!                                write (*,*), 'k=', k,'offset=',sgt_offset
                                sgt_offset=sgt_offset+(int(i/ntiskp_sgt/write_step)-1)*floatsize*write_step
!                                write (*,*), 'i/ntiskp_sgt/write_step=', int(i/ntiskp_sgt/write_step),'offset=',sgt_offset
			                    sgt_buf_buf=sgt_buf(j,k,1:write_step)
                                call MPI_FILE_WRITE_AT(sgt_mpi_file,sgt_offset,
     +                              sgt_buf_buf,write_step,
     +                              MPI_REAL,MPI_STATUS_IGNORE,err)
                            end do
			            end do
			            deallocate(sgt_buf_buf)
!********************   MPI Output   **************************
			            deallocate(sgt_buf)
                        sgt_count=0
                    endif
                endif
            endif !end (sgt_io_out .eq. 1)
cSGT --------end----------------------------------------
                ! --------------------------------------------------------------------------------
                ! wave propagation mode : output volume velocity every write_step2 timesteps
                !                         if isfcvlm = 1
                ! dynamic mode          : output surface velocity every write_step2 timesteps
                !                         if ifault = 4
                ! --------------------------------------------------------------------------------
                if ((isfcvlm==1 .or. ifault==4) .and. rank==recproc) then
                   it = ntiskp2*int(i/ntiskp2)
                   if (i==it) then
                      if (count2 == 0) then
                         allocate(bufx2(nprec*write_step2))
                         allocate(bufy2(nprec*write_step2))
                         allocate(bufz2(nprec*write_step2))
                         bufx2=-1.0
                         bufy2=-1.0
                         bufz2=-1.0
                      end if

                      count2 = count2 + 1
                      if ((count2 == write_step2) .or. (it==(ntiskp2*int(nt/ntiskp2)))) then
                         if(ivelocity==1) then
                            call open3dp2(1,MCR,fhx2,fhy2,fhz2,sxgro2_timestep,sygro2_timestep,
     +                           szgro2_timestep,fhx2_md5,fhy2_md5,fhz2_md5,sxgro2_md5,
     +                           sygro2_md5,szgro2_md5,imd5,ifault,1,filen)

                            call MPI_FILE_SET_VIEW(fhx2,disp2,MPI_REAL,filetype,'native',
     +                           MPI_INFO_NULL,err)

                            call MPI_FILE_SET_VIEW(fhy2,disp2,MPI_REAL,filetype,'native',
     +                           MPI_INFO_NULL,err)

                            call MPI_FILE_SET_VIEW(fhz2,disp2,MPI_REAL,filetype,'native',
     +                           MPI_INFO_NULL,err)
                         else
                            call MPI_FILE_SET_VIEW(fhx2,disp2,MPI_REAL,filetype,'native',
     +                           MPI_INFO_NULL,err)

                            call MPI_FILE_SET_VIEW(fhy2,disp2,MPI_REAL,filetype,'native',
     +                           MPI_INFO_NULL,err)

                            call MPI_FILE_SET_VIEW(fhz2,disp2,MPI_REAL,filetype,'native',
     +                           MPI_INFO_NULL,err)

                            disp2 = disp2 + write_step2*nrec2_l*nbyte_l
                         end if
                      end if

                      do j=1,nprec
                         xn = prec(j,1)
                         yn = prec(j,2)
                         zn = prec(j,3)
                         bufx2(nprec*(count2-1)+j) = xy(xn,yn,zn)
                         bufy2(nprec*(count2-1)+j) = yz(xn,yn,zn)
                         bufz2(nprec*(count2-1)+j) = xz(xn,yn,zn)
                      end do

                      if ((count2 == write_step2) .or. (it==(ntiskp2*int(nt/ntiskp2)))) then
                         t22=mpi_wtime()
                         call MPI_FILE_WRITE_AT_ALL(fhx2,offset,bufx2,nprec*write_step2,MPI_REAL,status,err)
                         t23=mpi_wtime()
                         call MPI_FILE_WRITE_AT_ALL(fhy2,offset,bufy2,nprec*write_step2,MPI_REAL,status,err)
                         call MPI_FILE_WRITE_AT_ALL(fhz2,offset,bufz2,nprec*write_step2,MPI_REAL,status,err)

                         if(ivelocity==1) then
                            call open3dp2(0,MCR,fhx2,fhy2,fhz2,sxgro2_timestep,sygro2_timestep,szgro2_timestep,
     +                           fhx2_md5,fhy2_md5,fhz2_md5,sxgro2_md5,sygro2_md5,
     +                           szgro2_md5,imd5,ifault,1,filen)
                         end if

                         deallocate(bufx2)
                         deallocate(bufy2)
                         deallocate(bufz2)
                         count2 = 0
                      end if
                   end if
                end if

                ! -------------------------------------------------------------------------------------------------------
                ! wave propagation mode : write surface velocity every write_step time steps
                ! dynamic mode          : write d_u, rate_u every write_step time steps; rupture data on entire fault
                !                         for the last time step
                ! dynamic mode (idyna=1): write additional 19 variables for sgsn ,mode
                ! -------------------------------------------------------------------------------------------------------
                if(rank==recproc) then
                   it = ntiskp*int(i/ntiskp)

                   if(i==it) then
                      if (count == 0) then
                         allocate(bufx(nprec*write_step))
                         allocate(bufy(nprec*write_step))
                         bufx=-1.0
                         bufy=-1.0

                         if (ifault .le. 2) then
                            allocate(bufz(nprec*write_step))
                            bufz=-1.0
                         end if

                         if (idyna==1 .and. ifault==4) then
                            allocate(bufxyz(nprec*write_step*19))
                            bufxyz=-1.0
                         end if
                      end if

                      if ((ifault==3 .or. ifault==4) .and. it==(ntiskp*int(nt/ntiskp))) then
                         allocate(bufz(nprec3))
                         bufz=-1.0
                      end if

                      count = count + 1
                      if ((count == write_step) .or. (it==(ntiskp*int(nt/ntiskp)))) then
                         if (ivelocity==1) then
                            if ((ifault ==3 .or. ifault == 4) .and. it .ne. (ntiskp*int(nt/ntiskp))) then
                               call open3dp(1,MCR,fhx,fhy,fhz,sxgro_timestep,sygro_timestep,
     +                              szgro_timestep,fhx_md5,fhy_md5,fhz_md5,sxgro_md5,
     +                              sygro_md5,szgro_md5,imd5,ifault,0,filen)
                            else
                               call open3dp(1,MCR,fhx,fhy,fhz,sxgro_timestep,sygro_timestep,
     +                              szgro_timestep,fhx_md5,fhy_md5,fhz_md5,sxgro_md5,
     +                              sygro_md5,szgro_md5,imd5,ifault,1,filen)
                            endif


                            call MPI_FILE_SET_VIEW(fhx,disp,MPI_REAL,filetype,'native',
     +                           MPI_INFO_NULL,err)
                            call MPI_FILE_SET_VIEW(fhy,disp,MPI_REAL,filetype,'native',
     +                           MPI_INFO_NULL,err)

                            if (ifault .le. 2) then
                               call MPI_FILE_SET_VIEW(fhz,disp,MPI_REAL,filetype,'native',
     +                              MPI_INFO_NULL,err)
                            else
                               if (it==(ntiskp*int(nt/ntiskp))) then
                                  call MPI_FILE_SET_VIEW(fhz,disp3,MPI_REAL,filetype3,'native',
     +                                 MPI_INFO_NULL,err)
                               end if
                            end if

                            if (idyna==1 .and. ifault==4) then
                               call MPI_FILE_OPEN(MCR,sgsn_timestep,MPI_MODE_WRONLY +
     +                              MPI_MODE_CREATE,MPI_INFO_NULL,fhsgsn,err)

                               call MPI_FILE_SET_VIEW(fhsgsn,disp4,MPI_REAL,filetype4,'native',
     +                              MPI_INFO_NULL,err)
                            end if

                         else

                            call MPI_FILE_SET_VIEW(fhx,disp,MPI_REAL,filetype,'native',
     +                           MPI_INFO_NULL,err)
                            call MPI_FILE_SET_VIEW(fhy,disp,MPI_REAL,filetype,'native',
     +                           MPI_INFO_NULL,err)

                            if (ifault .le. 2) then
                               call MPI_FILE_SET_VIEW(fhz,disp,MPI_REAL,filetype,'native',
     +                              MPI_INFO_NULL,err)
                            else
                               if (it==(ntiskp*int(nt/ntiskp))) then
                                  call MPI_FILE_SET_VIEW(fhz,disp3,MPI_REAL,filetype3,'native',
     +                                 MPI_INFO_NULL,err)
                               end if
                            end if

                            disp = disp + write_step*nrec_l*nbyte_l
                            if (idyna==1 .and. ifault==4) then
                               call MPI_FILE_SET_VIEW(fhsgsn,disp4,MPI_REAL,filetype4,'native',
     +                              MPI_INFO_NULL,err)
                               disp4 = disp4 + write_step*19*nrec_l*nbyte_l
                            end if
                         end if
                      end if

                      ! --------------------------------------------------------
                      ! write u2, u3, and rupture data(last time step only)
                      ! to buffer
                      ! --------------------------------------------------------
                      if (ifault==3 .or. ifault==4) then
                         do j=1,nprec
                            xn = prec(j,1)
                            yn = prec(j,2)
                            zn = prec(j,3)
                            bufx(nprec*(count-1)+j) = u2(xn,yn,zn)
                            bufy(nprec*(count-1)+j) = u3(xn,yn,zn)
                         end do

                         if (it==(ntiskp*int(nt/ntiskp))) then
                            do j=1,nprec3
                               xn = prec3(j,1)
                               yn = prec3(j,2)
                               zn = prec3(j,3)
                               bufz(j) = rupt(xn,yn,zn)
                            end do
                         end if
                      end if

                      ! ---------------------------------------------
                      ! write wave velocity in wave propagation mode
                      ! to buffers
                      ! ---------------------------------------------
                      if (ifault .le. 2) then
                         do j=1,nprec
                            xn = prec(j,1)
                            yn = prec(j,2)
                            zn = prec(j,3)

c changed output to velocities KBO 7/13/12

                            uav = (u1(xn,yn,zn)+u1(xn+1,yn,zn)+
     +                             u1(xn+1,yn+1,zn)+u1(xn,yn+1,zn)+
     +                             u1(xn,yn,zn+1)+u1(xn+1,yn,zn+1)+
     +                             u1(xn+1,yn+1,zn+1)+u1(xn,yn+1,zn+1))/8.

                            vav = (v1(xn,yn,zn)+v1(xn,yn,zn+1))/2.

                            wav = (w1(xn,yn,zn)+w1(xn,yn+1,zn))/2.

c                           bufx(nprec*(count-1)+j) = uav
c                           bufy(nprec*(count-1)+j) = vav
c                           bufz(nprec*(count-1)+j) = wav
                            bufx(nprec*(count-1)+j) = u1(xn,yn,zn)
                            bufy(nprec*(count-1)+j) = v1(xn,yn,zn)
                            bufz(nprec*(count-1)+j) = w1(xn,yn,zn)
c                           bufx(nprec*(count-1)+j) = xx(xn,yn,zn)
c                           bufy(nprec*(count-1)+j) = yy(xn,yn,zn)
c                           bufz(nprec*(count-1)+j) = zz(xn,yn,zn)
                         end do
                      endif

                      ! ----------------------------------------------
                      ! write additional 19 variables in sgsn mode to
                      ! buffers
                      ! ----------------------------------------------
                      if (idyna==1 .and. ifault==4) then
                         do j=1,nprec
                            xn = prec(j,1)
                            zn = prec(j,3)
                            tmp_count2 = (count-1)*19*nprec

                            bufxyz(tmp_count2+j) = rateu_u(xn,zn)
                            bufxyz(tmp_count2+nprec+j)= rateu_w(xn,zn)
                            bufxyz(tmp_count2+2*nprec+j) = slipu_u(xn,zn)
                            bufxyz(tmp_count2+3*nprec+j) = slipu_w(xn,zn)
                            bufxyz(tmp_count2+4*nprec+j) = tru1(xn,zn)
                            bufxyz(tmp_count2+5*nprec+j) = trw1(xn,zn)
                            bufxyz(tmp_count2+6*nprec+j) = trv1(xn,zn)
                            bufxyz(tmp_count2+7*nprec+j) = uplus(xn,zn)
                            bufxyz(tmp_count2+8*nprec+j) = vplus(xn,zn)
                            bufxyz(tmp_count2+9*nprec+j) = wplusi(xn,zn)
                            bufxyz(tmp_count2+10*nprec+j) = umin(xn,zn)
                            bufxyz(tmp_count2+11*nprec+j) = vmin(xn,zn)
                            bufxyz(tmp_count2+12*nprec+j) = wmin(xn,zn)
                            bufxyz(tmp_count2+13*nprec+j) = duplus(xn,zn)
                            bufxyz(tmp_count2+14*nprec+j) = dvplus(xn,zn)
                            bufxyz(tmp_count2+15*nprec+j) = dwplus(xn,zn)
                            bufxyz(tmp_count2+16*nprec+j) = dumin(xn,zn)
                            bufxyz(tmp_count2+17*nprec+j) = dvmin(xn,zn)
                            bufxyz(tmp_count2+18*nprec+j) = dwmin(xn,zn)
                         end do
                      end if

                      ! -------------------------------------------
                      ! Begin writing buffered data to output files
                      ! -------------------------------------------
                      if ((count == write_step) .or. (it==(ntiskp*int(nt/ntiskp)))) then
                         t13=mpi_wtime()
                         call MPI_FILE_WRITE_AT_ALL(fhx,offset,bufx,nprec*write_step,MPI_REAL,status,err)
                         t14=mpi_wtime()
                         call MPI_FILE_WRITE_AT_ALL(fhy,offset,bufy,nprec*write_step,MPI_REAL,status,err)

                         if (idyna==1 .and. ifault==4) then
                            call MPI_FILE_WRITE_AT_ALL(fhsgsn,offset,bufxyz,nprec*write_step*19,MPI_REAL,status,err)
                         end if
                      end if

c     -----------------------------------------------
c     Begin writing md5 files
c     NOTES:
c     1. XT3 Platform : don't use any line referencing
c                       md5
c     2. BG/L Platform : replace md5_sum with md5_sum_
c     ------------------------------------------------
                      t16=t16+(t14-t13)
c     if (imd5==1) then
c     call md5_sum(bufx,nprec,resblock)
c     disp_md5 = rank*32
c     call MPI_FILE_SEEK(fhx_md5,disp_md5,MPI_SEEK_SET,err)
c     call MPI_FILE_WRITE(fhx_md5,resblock,32,MPI_CHARACTER,status,err)

c     call md5_sum(bufy,nprec,resblock)
c     call MPI_FILE_SEEK(fhy_md5,disp_md5,MPI_SEEK_SET,err)
c     call MPI_FILE_WRITE(fhy_md5,resblock,32,MPI_CHARACTER,status,err)
c     end if

                      if ((ifault .ne. 3 .and. ifault .ne. 4) .or. (it==(ntiskp*int(nt/ntiskp)))) then
c     if (imd5==1) then
c     call md5_sum(bufz,nprec,resblock)
c     call MPI_FILE_SEEK(fhz_md5,disp_md5,MPI_SEEK_SET,err)
c     call MPI_FILE_WRITE(fhz_md5,resblock,32,MPI_CHARACTER,status,err)
c     end if

                         if ((count == write_step) .or. (it==(ntiskp*int(nt/ntiskp)))) then
                            if (ifault .le. 2) then
                               call MPI_FILE_WRITE_AT_ALL(fhz,offset,bufz,nprec*write_step,MPI_REAL,status,err)
                            else
                               if (it==(ntiskp*int(nt/ntiskp))) then
                                  call MPI_FILE_WRITE_AT_ALL(fhz,offset,bufz,nprec3,MPI_REAL,status,err)
                               end if
                            end if
                         end if
                      endif

                      ! -----------------------------------------------------------
                      ! close files and deallocate buffers if data not accumulated
                      ! -----------------------------------------------------------
                      if(ivelocity==1) then
                         if ((count == write_step) .or. (it==(ntiskp*int(nt/ntiskp)))) then
                            if ((ifault == 3 .or. ifault == 4) .and. it .ne. (ntiskp*int(nt/ntiskp))) then
                               call open3dp(0,MCR,fhx,fhy,fhz,sxgro_timestep,sygro_timestep,szgro_timestep,
     +                              fhx_md5,fhy_md5,fhz_md5,sxgro_md5,sygro_md5,
     +                              szgro_md5,imd5,ifault,0,filen)
                            else
                               call open3dp(0,MCR,fhx,fhy,fhz,sxgro_timestep,sygro_timestep,szgro_timestep,
     +                              fhx_md5,fhy_md5,fhz_md5,sxgro_md5,sygro_md5,
     +                              szgro_md5,imd5,ifault,1,filen)
                            endif

                            if (idyna==1 .and. ifault==4) then
                               call MPI_FILE_CLOSE(fhsgsn,err)
                            end if

                         endif
                      end if

                      if ((count == write_step) .or. (it==(ntiskp*int(nt/ntiskp)))) then
                         deallocate(bufx)
                         deallocate(bufy)
                         if ((ifault .le. 2) .or. (it==(ntiskp*int(nt/ntiskp)))) then
                            deallocate(bufz)
                         end if

                         if (idyna==1 .and. ifault==4) then
                            deallocate(bufxyz)
                         end if

                         count = 0
                      endif
                   endif ! if(i==it) then
               endif ! if(rank==recproc) then
!            endif ! if (sgt_io_out .eq. 1) then
        endif ! if (io_opt == 1) then

      ! --------------------------------------------------
      ! Write checkpoint files in case system crashes 
      ! during run
      ! --------------------------------------------------
      if(checkpoint.eq.0)  goto 333 
      if( mod( i, checkpoint ) .eq. 0 ) then
         t20 = t20 + 1.0
c     ckp_count = ckp_count + 1
         t18 = mpi_wtime()          

         ! --------------------------------------
         ! allow IOST number of readers
         ! --------------------------------------
         do iwrite=0,IOST
            if (mod(rank,IOST+1).eq.iwrite) then
               write( ckpfile, '(a,i7.7,i7.7,a)') 'output_ckp/ckp', rank, i, '.bin'
               if( rank .eq. 0 ) write( 0, * ) 'Writing checkpoint file'
               open( 9, file=ckpfile, form='unformatted', access='direct',
     $              recl=bsize*(nxt+4)*(nyt+4)*(nzt+4)*3 )
               write( 9, rec=1 ) u1, v1, w1
               write( 9, rec=2 ) xx, yy, zz
               write( 9, rec=3 ) xy, xz, yz
               write( 9, rec=4 ) r1, r2, r3
               write( 9, rec=5 ) r4, r5, r6
               close( 9 )
               
               write( ckpfile, '(a,i7.7,i7.7,a)') 'output_ckp/ckpdy', rank, i, '.bin'
               open( 9, file=ckpfile, form='unformatted', access='direct',
     $              recl=bsize*(nxt+4)*(nyt+4)*(nzt+4)*3 )
               write( 9, rec=1 ) u2,broken,mu
               close( 9 )
               
               if (npc==1) then
                  write( ckpfile, '(a,i7.7,i7.7,a)') 'output_ckp/ckpew', rank, i, '.bin'
              open( 9, file=ckpfile, form='unformatted', access='direct',
     $                 recl=bsize*ndt*nyt*nzt*6 )
                  write( 9, rec=1 ) ueT, ueH, veT, veH, weT, weH
                  write( 9, rec=2 ) xxeT, xxeH, yyeT, yyeH, zzeT, zzeH
                  write( 9, rec=3 ) xzeT, xzeH, yzeT, yzeH, xyeT, xyeH
                  write( 9, rec=4 ) uwT, uwH, vwT, vwH, wwT, wwH
                  write( 9, rec=5 ) xxwT, xxwH, yywT, yywH, zzwT, zzwH
                  write( 9, rec=6 ) xzwT, xzwH, yzwT, yzwH, xywT, xywH
                  close( 9 )
                  ! write(29,*) ckpfile
                  ! call flush(29)

                write( ckpfile, '(a,i7.7,i7.7,a)') 'output_ckp/ckpns', rank, i, '.bin'
              open( 9, file=ckpfile, form='unformatted', access='direct',
     $                 recl=bsize*nxt*ndt*nzt*6 )
                  write( 9, rec=1 ) unT, unH, vnT, vnH, wnT, wnH
                  write( 9, rec=2 ) xxnT, xxnH, yynT, yynH, zznT, zznH
                  write( 9, rec=3 ) xznT, xznH, yznT, yznH, xynT, xynH
                  write( 9, rec=4 ) usT, usH, vsT, vsH, wsT, wsH
                  write( 9, rec=5 ) xxsT, xxsH, yysT, yysH, zzsT, zzsH
                  write( 9, rec=6 ) xzsT, xzsH, yzsT, yzsH, xysT, xysH
                  close( 9 )
                  ! write(29,*) ckpfile
                  ! call flush(29)

                  write( ckpfile, '(a,i7.7,i7.7,a)') 'output_ckp/ckpb', rank, i, '.bin'
              open( 9, file=ckpfile, form='unformatted', access='direct',
     $                 recl=bsize*nxt*nyt*ndt*6 )
                  write( 9, rec=1 ) ubT, ubH, vbT, vbH, wbT, wbH
                  write( 9, rec=2 ) xxbT, xxbH, yybT, yybH, zzbT, zzbH
                  write( 9, rec=3 ) xzbT, xzbH, yzbT, yzbH, xybT, xybH
                  close( 9 )
                  ! write(29,*) ckpfile
                  ! call flush(29)
               end if           ! NPC
           write( ckpfile, '(a,i7.7,a)' ) 'output_ckp/ckp', rank, '.hdr'
               open( 9, file=ckpfile )
               write( 9, * ) i
               close( 9 )
            end if              ! IOST
            call MPI_BARRIER(MC1,err)
         end do                 ! iwrite

         t19=mpi_wtime()
         t21=t21+(t19-t18)
         ! write(29,*) ckpfile
         ! call flush(29)
      end if
 333  continue

      if (iperf == 1) then
         call MPI_BARRIER(MC1,err)
         t9=mpi_wtime()
         t10=t10+(t9-t7)
      end if
      
      ! ----------------------------------------------
      ! Get some more timing statistics
      ! ----------------------------------------------
      if ((iperf==1).and.(i== write_step*ntiskp)) then 
         dtime(6)=t8/t12
         dtime(8)=t10/t12
         do idtout=1,11
            dtout(idtout)=0.0d-00
         enddo
         call MPI_ALLREDUCE(dtime,dtout,11,MPI_REAL8,MPI_MAX,MC1,err)
         if (rank==master) then
            write (*,'(" computing time per time step =",f10.3," sec")') dtout(6)
            if (io_opt.ne.0) then
               write (*,'(" mpiio time per time step =",f10.3," sec")') dtout(8)
            endif
         endif
      end if
      end do                    ! END OF MAIN LOOP
      
C     ==============================================
C     CLOSE FILES WHEN DATA IS AGGREGATED
C     ==============================================
      if(ivelocity .ne. 1) then
         call open3dp(0,MCR,fhx,fhy,fhz,sxgro_timestep,sygro_timestep,szgro_timestep,
     +        fhx_md5,fhy_md5,fhz_md5,sxgro_md5,sygro_md5,
     +        szgro_md5,imd5,ifault,1,filen)

         if(isfcvlm==1 .or. ifault==4) then
            call open3dp2(0,MCR,fhx2,fhy2,fhz2,sxgro2_timestep,sygro2_timestep,szgro2_timestep,
     +           fhx2_md5,fhy2_md5,fhz2_md5,sxgro2_md5,sygro2_md5,
     +           szgro2_md5,imd5,ifault,1,filen)
         end if

         if (ifault==4 .and. idyna==1) then
            call MPI_FILE_CLOSE(fhsgsn,err)
         end if
      end if
      
C     ======================================
C     CLOSE REMAINING FILES 
C     ======================================
      if(rank==master) call open3d(0,nve,chkp,chkj,insrc,invel,filen,ifault,nxt,nst)
      if (iperf == 1) then
        call MPI_BARRIER(MC1,err)
      end if

!********************   MPI Output   **************************
      if (sgt_mpi_file .ne. 0) then
        call MPI_FILE_CLOSE(sgt_mpi_file,err)
      endif
!********************   MPI Output   **************************

C     ============================================
C     GET ELAPSED TIME, COMPUTATION TIME, MPI-IO
C     TIME
C     ============================================
      if (iperf==1) then
        t11=mpi_wtime()
        dtime(1)=t21/t20
c       dtime(1)=t21/real(ckp_count)
        dtime(2)=t23-t22
        dtime(3)=t15-t4
        dtime(4)=t5-t1
        dtime(5)=t3-t2
        dtime(6)=t8/t12
c       dtime(6)=t8/real(nt-t0)
        dtime(7)=t4-t1
        dtime(8)=t10/t12
c       dtime(8)=t10/real(nt-t0)
        dtime(9)=t11-t1
        dtime(10)=t14-t13
c       dtime(10)=t16/t12

c     -------------------------------------------------------
c     in the case of only slice, some nodes have only nprec=0
c     but master processor nprec.ne.0
c     -------------------------------------------------------
        xmbt= real(size1*4*nprec)/real(1024*1024)
        xmb2t= real(size1*4*nprec2)/real(1024*1024)
c       xmb=real(size1)*xmb1
c       print *, 'npts=',npts,'nprec=',nprec,'nprec2=',nprec2

        do idtout=1,11
          dtout(idtout)=0.0d-00
        enddo
c       print *, 'dtime=',dtime,rank
        call MPI_ALLREDUCE(dtime,dtout,11,MPI_REAL8,MPI_MAX,MC1,err)
        call MPI_ALLREDUCE(xmbt,xmb,1,MPI_REAL8,MPI_MAX,MC1,err)
        call MPI_ALLREDUCE(xmb2t,xmb2,1,MPI_REAL8,MPI_MAX,MC1,err)
c       if(rank==master) print *, 'dtout=',dtout

        if (rank==master) then
c         write (*,'("")')
          write (*,'(" Final---------")')
          write (*,'("    inialization time =",f10.3," sec")') dtout(4)
          write (*,'("       get station time =",f10.3," sec")') dtout(3)
          write (*,'("       read source time =",f10.3," sec")') dtout(7)
          write (*,'("       read media time  =",f10.3," sec")') dtout(5)
          write (*,'("    computing time per time step =",f10.3," sec")') dtout(6)  
          if (io_opt.ne.0) then
            write (*,'("    mpiio time per time step =",f10.3," sec")') dtout(8) 
          endif
          write (*,'("    total elapsed time is =",f10.3," sec")') dtout(9)
          if (checkpoint.ne.0) then 
            write (*,'("    checkpoint time on average =",f10.3," sec")') dtout(1)
          endif
c         if (dtout(2).gt.0.0001) then
c           write (*,'("single MPI write time (2) =",f6.2," sec")') dtout(2)
c           write (*,'("MPIIO Performance (2) =",f9.2,
c     $    " MB/sec")') xmb2/dtout(2)
c          endif
c          if (dtout(10).gt.0.0001) then
c            write (*,'("single MPI write time (1)=",f6.2," sec")') dtout(10)
c            write (*,'("MPIIO Performance (1) =",f9.2,
c     $    " MB/sec")') xmb/dtout(10)
c          endif
c     $    " MB/sec (total",f9.2," MB)")') xmb/dtout(10),xmb
c          write (*,'("MPIIO Performance per processor =",f9.2,
c     $    " MB/sec (",f9.2," MB/p)")') xmb1/dtout(10),xmb1
        endif
      end if

cPo SGT closing files
      if(igreen.ne.-1) then
         deallocate(sg1)
         deallocate(sg2)
         deallocate(sg3)
      endif
cSGT closed files
   
222   call MPI_FINALIZE(err)

      stop
      end
