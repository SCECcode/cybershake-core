CCC   'memory.f' ANNOUNCES AND ALLOCATES MEMORY FOR ARRAYS

CC
      module parstat 
      include "mpif.h"

C     32-bit or 64-bit depending on the machine

      integer,parameter :: b8  = selected_real_kind(15) 
C     On 64-bit machine, using
      integer,parameter :: i8 = selected_int_kind(10)
C     On 32-bit machine, using 
c      integer,parameter :: i8 = selected_int_kind(7)

      integer, parameter :: kblock = 16
      integer, parameter :: jblock = 8
      integer :: requests_v(36), num_of_requests_v
CCC   GLOBAL FIELD VARIABLE AND PML FIELD VARIABLE MEMORY MODULE

C     SRC, RECVR, BOUNDS, WAVEFIELD VARIABLES AND MEDIA PARAMETERS

C      integer, dimension (:,:), allocatable :: tpsrc,psrc
      integer, dimension (:,:), allocatable :: tprec,prec,prec2,prec3,prec4,sgt_prec
      integer :: sgt_io_out,sgt_numsta,sgt_comm
      integer, dimension (:,:), allocatable :: sgt_sta_coord, sgt_local_sta_coord
      real, dimension (:,:), allocatable :: taxx,tayy,tazz,taxz,tayz,taxy
      real, dimension (:,:), allocatable :: axx,ayy,azz,axz,ayz,axy
      integer(kind=MPI_ADDRESS_KIND), dimension (:), allocatable :: tpmap,pmap,pmap2,pmap3,pmap4
      real, dimension (:), allocatable :: strisxx,strinxx,dmaxx
      real, dimension (:), allocatable :: tstrisxx,tstrinxx,tdmaxx
      real, dimension (:,:,:), allocatable :: u2,u3,rupt
      integer, dimension (:,:,:), allocatable :: broken

      integer, dimension (:), allocatable :: xi,xf,yi,yf,zi,zf
 
      real, dimension (:,:,:), allocatable :: u1,v1,w1
      real, dimension (:,:,:), allocatable :: xx,yy,zz,xy,xz,yz
      real, dimension (:,:,:), allocatable :: r1,r2,r3,r4,r5,r6
      real, dimension (:,:,:), allocatable :: qp,qs
      real, dimension (:,:,:), allocatable :: tmpvp,tmpvs,tmpdd,tmppq,tmpsq

      real, dimension (:,:,:), allocatable :: d1,mu,lam

      real, dimension (2,2,2) :: tau

C     PML DAMPING AND FIELD VARIABLES
 
      real, dimension (:), allocatable :: rh,rf
      real :: pht !damping ratio

      real, dimension (:,:,:), allocatable :: usT,usH,vsT,vsH,wsT,wsH
      real, dimension (:,:,:), allocatable :: uwT,uwH,vwT,vwH,wwT,wwH
      real, dimension (:,:,:), allocatable :: unT,unH,vnT,vnH,wnT,wnH
      real, dimension (:,:,:), allocatable :: ueT,ueH,veT,veH,weT,weH
      real, dimension (:,:,:), allocatable :: ubT,ubH,vbT,vbH,wbT,wbH

      real, dimension (:,:,:), allocatable :: xxsT,xxsH,yysT,yysH,zzsT,zzsH
      real, dimension (:,:,:), allocatable :: xxwT,xxwH,yywT,yywH,zzwT,zzwH
      real, dimension (:,:,:), allocatable :: xxnT,xxnH,yynT,yynH,zznT,zznH
      real, dimension (:,:,:), allocatable :: xxeT,xxeH,yyeT,yyeH,zzeT,zzeH
      real, dimension (:,:,:), allocatable :: xxbT,xxbH,yybT,yybH,zzbT,zzbH

      real, dimension (:,:,:), allocatable :: xzsT,xzsH,yzsT,yzsH,xysT,xysH
      real, dimension (:,:,:), allocatable :: xzwT,xzwH,yzwT,yzwH,xywT,xywH
      real, dimension (:,:,:), allocatable :: xznT,xznH,yznT,yznH,xynT,xynH
      real, dimension (:,:,:), allocatable :: xzeT,xzeH,yzeT,yzeH,xyeT,xyeH
      real, dimension (:,:,:), allocatable :: xzbT,xzbH,yzbT,yzbH,xybT,xybH

C     CERJAN DAMPING ARRAYS

      real, dimension (:,:,:), allocatable :: dcrjx, dcrjy, dcrjz
      integer :: flag_cerjanx, flag_cerjany, flag_cerjanz
      integer :: xs_crj(3),xe_crj(3),ys_crj(3),ye_crj(3),zs_crj(3),ze_crj(3) 
C
CC    FOR SGSN DYNAMIC FAULT MODEL
C
      integer, dimension (:,:), allocatable :: brokenu,brokenw
      integer, dimension (:,:), allocatable :: tpsrc,psrc      
      
      real, dimension (:,:), allocatable :: trupt,strisx,strisz,striny
      real, dimension (:), allocatable :: tstrisx,tstrisz,tstriny      
      real, dimension (:,:), allocatable :: tru1,tru2,trw1,trw2      
      real, dimension (:,:), allocatable :: d_u,d_w,d0,mu_d,mu_s      
      real, dimension (:), allocatable :: tdmax,tmu_d,tmu_s            
      real, dimension (:,:), allocatable :: rateu_u,rateu_w,rateu       
      real, dimension (:,:), allocatable :: ratew_u,ratew_w,ratew      
      real, dimension (:,:), allocatable :: Rplusu,Rplusw,Rminusu,Rminusw
      real, dimension (:,:), allocatable :: Rplusue,Rpluswe,Rminusue,Rminuswe      
      real, dimension (:,:), allocatable :: uplus,wplus,xxplus,zzplus,xzplus

      real, dimension (:,:,:), allocatable :: u_f,v_f,w_f
      real, dimension (:,:,:), allocatable :: xx_f,yy_f,zz_f,xy_f,xz_f,yz_f

      real, dimension (:,:), allocatable :: slipu_u,slipu_w
      real, dimension (:,:), allocatable :: trv1
      real, dimension (:,:), allocatable :: vplus,wplusi
      real, dimension (:,:), allocatable :: umin,vmin,wmin
      real, dimension (:,:), allocatable :: duplus,dvplus,dwplus
      real, dimension (:,:), allocatable :: dumin,dvmin,dwmin
c
      end module parstat
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
     
      subroutine bndmem(maxzone)

CCC   ALLOCATE MEMORY FOR COMPUTATIONAL ZONE BOUNDS

      use parstat

      integer :: i,maxzone

      allocate(xi(maxzone),xf(maxzone))
      allocate(yi(maxzone),yf(maxzone))
      allocate(zi(maxzone),zf(maxzone))

C     SET INITIAL BOUNDS TO ALLOW NO LOOPING - KEY CONCEPT
C     EVERY PROC CONSIDERS EVERY ZONE - UNEEDED ZONES KEEP
C     BOUNDS THAT DON"T ALLOW LOOPING

      do i=1,maxzone
        xi(i)=2
        yi(i)=2
        zi(i)=2
        xf(i)=1
        yf(i)=1
        zf(i)=1
      end do

      return
      end




      subroutine vsemem(nve,nxt,nyt,nzt)

CCC   DYNAMICALLY ALLOCATE FIELD VARIABLE AND MEDIA MEMORY

      use parstat

      integer :: nxt,nyt,nzt,nve

C     SET TO OVERLAP TWO NODES BETWEEN PROCS
 
      allocate(u1(-1:nxt+2,-1:nyt+2,-1:nzt+2))
      allocate(v1(-1:nxt+2,-1:nyt+2,-1:nzt+2))
      allocate(w1(-1:nxt+2,-1:nyt+2,-1:nzt+2))

      allocate(xx(-1:nxt+2,-1:nyt+2,-1:nzt+2))
      allocate(yy(-1:nxt+2,-1:nyt+2,-1:nzt+2))
      allocate(zz(-1:nxt+2,-1:nyt+2,-1:nzt+2))
      allocate(xy(-1:nxt+2,-1:nyt+2,-1:nzt+2))
      allocate(yz(-1:nxt+2,-1:nyt+2,-1:nzt+2))
      allocate(xz(-1:nxt+2,-1:nyt+2,-1:nzt+2))

      if (ifault==3 .or. ifault==4) then
        allocate(u2(-1:nxt+2,-1:nyt+2,-1:nzt+2))
        allocate(u3(-1:nxt+2,-1:nyt+2,-1:nzt+2))
        allocate(rupt(-1:nxt+2,-1:nyt+2,-1:nzt+2))
        allocate(broken(-1:nxt+2,-1:nyt+2,-1:nzt+2))
        u2=0
        u3=0
        rupt=0
        broken=0
      end if

      u1=0
      v1=0
      w1=0

      xx=0
      yy=0
      zz=0
      xy=0
      yz=0
      xz=0

      if(nve==1) then
        allocate(r1(1:nxt,1:nyt,1:nzt))
        allocate(r2(1:nxt,1:nyt,1:nzt))
        allocate(r3(1:nxt,1:nyt,1:nzt))
        allocate(r4(1:nxt,1:nyt,1:nzt))
        allocate(r5(1:nxt,1:nyt,1:nzt))
        allocate(r6(1:nxt,1:nyt,1:nzt))

        allocate(qp(0:nxt+1,0:nyt+1,0:nzt+1))
        allocate(qs(0:nxt+1,0:nyt+1,0:nzt+1))

        r1=0
        r2=0
        r3=0
        r4=0
        r5=0
        r6=0

        qp=0
        qs=0
      end if

C     IN ELASTIC PARAMETERS TOO FOR MEDIA AVERAGING

      allocate(d1(0:nxt+1,0:nyt+1,0:nzt+1))
      allocate(mu(0:nxt+1,0:nyt+1,0:nzt+1))
      allocate(lam(0:nxt+1,0:nyt+1,0:nzt+1))

      d1=0
      mu=0
      lam=0

      return
      end 


 
      subroutine pmlmem(dims,coords,maxdim,nxt,nyt,nzt,ndt,fs,fw,fn,fe,fb)

CCC   DYNAMICALLY ALLOCATE PML VARIABLE MEMORY

      use parstat

      integer :: i,j,k

      integer :: nxt,nyt,nzt,ndt
      integer :: nxp,nyp,xc,yc,zc

      integer :: fw,fs,fb,fe,fn

      integer, dimension(maxdim) :: dims,coords

C     NUMBER OF PROCS IN EACH DIMENSION AND COORDINATE OF PROC

      nxp = dims(1)
      nyp = dims(2)

      xc = coords(1)
      yc = coords(2)
      zc = coords(3)

C     ALLOCATION CALLS FOR EACH REQUIRED PLANE (S,W,N,E,B), AND
C     SET FLAGS FOR EACH PLANE (ON=1, INITIALIZED OFF=0) 

      call allocs(nxt,nzt,ndt)
      if(yc==0) then
        fs=1
      end if

      call allocw(nyt,nzt,ndt)
      if(xc==0) then
        fw=1
      end if
       
      call allocn(nxt,nzt,ndt)
      if(mod(yc+1,nyp)==0) then
        fn=1
      end if

      call alloce(nyt,nzt,ndt)
      if(mod(xc+1,nxp)==0) then
        fe=1
      end if

      call allocb(nxt,nyt,ndt)
      if(zc==0) then
        fb=1
      end if

C     ALLOCATE MEMORY FOR DAMPING ARRAYS

      allocate(rh(ndt))
      allocate(rf(ndt))

      rh=0
      rf=0

      return
      end



      subroutine allocs(nxt,nzt,ndt)

CCC   ALLOCATE S-PLANE MEMORY

      use parstat

      integer :: nxt,nyt,nzt,ndt

      allocate(usT(nxt,ndt,nzt),usH(nxt,ndt,nzt))
      allocate(vsT(nxt,ndt,nzt),vsH(nxt,ndt,nzt))
      allocate(wsT(nxt,ndt,nzt),wsH(nxt,ndt,nzt))  
      allocate(xxsT(nxt,ndt,nzt),xxsH(nxt,ndt,nzt))
      allocate(yysT(nxt,ndt,nzt),yysH(nxt,ndt,nzt))
      allocate(zzsT(nxt,ndt,nzt),zzsH(nxt,ndt,nzt))
      allocate(xzsT(nxt,ndt,nzt),xzsH(nxt,ndt,nzt))
      allocate(yzsT(nxt,ndt,nzt),yzsH(nxt,ndt,nzt))
      allocate(xysT(nxt,ndt,nzt),xysH(nxt,ndt,nzt))

      usT=0
      usH=0
      vsT=0
      vsH=0
      wsT=0
      wsH=0

      xxsT=0
      xxsH=0
      yysT=0
      yysH=0
      zzsT=0
      zzsH=0

      xzsT=0
      xzsH=0
      yzsT=0
      yzsH=0
      xysT=0
      xysH=0
  
      return 
      end



      subroutine allocw(nyt,nzt,ndt)

CCC   ALLOCATE W-PLANE MEMORY

      use parstat

      integer :: nxt,nyt,nzt,ndt

      allocate(uwT(ndt,nyt,nzt),uwH(ndt,nyt,nzt))
      allocate(vwT(ndt,nyt,nzt),vwH(ndt,nyt,nzt))
      allocate(wwT(ndt,nyt,nzt),wwH(ndt,nyt,nzt))
      allocate(xxwT(ndt,nyt,nzt),xxwH(ndt,nyt,nzt))
      allocate(yywT(ndt,nyt,nzt),yywH(ndt,nyt,nzt))
      allocate(zzwT(ndt,nyt,nzt),zzwH(ndt,nyt,nzt))
      allocate(xzwT(ndt,nyt,nzt),xzwH(ndt,nyt,nzt))
      allocate(yzwT(ndt,nyt,nzt),yzwH(ndt,nyt,nzt))
      allocate(xywT(ndt,nyt,nzt),xywH(ndt,nyt,nzt))

      uwT=0
      uwH=0
      vwT=0
      vwH=0
      wwT=0
      wwH=0

      xxwT=0
      xxwH=0
      yywT=0
      yywH=0
      zzwT=0
      zzwH=0

      xzwT=0
      xzwH=0
      yzwT=0
      yzwH=0
      xywT=0
      xywH=0

      return
      end



      subroutine allocn(nxt,nzt,ndt)

CCC   ALLOCATE N-PLANE MEMORY

      use parstat

      integer :: nxt,nzt,ndt

      allocate(unT(nxt,ndt,nzt),unH(nxt,ndt,nzt))
      allocate(vnT(nxt,ndt,nzt),vnH(nxt,ndt,nzt))
      allocate(wnT(nxt,ndt,nzt),wnH(nxt,ndt,nzt))  
      allocate(xxnT(nxt,ndt,nzt),xxnH(nxt,ndt,nzt))
      allocate(yynT(nxt,ndt,nzt),yynH(nxt,ndt,nzt))
      allocate(zznT(nxt,ndt,nzt),zznH(nxt,ndt,nzt))
      allocate(xznT(nxt,ndt,nzt),xznH(nxt,ndt,nzt))
      allocate(yznT(nxt,ndt,nzt),yznH(nxt,ndt,nzt))
      allocate(xynT(nxt,ndt,nzt),xynH(nxt,ndt,nzt))

      unT=0
      unH=0
      vnT=0
      vnH=0
      wnT=0
      wnH=0

      xxnT=0
      xxnH=0
      yynT=0
      yynH=0
      zznT=0
      zznH=0
      
      xznT=0
      xznH=0
      yznT=0
      yznH=0
      xynT=0
      xynH=0
  
      return 
      end



      subroutine alloce(nyt,nzt,ndt)

CCC   ALLOCATE E-PLANE MEMORY

      use parstat

      integer :: nxt,nyt,nzt,ndt

      allocate(ueT(ndt,nyt,nzt),ueH(ndt,nyt,nzt))
      allocate(veT(ndt,nyt,nzt),veH(ndt,nyt,nzt))
      allocate(weT(ndt,nyt,nzt),weH(ndt,nyt,nzt))
      allocate(xxeT(ndt,nyt,nzt),xxeH(ndt,nyt,nzt))
      allocate(yyeT(ndt,nyt,nzt),yyeH(ndt,nyt,nzt))
      allocate(zzeT(ndt,nyt,nzt),zzeH(ndt,nyt,nzt))
      allocate(xzeT(ndt,nyt,nzt),xzeH(ndt,nyt,nzt))
      allocate(yzeT(ndt,nyt,nzt),yzeH(ndt,nyt,nzt))
      allocate(xyeT(ndt,nyt,nzt),xyeH(ndt,nyt,nzt))

      ueT=0
      ueH=0
      veT=0
      veH=0
      weT=0
      weH=0

      xxeT=0
      xxeH=0
      yyeT=0
      yyeH=0
      zzeT=0
      zzeH=0

      xzeT=0
      xzeH=0
      yzeT=0
      yzeH=0
      xyeT=0
      xyeH=0

      return
      end



      subroutine allocb(nxt,nyt,ndt)

CCC   ALLOCATE B-PLANE MEMORY

      use parstat

      integer :: nxt,nyt,nzt,ndt

      allocate(ubT(nxt,nyt,ndt),ubH(nxt,nyt,ndt))
      allocate(vbT(nxt,nyt,ndt),vbH(nxt,nyt,ndt))
      allocate(wbT(nxt,nyt,ndt),wbH(nxt,nyt,ndt))
      allocate(xxbT(nxt,nyt,ndt),xxbH(nxt,nyt,ndt))
      allocate(yybT(nxt,nyt,ndt),yybH(nxt,nyt,ndt))
      allocate(zzbT(nxt,nyt,ndt),zzbH(nxt,nyt,ndt))
      allocate(xzbT(nxt,nyt,ndt),xzbH(nxt,nyt,ndt))
      allocate(yzbT(nxt,nyt,ndt),yzbH(nxt,nyt,ndt))
      allocate(xybT(nxt,nyt,ndt),xybH(nxt,nyt,ndt))

      ubT=0
      ubH=0
      vbT=0
      vbH=0
      wbT=0
      wbH=0

      xxbT=0
      xxbH=0
      yybT=0
      yybH=0
      zzbT=0
      zzbH=0

      xzbT=0
      xzbH=0
      yzbT=0
      yzbH=0
      xybT=0
      xybH=0

      return
      end




