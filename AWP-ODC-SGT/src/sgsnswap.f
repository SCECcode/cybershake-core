CCC   'sgsnswap.f' EXCHANGES SGSN VARIABLES BETWEEN PROCS OF DYNAMIC FAULT MODEL
CCC    of Dalguer and Day (2006)
CCC   
CCC           
CCC

      subroutine srcprmswp(rank,comm,nx,nz,fx_d,fz_d,ix_d,iz_d,
     +                    toe,fromw,tow,frome,tou,fromd,tod,fromu,
     +                    wu2ed,ed2wu,wd2eu,eu2wd)
     
     
CCC   EXCHANGES DYNAMIC SOURCE PARAMETERS BETWEEN PROCS

      use parstat

      integer :: rank,comm
      integer :: nx,nz,fx_d,fz_d,ix_d,iz_d

      integer :: fromw,toe,frome,tow
      integer :: fromd,tou,fromu,tod

      integer :: wu2ed,ed2wu,wd2eu,eu2wd
     

C     EAST AND WEST BORDERS OF FAULT

      call w2e_f(rank,comm,d0,ix_d,fx_d,fz_d,nx,nz,toe,fromw)
      call e2w_f(rank,comm,d0,ix_d,fx_d,fz_d,nx,nz,tow,frome)

      call w2e_f(rank,comm,mu_s,ix_d,fx_d,fz_d,nx,nz,toe,fromw)
      call e2w_f(rank,comm,mu_s,ix_d,fx_d,fz_d,nx,nz,tow,frome)            

      call w2e_f(rank,comm,mu_d,ix_d,fx_d,fz_d,nx,nz,toe,fromw)
      call e2w_f(rank,comm,mu_d,ix_d,fx_d,fz_d,nx,nz,tow,frome)      
      
      
C     UP AND DOWN BORDERS OF FAULT

      call d2u_f(rank,comm,d0,fx_d,iz_d,fz_d,nx,nz,tou,fromd)
      call u2d_f(rank,comm,d0,fx_d,iz_d,fz_d,nx,nz,tod,fromu)

      call d2u_f(rank,comm,mu_s,fx_d,iz_d,fz_d,nx,nz,tou,fromd)
      call u2d_f(rank,comm,mu_s,fx_d,iz_d,fz_d,nx,nz,tod,fromu)

      call d2u_f(rank,comm,mu_d,fx_d,iz_d,fz_d,nx,nz,tou,fromd)
      call u2d_f(rank,comm,mu_d,fx_d,iz_d,fz_d,nx,nz,tod,fromu)

C CORNERS OF THE FAULT
       
      call wu2ed_f(rank,comm,d0,ix_d,fx_d,iz_d,fz_d,nx,nz,wu2ed,ed2wu)
      call ed2wu_f(rank,comm,d0,ix_d,fx_d,iz_d,fz_d,nx,nz,ed2wu,wu2ed)

      call wu2ed_f(rank,comm,mu_s,ix_d,fx_d,iz_d,fz_d,nx,nz,wu2ed,ed2wu)
      call ed2wu_f(rank,comm,mu_s,ix_d,fx_d,iz_d,fz_d,nx,nz,ed2wu,wu2ed)

      call wu2ed_f(rank,comm,mu_d,ix_d,fx_d,iz_d,fz_d,nx,nz,wu2ed,ed2wu)
      call ed2wu_f(rank,comm,mu_d,ix_d,fx_d,iz_d,fz_d,nx,nz,ed2wu,wu2ed)

      
      call wd2eu_f(rank,comm,d0,ix_d,fx_d,iz_d,fz_d,nx,nz,wd2eu,eu2wd)
      call eu2wd_f(rank,comm,d0,ix_d,fx_d,iz_d,fz_d,nx,nz,eu2wd,wd2eu)

      call wd2eu_f(rank,comm,mu_s,ix_d,fx_d,iz_d,fz_d,nx,nz,wd2eu,eu2wd)
      call eu2wd_f(rank,comm,mu_s,ix_d,fx_d,iz_d,fz_d,nx,nz,eu2wd,wd2eu)

      call wd2eu_f(rank,comm,mu_d,ix_d,fx_d,iz_d,fz_d,nx,nz,wd2eu,eu2wd)
      call eu2wd_f(rank,comm,mu_d,ix_d,fx_d,iz_d,fz_d,nx,nz,eu2wd,wd2eu)


      return
      end      

      subroutine splitvelswp(rank,comm,nx,nz,fx_d,y_d,fz_d,ix_d,iz_d,
     +                    toe,fromw,tow,frome,tou,fromd,tod,fromu,
     +                    wu2ed,ed2wu,wd2eu,eu2wd)
     
     
CCC   EXCHANGES SPLIT VELOCITIES BETWEEN PROCS

      use parstat

      integer :: rank,comm
      integer :: nx,nz,fx_d,y_d,fz_d,ix_d,iz_d

      integer :: fromw,toe,frome,tow
      integer :: fromd,tou,fromu,tod

      integer :: wu2ed,ed2wu,wd2eu,eu2wd
     

C     EAST AND WEST BORDERS OF FAULT

      call w2e_f(rank,comm,uplus,ix_d,fx_d,fz_d,nx,nz,toe,fromw)
      call e2w_f(rank,comm,uplus,ix_d,fx_d,fz_d,nx,nz,tow,frome)

      call w2e_f(rank,comm,wplus,ix_d,fx_d,fz_d,nx,nz,toe,fromw)
      call e2w_f(rank,comm,wplus,ix_d,fx_d,fz_d,nx,nz,tow,frome)            
      
      
C     UP AND DOWN BORDERS OF FAULT

      call d2u_f(rank,comm,uplus,fx_d,iz_d,fz_d,nx,nz,tou,fromd)
      call u2d_f(rank,comm,uplus,fx_d,iz_d,fz_d,nx,nz,tod,fromu)

      call d2u_f(rank,comm,wplus,fx_d,iz_d,fz_d,nx,nz,tou,fromd)
      call u2d_f(rank,comm,wplus,fx_d,iz_d,fz_d,nx,nz,tod,fromu)

      
      
C CORNERS OF THE FAULT
       
      call wu2ed_f(rank,comm,uplus,ix_d,fx_d,iz_d,fz_d,nx,nz,wu2ed,ed2wu)
      call ed2wu_f(rank,comm,uplus,ix_d,fx_d,iz_d,fz_d,nx,nz,ed2wu,wu2ed)

      call wu2ed_f(rank,comm,wplus,ix_d,fx_d,iz_d,fz_d,nx,nz,wu2ed,ed2wu)
      call ed2wu_f(rank,comm,wplus,ix_d,fx_d,iz_d,fz_d,nx,nz,ed2wu,wu2ed)

      call wu2ed_f(rank,comm,u1(-1:nx+2,y_d,-1:nz+2),
     +            ix_d,fx_d,iz_d,fz_d,nx,nz,wu2ed,ed2wu)
      call ed2wu_f(rank,comm,u1(-1:nx+2,y_d,-1:nz+2),
     +            ix_d,fx_d,iz_d,fz_d,nx,nz,ed2wu,wu2ed)

      call wu2ed_f(rank,comm,w1(-1:nx+2,y_d,-1:nz+2),
     +            ix_d,fx_d,iz_d,fz_d,nx,nz,wu2ed,ed2wu)
      call ed2wu_f(rank,comm,w1(-1:nx+2,y_d,-1:nz+2),
     +            ix_d,fx_d,iz_d,fz_d,nx,nz,ed2wu,wu2ed)
      
      
      call wd2eu_f(rank,comm,uplus,ix_d,fx_d,iz_d,fz_d,nx,nz,wd2eu,eu2wd)
      call eu2wd_f(rank,comm,uplus,ix_d,fx_d,iz_d,fz_d,nx,nz,eu2wd,wd2eu)

      call wd2eu_f(rank,comm,wplus,ix_d,fx_d,iz_d,fz_d,nx,nz,wd2eu,eu2wd)
      call eu2wd_f(rank,comm,wplus,ix_d,fx_d,iz_d,fz_d,nx,nz,eu2wd,wd2eu)

      call wd2eu_f(rank,comm,u1(-1:nx+2,y_d,-1:nz+2),
     +             ix_d,fx_d,iz_d,fz_d,nx,nz,wd2eu,eu2wd)
      call eu2wd_f(rank,comm,u1(-1:nx+2,y_d,-1:nz+2),
     +             ix_d,fx_d,iz_d,fz_d,nx,nz,eu2wd,wd2eu)

      call wd2eu_f(rank,comm,w1(-1:nx+2,y_d,-1:nz+2),
     +             ix_d,fx_d,iz_d,fz_d,nx,nz,wd2eu,eu2wd)
      call eu2wd_f(rank,comm,w1(-1:nx+2,y_d,-1:nz+2),
     +             ix_d,fx_d,iz_d,fz_d,nx,nz,eu2wd,wd2eu)
     
      return
      end


      subroutine splitstrswp(rank,comm,nx,nz,fx_d,y_d,fz_d,ix_d,iz_d,
     +              toe,fromw,tow,frome,tou,fromd,tod,fromu)
      
     
CCC   EXCHANGES SPLIT STRESSES  BETWEEN PROCS

      use parstat

      integer :: rank,comm
      integer :: nx,nz,fx_d,y_d,fz_d,ix_d,iz_d

      integer :: fromw,toe,frome,tow
      integer :: fromd,tou,fromu,tod

     

C     EAST AND WEST BORDERS OF FAULT

      call w2e_f(rank,comm,xxplus,ix_d,fx_d,fz_d,nx,nz,toe,fromw)
      call e2w_f(rank,comm,xxplus,ix_d,fx_d,fz_d,nx,nz,tow,frome)

      call w2e_f(rank,comm,zzplus,ix_d,fx_d,fz_d,nx,nz,toe,fromw)
      call e2w_f(rank,comm,zzplus,ix_d,fx_d,fz_d,nx,nz,tow,frome)            
      
      call w2e_f(rank,comm,xzplus,ix_d,fx_d,fz_d,nx,nz,toe,fromw)
      call e2w_f(rank,comm,xzplus,ix_d,fx_d,fz_d,nx,nz,tow,frome)            
      
C     UP AND DOWN BORDERS OF FAULT

      call d2u_f(rank,comm,xxplus,fx_d,iz_d,fz_d,nx,nz,tou,fromd)
      call u2d_f(rank,comm,xxplus,fx_d,iz_d,fz_d,nx,nz,tod,fromu)

      call d2u_f(rank,comm,zzplus,fx_d,iz_d,fz_d,nx,nz,tou,fromd)
      call u2d_f(rank,comm,zzplus,fx_d,iz_d,fz_d,nx,nz,tod,fromu)

      call d2u_f(rank,comm,xzplus,fx_d,iz_d,fz_d,nx,nz,tou,fromd)
      call u2d_f(rank,comm,xzplus,fx_d,iz_d,fz_d,nx,nz,tod,fromu)
      
      
      return
      end

      subroutine tracswp(rank,comm,nx,nz,fx_d,y_d,fz_d,ix_d,iz_d,
     +                    toe,fromw,tow,frome,tou,fromd,tod,fromu,
     +                    wu2ed,ed2wu,wd2eu,eu2wd)
     
     
CCC   EXCHANGES TRACTIONS BETWEEN PROCS

      use parstat

      integer :: rank,comm
      integer :: nx,nz,fx_d,y_d,fz_d,ix_d,iz_d

      integer :: fromw,toe,frome,tow
      integer :: fromd,tou,fromu,tod

      integer :: wu2ed,ed2wu,wd2eu,eu2wd
     

C     EAST AND WEST BORDERS OF FAULT

      call w2e_f(rank,comm,tru1,ix_d,fx_d,fz_d,nx,nz,toe,fromw)
      call e2w_f(rank,comm,tru1,ix_d,fx_d,fz_d,nx,nz,tow,frome)

      call w2e_f(rank,comm,trw1,ix_d,fx_d,fz_d,nx,nz,toe,fromw)
      call e2w_f(rank,comm,trw1,ix_d,fx_d,fz_d,nx,nz,tow,frome)            
      
      
C     UP AND DOWN BORDERS OF FAULT

      call d2u_f(rank,comm,tru1,fx_d,iz_d,fz_d,nx,nz,tou,fromd)
      call u2d_f(rank,comm,tru1,fx_d,iz_d,fz_d,nx,nz,tod,fromu)

      call d2u_f(rank,comm,trw1,fx_d,iz_d,fz_d,nx,nz,tou,fromd)
      call u2d_f(rank,comm,trw1,fx_d,iz_d,fz_d,nx,nz,tod,fromu)

      
      
C CORNERS OF THE FAULT
       
      call wu2ed_f(rank,comm,tru1,ix_d,fx_d,iz_d,fz_d,nx,nz,wu2ed,ed2wu)
      call ed2wu_f(rank,comm,tru1,ix_d,fx_d,iz_d,fz_d,nx,nz,ed2wu,wu2ed)

      call wu2ed_f(rank,comm,trw1,ix_d,fx_d,iz_d,fz_d,nx,nz,wu2ed,ed2wu)
      call ed2wu_f(rank,comm,trw1,ix_d,fx_d,iz_d,fz_d,nx,nz,ed2wu,wu2ed)
      
      
      call wd2eu_f(rank,comm,tru1,ix_d,fx_d,iz_d,fz_d,nx,nz,wd2eu,eu2wd)
      call eu2wd_f(rank,comm,tru1,ix_d,fx_d,iz_d,fz_d,nx,nz,eu2wd,wd2eu)

      call wd2eu_f(rank,comm,trw1,ix_d,fx_d,iz_d,fz_d,nx,nz,wd2eu,eu2wd)
      call eu2wd_f(rank,comm,trw1,ix_d,fx_d,iz_d,fz_d,nx,nz,eu2wd,wd2eu)

      return
      end
      
     
      subroutine w2e_f(rank,comm,b,ix_d,fx_d,fz_d,nx,nz,dst,src)

CCC   TRANSFER FAULT BORDERS FROM W TO E PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,nz,fx_d,fz_d,ix_d
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(-1:nx+2,-1:nz+2) :: b

      real, dimension(fx_d-1:fx_d,1:fz_d) :: bto
      real, dimension(-1:0,1:fz_d) :: bfrom

      npts = 2*fz_d

      if(dst /= MPI_PROC_NULL.and.fx_d==nx) then
        bto = b(fx_d-1:fx_d,1:fz_d)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL.and.ix_d==1) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(-1:0,1:fz_d) = bfrom
      end if

      return
      end
      

      subroutine e2w_f(rank,comm,b,ix_d,fx_d,fz_d,nx,nz,dst,src)

CCC   TRANSFER FAULT BORDERS FROM E TO W PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,nz,fx_d,fz_d,ix_d
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(-1:nx+2,-1:nz+2) :: b

      real, dimension(1:2,1:fz_d) :: bto
      real, dimension(fx_d+1:fx_d+2,1:fz_d) :: bfrom

      npts = 2*fz_d

      if(dst /= MPI_PROC_NULL.and.ix_d==1) then
        bto = b(1:2,1:fz_d)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL.and.fx_d==nx) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(fx_d+1:fx_d+2,1:fz_d) = bfrom
      end if

      return
      end
      

      subroutine d2u_f(rank,comm,b,fx_d,iz_d,fz_d,nx,nz,dst,src)

CCC   TRANSFER FAULT BORDERS FROM DOWN TO UP PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,nz,fx_d,fz_d,iz_d
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(-1:nx+2,-1:nz+2) :: b

      real, dimension(1:fx_d,fz_d-1:fz_d) :: bto
      real, dimension(1:fx_d,-1:0) :: bfrom

      npts = 2*fx_d

      if(dst /= MPI_PROC_NULL.and.fz_d==nz) then
        bto = b(1:fx_d,fz_d-1:fz_d)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL.and.iz_d==1) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(1:fx_d,-1:0) = bfrom
      end if

      return
      end


      subroutine u2d_f(rank,comm,b,fx_d,iz_d,fz_d,nx,nz,dst,src)

CCC   TRANSFER FAULT BORDERS FROM UP TO DOWN PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,nz,fx_d,fz_d,iz_d
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(-1:nx+2,-1:nz+2) :: b

      real, dimension(1:fx_d,1:2) :: bto
      real, dimension(1:fx_d,fz_d+1:fz_d+2) :: bfrom

      npts = 2*fx_d

      if(dst /= MPI_PROC_NULL.and.iz_d==1) then
        bto = b(1:fx_d,1:2)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL.and.fz_d==nz) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(1:fx_d,fz_d+1:fz_d+2) = bfrom
      end if

      return
      end

      subroutine wu2ed_f(rank,comm,b,ix_d,fx_d,iz_d,fz_d,nx,nz,dst,src)      

CCC   TRANSFER FAULT CORNERS FROM W-UP TO E-DOWN PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,nz,fx_d,fz_d,iz_d,ix_d
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(-1:nx+2,-1:nz+2) :: b

      real, dimension(fx_d-1:fx_d,1:2) :: bto
      real, dimension(-1:0,fz_d+1:fz_d+2) :: bfrom

      npts = 4

      if(dst /= MPI_PROC_NULL.and.fx_d==nx.and.iz_d==1) then
        bto = b(fx_d-1:fx_d,1:2)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL.and.ix_d==1.and.fz_d==nz) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(-1:0,fz_d+1:fz_d+2) = bfrom
      end if

      return
      end


      subroutine ed2wu_f(rank,comm,b,ix_d,fx_d,iz_d,fz_d,nx,nz,dst,src)      

CCC   TRANSFER FAULT CORNERS FROM E-DOWN TO W-UP PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,nz,fx_d,fz_d,iz_d,ix_d
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(-1:nx+2,-1:nz+2) :: b

      real, dimension(1:2,fz_d-1:fz_d) :: bto
      real, dimension(fx_d+1:fx_d+2,-1:0) :: bfrom

      npts = 4

      if(dst /= MPI_PROC_NULL.and.ix_d==1.and.fz_d==nz) then
        bto = b(1:2,fz_d-1:fz_d)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL.and.fx_d==nx.and.iz_d==1) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(fx_d+1:fx_d+2,-1:0) = bfrom
      end if

      return
      end
      
      
      subroutine wd2eu_f(rank,comm,b,ix_d,fx_d,iz_d,fz_d,nx,nz,dst,src)      

CCC   TRANSFER FAULT CORNERS FROM W-DOWN TO E-UP PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,nz,fx_d,fz_d,iz_d,ix_d
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(-1:nx+2,-1:nz+2) :: b

      real, dimension(fx_d-1:fx_d,fz_d-1:fz_d) :: bto
      real, dimension(-1:0,-1:0) :: bfrom

      npts = 4

      if(dst /= MPI_PROC_NULL.and.fx_d==nx.and.fz_d==nz) then
        bto = b(fx_d-1:fx_d,fz_d-1:fz_d)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL.and.ix_d==1.and.iz_d==1) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(-1:0,-1:0) = bfrom
      end if

      return
      end


      subroutine eu2wd_f(rank,comm,b,ix_d,fx_d,iz_d,fz_d,nx,nz,dst,src)      

CCC   TRANSFER FAULT CORNERS FROM E-UP TO W-DOWN PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,nz,fx_d,fz_d,iz_d,ix_d
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(-1:nx+2,-1:nz+2) :: b

      real, dimension(1:2,1:2) :: bto
      real, dimension(fx_d+1:fx_d+2,fz_d+1:fz_d+2) :: bfrom

      npts = 4

      if(dst /= MPI_PROC_NULL.and.ix_d==1.and.iz_d==1) then
        bto = b(1:2,1:2)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL.and.fx_d==nx.and.fz_d==nz) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(fx_d+1:fx_d+2,fz_d+1:fz_d+2) = bfrom
      end if

      return
      end
     
