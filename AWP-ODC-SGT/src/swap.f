CCC   'swap.f' EXCHANGES VELOCITY/STRESS PLANES BETWEEN PROCS

CCC   PERFORMANCE ISSUES HERE - OPTIMIZATION IN PROGRESS
CCC   SOME SUBROUTINE DESCRIPTIONS INCORRECT
      subroutine mediaswp(rank,nve,comm,nx,ny,nz,ton,froms,tos,fromn,
     +                    toe,fromw,tow,frome,tou,fromd,tod,fromu,
     +                    swd2neu,neu2swd,swu2ned,ned2swu,
     +                    sed2nwu,nwu2sed,seu2nwd,nwd2seu,
     +                    sw2ne,ne2sw,se2nw,nw2se,su2nd,nd2su,
     +                    sd2nu,nu2sd,wu2ed,ed2wu,wd2eu,eu2wd)

CCC   EXCHANGES MEDIA VALUES BETWEEN PROCS

      use parstat

      integer :: rank,comm,nve
      integer :: nx,ny,nz

      integer :: fromw,toe,frome,tow
      integer :: froms,ton,fromn,tos
      integer :: fromd,tou,fromu,tod

      integer :: swd2neu,neu2swd,swu2ned,ned2swu
      integer :: sed2nwu,nwu2sed,seu2nwd,nwd2seu

      integer :: sw2ne,ne2sw,se2nw,nw2se,su2nd,nd2su
      integer :: sd2nu,nu2sd,wu2ed,ed2wu,wd2eu,eu2wd
      

C     NORTH AND SOUTH PLANES

      call s2n(rank,comm,d1,nx,ny,nz,ton,froms,10,1,1)
      call n2s(rank,comm,d1,nx,ny,nz,tos,fromn,10,1,1)

      call s2n(rank,comm,mu,nx,ny,nz,ton,froms,11,1,1)
      call n2s(rank,comm,mu,nx,ny,nz,tos,fromn,11,1,1)

      call s2n(rank,comm,lam,nx,ny,nz,ton,froms,12,1,1)
      call n2s(rank,comm,lam,nx,ny,nz,tos,fromn,12,1,1)

      if (nve==1) then
        call s2n(rank,comm,qp,nx,ny,nz,ton,froms,13,1,1)
        call n2s(rank,comm,qp,nx,ny,nz,tos,fromn,13,1,1)

        call s2n(rank,comm,qs,nx,ny,nz,ton,froms,14,1,1)
        call n2s(rank,comm,qs,nx,ny,nz,tos,fromn,14,1,1)
      end if
      


C     EAST AND WEST PLANES

      call w2e(rank,comm,d1,nx,ny,nz,toe,fromw,10,1,1)
      call e2w(rank,comm,d1,nx,ny,nz,tow,frome,10,1,1)

      call w2e(rank,comm,mu,nx,ny,nz,toe,fromw,11,1,1)
      call e2w(rank,comm,mu,nx,ny,nz,tow,frome,11,1,1)

      call w2e(rank,comm,lam,nx,ny,nz,toe,fromw,12,1,1)
      call e2w(rank,comm,lam,nx,ny,nz,tow,frome,12,1,1)

      if (nve==1) then
        call w2e(rank,comm,qp,nx,ny,nz,toe,fromw,13,1,1)
        call e2w(rank,comm,qp,nx,ny,nz,tow,frome,13,1,1)

        call w2e(rank,comm,qs,nx,ny,nz,toe,fromw,14,1,1)
        call e2w(rank,comm,qs,nx,ny,nz,tow,frome,14,1,1)
      end if
      


C     UP AND DOWN PLANES

      call d2u(rank,comm,d1,nx,ny,nz,tou,fromd,10,1,1)
      call u2d(rank,comm,d1,nx,ny,nz,tod,fromu,10,1,1)

      call d2u(rank,comm,mu,nx,ny,nz,tou,fromd,11,1,1)
      call u2d(rank,comm,mu,nx,ny,nz,tod,fromu,12,1,1)

      call d2u(rank,comm,lam,nx,ny,nz,tou,fromd,12,1,1)
      call u2d(rank,comm,lam,nx,ny,nz,tod,fromu,12,1,1)

      if (nve==1) then
        call d2u(rank,comm,qp,nx,ny,nz,tou,fromd,13,1,1)
        call u2d(rank,comm,qp,nx,ny,nz,tod,fromu,13,1,1)

        call d2u(rank,comm,qs,nx,ny,nz,tou,fromd,14,1,1)
        call u2d(rank,comm,qs,nx,ny,nz,tod,fromu,14,1,1)
      end if

C     CORNER CUBES

      call swd2neu1(rank,comm,d1,nx,ny,nz,swd2neu,neu2swd)
      call neu2swd1(rank,comm,d1,nx,ny,nz,neu2swd,swd2neu)

      call swd2neu1(rank,comm,mu,nx,ny,nz,swd2neu,neu2swd)
      call neu2swd1(rank,comm,mu,nx,ny,nz,neu2swd,swd2neu)

      call swd2neu1(rank,comm,lam,nx,ny,nz,swd2neu,neu2swd)
      call neu2swd1(rank,comm,lam,nx,ny,nz,neu2swd,swd2neu)



      if (nve==1) then
        call swd2neu1(rank,comm,qp,nx,ny,nz,swd2neu,neu2swd)
        call neu2swd1(rank,comm,qp,nx,ny,nz,neu2swd,swd2neu)

        call swd2neu1(rank,comm,qs,nx,ny,nz,swd2neu,neu2swd)
        call neu2swd1(rank,comm,qs,nx,ny,nz,neu2swd,swd2neu)
      end if

      call swu2ned1(rank,comm,d1,nx,ny,nz,swu2ned,ned2swu)
      call ned2swu1(rank,comm,d1,nx,ny,nz,ned2swu,swu2ned) 

      call swu2ned1(rank,comm,mu,nx,ny,nz,swu2ned,ned2swu)
      call ned2swu1(rank,comm,mu,nx,ny,nz,ned2swu,swu2ned)

      call swu2ned1(rank,comm,lam,nx,ny,nz,swu2ned,ned2swu)
      call ned2swu1(rank,comm,lam,nx,ny,nz,ned2swu,swu2ned)


      if (nve==1) then
        call swu2ned1(rank,comm,qp,nx,ny,nz,swu2ned,ned2swu)
        call ned2swu1(rank,comm,qp,nx,ny,nz,ned2swu,swu2ned)

        call swu2ned1(rank,comm,qs,nx,ny,nz,swu2ned,ned2swu)
        call ned2swu1(rank,comm,qs,nx,ny,nz,ned2swu,swu2ned)
      end if

      call sed2nwu1(rank,comm,d1,nx,ny,nz,sed2nwu,nwu2sed)
      call nwu2sed1(rank,comm,d1,nx,ny,nz,nwu2sed,sed2nwu)

      call sed2nwu1(rank,comm,mu,nx,ny,nz,sed2nwu,nwu2sed)
      call nwu2sed1(rank,comm,mu,nx,ny,nz,nwu2sed,sed2nwu)

      call sed2nwu1(rank,comm,lam,nx,ny,nz,sed2nwu,nwu2sed)
      call nwu2sed1(rank,comm,lam,nx,ny,nz,nwu2sed,sed2nwu)

      if (nve==1) then
        call sed2nwu1(rank,comm,qp,nx,ny,nz,sed2nwu,nwu2sed)
        call nwu2sed1(rank,comm,qp,nx,ny,nz,nwu2sed,sed2nwu)

        call sed2nwu1(rank,comm,qs,nx,ny,nz,sed2nwu,nwu2sed)
        call nwu2sed1(rank,comm,qs,nx,ny,nz,nwu2sed,sed2nwu)
      end if

      call seu2nwd1(rank,comm,d1,nx,ny,nz,seu2nwd,nwd2seu)
      call nwd2seu1(rank,comm,d1,nx,ny,nz,nwd2seu,seu2nwd)

      call seu2nwd1(rank,comm,mu,nx,ny,nz,seu2nwd,nwd2seu)
      call nwd2seu1(rank,comm,mu,nx,ny,nz,nwd2seu,seu2nwd)

      call seu2nwd1(rank,comm,lam,nx,ny,nz,seu2nwd,nwd2seu)
      call nwd2seu1(rank,comm,lam,nx,ny,nz,nwd2seu,seu2nwd)

      if (nve==1) then
        call seu2nwd1(rank,comm,qp,nx,ny,nz,seu2nwd,nwd2seu)
        call nwd2seu1(rank,comm,qp,nx,ny,nz,nwd2seu,seu2nwd)

        call seu2nwd1(rank,comm,qs,nx,ny,nz,seu2nwd,nwd2seu)
        call nwd2seu1(rank,comm,qs,nx,ny,nz,nwd2seu,seu2nwd)
      end if
      


C     EDGE COLUMNS

      call sw2ne1(rank,comm,d1,nx,ny,nz,sw2ne,ne2sw)
      call ne2sw1(rank,comm,d1,nx,ny,nz,ne2sw,sw2ne)

      call sw2ne1(rank,comm,mu,nx,ny,nz,sw2ne,ne2sw)
      call ne2sw1(rank,comm,mu,nx,ny,nz,ne2sw,sw2ne)

      call sw2ne1(rank,comm,lam,nx,ny,nz,sw2ne,ne2sw)
      call ne2sw1(rank,comm,lam,nx,ny,nz,ne2sw,sw2ne)

      if (nve==1) then
        call sw2ne1(rank,comm,qp,nx,ny,nz,sw2ne,ne2sw)
        call ne2sw1(rank,comm,qp,nx,ny,nz,ne2sw,sw2ne)

        call sw2ne1(rank,comm,qs,nx,ny,nz,sw2ne,ne2sw)
        call ne2sw1(rank,comm,qs,nx,ny,nz,ne2sw,sw2ne)
      end if

      call se2nw1(rank,comm,d1,nx,ny,nz,se2nw,nw2se)
      call nw2se1(rank,comm,d1,nx,ny,nz,nw2se,se2nw)

      call se2nw1(rank,comm,mu,nx,ny,nz,se2nw,nw2se)
      call nw2se1(rank,comm,mu,nx,ny,nz,nw2se,se2nw)

      call se2nw1(rank,comm,lam,nx,ny,nz,se2nw,nw2se)
      call nw2se1(rank,comm,lam,nx,ny,nz,nw2se,se2nw)

      if (nve==1) then
        call se2nw1(rank,comm,qp,nx,ny,nz,se2nw,nw2se)
        call nw2se1(rank,comm,qp,nx,ny,nz,nw2se,se2nw)

        call se2nw1(rank,comm,qs,nx,ny,nz,se2nw,nw2se)
        call nw2se1(rank,comm,qs,nx,ny,nz,nw2se,se2nw)
      end if

      call su2nd1(rank,comm,d1,nx,ny,nz,su2nd,nd2su)
      call nd2su1(rank,comm,d1,nx,ny,nz,nd2su,su2nd)

      call su2nd1(rank,comm,mu,nx,ny,nz,su2nd,nd2su)
      call nd2su1(rank,comm,mu,nx,ny,nz,nd2su,su2nd)

      call su2nd1(rank,comm,lam,nx,ny,nz,su2nd,nd2su)
      call nd2su1(rank,comm,lam,nx,ny,nz,nd2su,su2nd)

      if (nve==1) then
        call su2nd1(rank,comm,qp,nx,ny,nz,su2nd,nd2su)
        call nd2su1(rank,comm,qp,nx,ny,nz,nd2su,su2nd)

        call su2nd1(rank,comm,qs,nx,ny,nz,su2nd,nd2su)
        call nd2su1(rank,comm,qs,nx,ny,nz,nd2su,su2nd)
      end if

      call sd2nu1(rank,comm,d1,nx,ny,nz,sd2nu,nu2sd)
      call nu2sd1(rank,comm,d1,nx,ny,nz,nu2sd,sd2nu)

      call sd2nu1(rank,comm,mu,nx,ny,nz,sd2nu,nu2sd)
      call nu2sd1(rank,comm,mu,nx,ny,nz,nu2sd,sd2nu)

      call sd2nu1(rank,comm,lam,nx,ny,nz,sd2nu,nu2sd)
      call nu2sd1(rank,comm,lam,nx,ny,nz,nu2sd,sd2nu)

      if (nve==1) then
        call sd2nu1(rank,comm,qp,nx,ny,nz,sd2nu,nu2sd)
        call nu2sd1(rank,comm,qp,nx,ny,nz,nu2sd,sd2nu)

        call sd2nu1(rank,comm,qs,nx,ny,nz,sd2nu,nu2sd)
        call nu2sd1(rank,comm,qs,nx,ny,nz,nu2sd,sd2nu)
      end if

      call wu2ed1(rank,comm,d1,nx,ny,nz,wu2ed,ed2wu)
      call ed2wu1(rank,comm,d1,nx,ny,nz,ed2wu,wu2ed)

      call wu2ed1(rank,comm,mu,nx,ny,nz,wu2ed,ed2wu)
      call ed2wu1(rank,comm,mu,nx,ny,nz,ed2wu,wu2ed)

      call wu2ed1(rank,comm,lam,nx,ny,nz,wu2ed,ed2wu)
      call ed2wu1(rank,comm,lam,nx,ny,nz,ed2wu,wu2ed)
 
      if (nve==1) then
        call wu2ed1(rank,comm,qp,nx,ny,nz,wu2ed,ed2wu)
        call ed2wu1(rank,comm,qp,nx,ny,nz,ed2wu,wu2ed)

        call wu2ed1(rank,comm,qs,nx,ny,nz,wu2ed,ed2wu)
        call ed2wu1(rank,comm,qs,nx,ny,nz,ed2wu,wu2ed)
      end if

      call wd2eu1(rank,comm,d1,nx,ny,nz,wd2eu,eu2wd)
      call eu2wd1(rank,comm,d1,nx,ny,nz,eu2wd,wd2eu)

      call wd2eu1(rank,comm,mu,nx,ny,nz,wd2eu,eu2wd)
      call eu2wd1(rank,comm,mu,nx,ny,nz,eu2wd,wd2eu)

      call wd2eu1(rank,comm,lam,nx,ny,nz,wd2eu,eu2wd)
      call eu2wd1(rank,comm,lam,nx,ny,nz,eu2wd,wd2eu)

      if (nve==1) then
        call wd2eu1(rank,comm,qp,nx,ny,nz,wd2eu,eu2wd)
        call eu2wd1(rank,comm,qp,nx,ny,nz,eu2wd,wd2eu)

        call wd2eu1(rank,comm,qs,nx,ny,nz,wd2eu,eu2wd)
        call eu2wd1(rank,comm,qs,nx,ny,nz,eu2wd,wd2eu)
      end if
      

 
      return
      end



      subroutine velswp(rank,comm,nx,ny,nz,ton,froms,tos,fromn,
     +                  toe,fromw,tow,frome,tou,fromd,tod,fromu)

CCC   EXCHANGES VELOCITIES BETWEEN PROCS

      use parstat

      integer :: rank,comm
      integer :: nx,ny,nz

      integer :: fromw,toe,frome,tow
      integer :: froms,ton,fromn,tos
      integer :: fromd,tou,fromu,tod

C     NORTH AND SOUTH

      call s2n(rank,comm,u1,nx,ny,nz,ton,froms,1,1,2)
      call n2s(rank,comm,u1,nx,ny,nz,tos,fromn,1,2,2)

      call s2n(rank,comm,v1,nx,ny,nz,ton,froms,2,2,2)
      call n2s(rank,comm,v1,nx,ny,nz,tos,fromn,2,1,2)

      call s2n(rank,comm,w1,nx,ny,nz,ton,froms,3,1,2)
      call n2s(rank,comm,w1,nx,ny,nz,tos,fromn,3,2,2)

C     EAST AND WEST

      call w2e(rank,comm,u1,nx,ny,nz,toe,fromw,1,1,2)
      call e2w(rank,comm,u1,nx,ny,nz,tow,frome,1,2,2)

      call w2e(rank,comm,v1,nx,ny,nz,toe,fromw,2,2,2)
      call e2w(rank,comm,v1,nx,ny,nz,tow,frome,2,1,2)

      call w2e(rank,comm,w1,nx,ny,nz,toe,fromw,3,2,2)
      call e2w(rank,comm,w1,nx,ny,nz,tow,frome,3,1,2)

C     UP AND DOWN

      call d2u(rank,comm,u1,nx,ny,nz,tou,fromd,1,1,2)
      call u2d(rank,comm,u1,nx,ny,nz,tod,fromu,1,2,2)

      call d2u(rank,comm,v1,nx,ny,nz,tou,fromd,2,1,2)
      call u2d(rank,comm,v1,nx,ny,nz,tod,fromu,2,2,2)

      call d2u(rank,comm,w1,nx,ny,nz,tou,fromd,3,2,2)
      call u2d(rank,comm,w1,nx,ny,nz,tod,fromu,3,1,2)

      return
      end


      subroutine fvelswpxy(rank,comm,nx,ny,nz,ton,froms,tow,frome)
      use parstat
      integer :: rank,comm
      integer :: nx,ny,nz
      integer :: frome,tow,froms,ton,fromn

      call fs2n(rank,comm,v1,nx,ny,nz,ton,froms,2)
      call fe2w(rank,comm,u1,nx,ny,nz,tow,frome,1)
      return
      end

      subroutine fs2n(rank,comm,b,nx,ny,nz,dst,src,ident)
      include "mpif.h"

      INTEGER :: requests(2), num_of_requests
      INTEGER :: mpistatus(MPI_STATUS_SIZE,2)
      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src,ident
      integer :: npts
      real, dimension(-1:nx+2,-1:ny+2,-1:nz+2) :: b
      real, dimension(1:nx) :: bto
      real, dimension(1:nx) :: bfrom

      npts = nx
      num_of_requests=0

      if(dst /= MPI_PROC_NULL) then
        bto = b(1:nx,ny,nz+1)
        num_of_requests=num_of_requests+1
        call MPI_ISEND(bto,npts,MPI_REAL,dst,11000+ident,comm,requests(num_of_requests),err)
      end if

      if(src /= MPI_PROC_NULL) then
        num_of_requests=num_of_requests+1
        call MPI_IRECV(bfrom,npts,MPI_REAL,src,11000+ident,comm,requests(num_of_requests),err)
      end if

      CALL MPI_WAITALL(num_of_requests,requests,mpistatus,err)
      if(src /= MPI_PROC_NULL) then
        b(1:nx,0,nz+1) = bfrom
      end if

      return
      end

      subroutine fe2w(rank,comm,b,nx,ny,nz,dst,src,ident)
      include "mpif.h"

      INTEGER :: requests(2), num_of_requests
      INTEGER :: mpistatus(MPI_STATUS_SIZE,2)
      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src,ident
      integer :: npts
      real, dimension(-1:nx+2,-1:ny+2,-1:nz+2) :: b
      real, dimension(1:ny) :: bto
      real, dimension(1:ny) :: bfrom

      npts = ny
      num_of_requests=0

      if(dst /= MPI_PROC_NULL) then
        bto = b(1,1:ny,nz+1)
        num_of_requests=num_of_requests+1
        call MPI_ISEND(bto,npts,MPI_REAL,dst,11060+ident,comm,requests(num_of_requests),err)
      end if

      if(src /= MPI_PROC_NULL) then
        num_of_requests=num_of_requests+1
        call MPI_IRECV(bfrom,npts,MPI_REAL,src,11060+ident,comm,requests(num_of_requests),err)
      end if

      CALL MPI_WAITALL(num_of_requests,requests,mpistatus,err)
      if(src /= MPI_PROC_NULL) then
        b(nx+1,1:ny,nz+1) = bfrom
      end if

      return
      end


      subroutine strswp(rank,comm,nx,ny,nz,ton,froms,tos,fromn,
     +                  toe,fromw,tow,frome,tou,fromd,tod,fromu)

CCC   EXCHANGE STRESSES BETWEEN PROCS

      use parstat

      integer :: rank,comm
      integer :: nx,ny,nz

      integer :: fromw,toe,frome,tow
      integer :: froms,ton,fromn,tos
      integer :: fromd,tou,fromu,tod

C     NORTH AND SOUTH

      call s2n(rank,comm,yy,nx,ny,nz,ton,froms,5,1,2)
      call n2s(rank,comm,yy,nx,ny,nz,tos,fromn,5,2,2)

      call s2n(rank,comm,yz,nx,ny,nz,ton,froms,8,2,2)
      call n2s(rank,comm,yz,nx,ny,nz,tos,fromn,8,1,2)

      call s2n(rank,comm,xy,nx,ny,nz,ton,froms,9,2,2)
      call n2s(rank,comm,xy,nx,ny,nz,tos,fromn,9,1,2)


C     EAST AND WEST

      call w2e(rank,comm,xx,nx,ny,nz,toe,fromw,4,2,2)
      call e2w(rank,comm,xx,nx,ny,nz,tow,frome,4,1,2)

      call w2e(rank,comm,xz,nx,ny,nz,toe,fromw,7,1,2)
      call e2w(rank,comm,xz,nx,ny,nz,tow,frome,7,2,2)

      call w2e(rank,comm,xy,nx,ny,nz,toe,fromw,9,1,2)
      call e2w(rank,comm,xy,nx,ny,nz,tow,frome,9,2,2)


C     UP AND DOWN

      call d2u(rank,comm,zz,nx,ny,nz,tou,fromd,6,1,2)
      call u2d(rank,comm,zz,nx,ny,nz,tod,fromu,6,2,2)

      call d2u(rank,comm,xz,nx,ny,nz,tou,fromd,7,2,2)
      call u2d(rank,comm,xz,nx,ny,nz,tod,fromu,7,1,2)

      call d2u(rank,comm,yz,nx,ny,nz,tou,fromd,8,2,2)
      call u2d(rank,comm,yz,nx,ny,nz,tod,fromu,8,1,2)

      return
      end



CCC   START SWAP ROUTINES

      subroutine se2nw1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts
      
      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real, dimension(1:nz) :: bto
      real, dimension(1:nz) :: bfrom

      npts = nz

      if(dst /= MPI_PROC_NULL) then
        bto = b(1,ny,1:nz)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(nx+1,0,1:nz) = bfrom
      end if

      return
      end



      subroutine nw2se1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real, dimension(1:nz) :: bto
      real, dimension(1:nz) :: bfrom

      npts = nz

      if(dst /= MPI_PROC_NULL) then
        bto = b(nx,1,1:nz)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(0,ny+1,1:nz) = bfrom
      end if

      return
      end



      subroutine sw2ne1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real, dimension(1:nz) :: bto
      real, dimension(1:nz) :: bfrom

      npts = nz

      if(dst /= MPI_PROC_NULL) then
        bto = b(nx,ny,1:nz)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(0,0,1:nz) = bfrom
      end if

      return
      end



      subroutine ne2sw1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real, dimension(1:nz) :: bto
      real, dimension(1:nz) :: bfrom

      npts = nz 

      if(dst /= MPI_PROC_NULL) then
        bto = b(1,1,1:nz)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(nx+1,ny+1,1:nz) = bfrom
      end if

      return
      end



      subroutine su2nd1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real, dimension(1:nx) :: bto
      real, dimension(1:nx) :: bfrom

      npts = nx

      if(dst /= MPI_PROC_NULL) then
        bto = b(1:nx,ny,1)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(1:nx,0,nz+1) = bfrom
      end if

      return
      end



      subroutine nd2su1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real, dimension(1:nx) :: bto
      real, dimension(1:nx) :: bfrom

      npts = nx

      if(dst /= MPI_PROC_NULL) then
        bto = b(1:nx,1,nz)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(1:nx,ny+1,0) = bfrom
      end if

      return
      end



      subroutine sd2nu1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real, dimension(1:nx) :: bto
      real, dimension(1:nx) :: bfrom

      npts = nx

      if(dst /= MPI_PROC_NULL) then
        bto = b(1:nx,ny,nz)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(1:nx,0,0) = bfrom
      end if

      return
      end



      subroutine nu2sd1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real, dimension(1:nx) :: bto
      real, dimension(1:nx) :: bfrom

      npts = nx

      if(dst /= MPI_PROC_NULL) then
        bto = b(1:nx,1,1)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(1:nx,ny+1,nz+1) = bfrom
      end if

      return
      end



      subroutine ed2wu1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real, dimension(1:ny) :: bto
      real, dimension(1:ny) :: bfrom

      npts = ny

      if(dst /= MPI_PROC_NULL) then
        bto = b(1,1:ny,nz)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(nx+1,1:ny,0) = bfrom
      end if

      return
      end



      subroutine wu2ed1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real, dimension(1:ny) :: bto
      real, dimension(1:ny) :: bfrom

      npts = ny

      if(dst /= MPI_PROC_NULL) then
        bto = b(nx,1:ny,1)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(0,1:ny,nz+1) = bfrom
      end if

      return
      end



      subroutine wd2eu1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real, dimension(1:ny) :: bto
      real, dimension(1:ny) :: bfrom

      npts = ny

      if(dst /= MPI_PROC_NULL) then
        bto = b(nx,1:ny,nz)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(0,1:ny,0) = bfrom
      end if

      return
      end



      subroutine eu2wd1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real, dimension(1:ny) :: bto
      real, dimension(1:ny) :: bfrom

      npts = ny

      if(dst /= MPI_PROC_NULL) then
        bto = b(1,1:ny,1)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(nx+1,1:ny,nz+1) = bfrom
      end if

      return
      end




      subroutine sed2nwu1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real :: bto
      real :: bfrom

      npts = 1

      if(dst /= MPI_PROC_NULL) then
        bto = b(1,ny,nz)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(nx+1,0,0) = bfrom
      end if

      return
      end



      subroutine nwu2sed1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real :: bto
      real :: bfrom

      npts = 1

      if(dst /= MPI_PROC_NULL) then
        bto = b(nx,1,1)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(0,ny+1,nz+1) = bfrom
      end if

      return
      end



      subroutine seu2nwd1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real :: bto
      real :: bfrom

      npts = 1

      if(dst /= MPI_PROC_NULL) then
        bto = b(1,ny,1)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(nx+1,0,nz+1) = bfrom
      end if

      return
      end



      subroutine nwd2seu1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real :: bto
      real :: bfrom

      npts = 1

      if(dst /= MPI_PROC_NULL) then
        bto = b(nx,1,nz)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(0,ny+1,0) = bfrom
      end if

      return
      end



      subroutine swd2neu1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real :: bto
      real :: bfrom

      npts = 1

      if(dst /= MPI_PROC_NULL) then
        bto = b(nx,ny,nz)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(0,0,0) = bfrom
      end if

      return
      end

 

      subroutine neu2swd1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real :: bto
      real :: bfrom

      npts = 1

      if(dst /= MPI_PROC_NULL) then
        bto = b(1,1,1) 
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(nx+1,ny+1,nz+1) = bfrom
      end if

      return
      end



      subroutine swu2ned1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real :: bto
      real :: bfrom

      npts = 1

      if(dst /= MPI_PROC_NULL) then
        bto = b(nx,ny,1)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(0,0,nz+1) = bfrom
      end if

      return
      end



      subroutine ned2swu1(rank,comm,b,nx,ny,nz,dst,src)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"

      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src
      integer :: npts

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(0:nx+1,0:ny+1,0:nz+1) :: b

      real :: bto
      real :: bfrom

      npts = 1

      if(dst /= MPI_PROC_NULL) then
        bto = b(1,1,nz)
        call MPI_SEND(bto,npts,MPI_REAL,dst,rank,comm,err)
      end if

      if(src /= MPI_PROC_NULL) then
        call MPI_RECV(bfrom,npts,MPI_REAL,src,MPI_ANY_TAG,comm,status,err)
        b(nx+1,ny+1,0) = bfrom
      end if

      return
      end



      subroutine s2n(rank,comm,b,nx,ny,nz,dst,src,ident,layer,width)

CCC   TRANSFER X-Z PLANES FROM S TO N PROCS

      include "mpif.h"
      INTEGER :: requests(2), num_of_requests
      INTEGER :: mpistatus(MPI_STATUS_SIZE,2)
      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src,ident
      integer :: npts,width,layer

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(1-width:nx+width,1-width:ny+width,1-width:nz+width) :: b

      real, dimension(1:nx,ny-layer+1:ny,1:nz) :: bto
      real, dimension(1:nx,1-layer:0,1:nz) :: bfrom

      npts = layer*nx*nz
      num_of_requests=0

      if(dst /= MPI_PROC_NULL) then
        bto = b(1:nx,ny-layer+1:ny,1:nz)
        num_of_requests=num_of_requests+1
        call MPI_ISEND(bto,npts,MPI_REAL,dst,10000+ident,comm,requests(num_of_requests),err)
      end if

      if(src /= MPI_PROC_NULL) then
        num_of_requests=num_of_requests+1
        call MPI_IRECV(bfrom,npts,MPI_REAL,src,10000+ident,comm,requests(num_of_requests),err)
      end if
      
      CALL MPI_WAITALL(num_of_requests,requests,mpistatus,err)
      if(src /= MPI_PROC_NULL) then
        b(1:nx,1-layer:0,1:nz) = bfrom
      end if

      return
      end



      subroutine n2s(rank,comm,b,nx,ny,nz,dst,src,ident,layer,width)

CCC   TRANSFER X-Z PLANES FROM N TO S PROCS

      include "mpif.h"

      INTEGER :: requests(2), num_of_requests
      INTEGER :: mpistatus(MPI_STATUS_SIZE,2)
      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src,ident
      integer :: npts,layer,width

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(1-width:nx+width,1-width:ny+width,1-width:nz+width) :: b

      real, dimension(1:nx,1:layer,1:nz) :: bto
      real, dimension(1:nx,ny+1:ny+layer,1:nz) :: bfrom

      npts = layer*nx*nz
      num_of_requests=0

      if(dst /= MPI_PROC_NULL) then
        bto = b(1:nx,1:layer,1:nz)
        num_of_requests=num_of_requests+1
        call MPI_ISEND(bto,npts,MPI_REAL,dst,10020+ident,comm,requests(num_of_requests),err)
      end if

      if(src /= MPI_PROC_NULL) then
        num_of_requests=num_of_requests+1
        call MPI_IRECV(bfrom,npts,MPI_REAL,src,10020+ident,comm,requests(num_of_requests),err)
      end if

      CALL MPI_WAITALL(num_of_requests,requests,mpistatus,err)
      if(src /= MPI_PROC_NULL) then
        b(1:nx,ny+1:ny+layer,1:nz) = bfrom
      end if

      return
      end



      subroutine w2e(rank,comm,b,nx,ny,nz,dst,src,ident,layer,width)

CCC   TRANSFER Y-Z PLANES FROM W TO E PROCS

      include "mpif.h"

      INTEGER :: requests(2), num_of_requests
      INTEGER :: mpistatus(MPI_STATUS_SIZE,2)
      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src,ident
      integer :: npts,layer,width

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(1-width:nx+width,1-width:ny+width,1-width:nz+width) :: b

      real, dimension(nx-layer+1:nx,1:ny,1:nz) :: bto
      real, dimension(1-layer:0,1:ny,1:nz) :: bfrom

      npts = layer*ny*nz
      num_of_requests=0

      if(dst /= MPI_PROC_NULL) then
        bto = b(nx-layer+1:nx,1:ny,1:nz)
        num_of_requests=num_of_requests+1
        call MPI_ISEND(bto,npts,MPI_REAL,dst,10040+ident,comm,requests(num_of_requests),err)
      end if

      if(src /= MPI_PROC_NULL) then
        num_of_requests=num_of_requests+1
        call MPI_IRECV(bfrom,npts,MPI_REAL,src,10040+ident,comm,requests(num_of_requests),err)
      end if

      CALL MPI_WAITALL(num_of_requests,requests,mpistatus,err)
      if(src /= MPI_PROC_NULL) then
        b(1-layer:0,1:ny,1:nz) = bfrom
      end if

      return
      end



      subroutine e2w(rank,comm,b,nx,ny,nz,dst,src,ident,layer,width)

CCC   TRANSFER Y-Z PLANES FROM E TO W PROCS

      include "mpif.h"

      INTEGER :: requests(2), num_of_requests
      INTEGER :: mpistatus(MPI_STATUS_SIZE,2)
      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src,ident
      integer :: npts,layer,width

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(1-width:nx+width,1-width:ny+width,1-width:nz+width) :: b

      real, dimension(1:layer,1:ny,1:nz) :: bto
      real, dimension(nx+1:nx+layer,1:ny,1:nz) :: bfrom

      npts = layer*ny*nz
      num_of_requests=0

      if(dst /= MPI_PROC_NULL) then
        bto = b(1:layer,1:ny,1:nz)
        num_of_requests=num_of_requests+1
        call MPI_ISEND(bto,npts,MPI_REAL,dst,10060+ident,comm,requests(num_of_requests),err)
      end if

      if(src /= MPI_PROC_NULL) then
        num_of_requests=num_of_requests+1
        call MPI_IRECV(bfrom,npts,MPI_REAL,src,10060+ident,comm,requests(num_of_requests),err)
      end if

      CALL MPI_WAITALL(num_of_requests,requests,mpistatus,err)
      if(src /= MPI_PROC_NULL) then
        b(nx+1:nx+layer,1:ny,1:nz) = bfrom
      end if

      return
      end



      subroutine d2u(rank,comm,b,nx,ny,nz,dst,src,ident,layer,width)

CCC   TRANSFER X-Y PLANES FROM D TO U PROCS

      include "mpif.h"

      INTEGER :: requests(2), num_of_requests
      INTEGER :: mpistatus(MPI_STATUS_SIZE,2)
      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src,ident
      integer :: npts,layer,width

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(1-width:nx+width,1-width:ny+width,1-width:nz+width) :: b

      real, dimension(1:nx,1:ny,nz-layer+1:nz) :: bto
      real, dimension(1:nx,1:ny,1-layer:0) :: bfrom

      npts = layer*nx*ny
      num_of_requests=0

      if(dst /= MPI_PROC_NULL) then
        bto = b(1:nx,1:ny,nz-layer+1:nz)
        num_of_requests=num_of_requests+1
        call MPI_ISEND(bto,npts,MPI_REAL,dst,10080+ident,comm,requests(num_of_requests),err)
      end if

      if(src /= MPI_PROC_NULL) then
        num_of_requests=num_of_requests+1
        call MPI_IRECV(bfrom,npts,MPI_REAL,src,10080+ident,comm,requests(num_of_requests),err)
      end if

      CALL MPI_WAITALL(num_of_requests,requests,mpistatus,err)
      if(src /= MPI_PROC_NULL) then
        b(1:nx,1:ny,1-layer:0) = bfrom
      end if

      return
      end



      subroutine u2d(rank,comm,b,nx,ny,nz,dst,src,ident,layer,width)

CCC   TRANSFER X-Y PLANES FROM U TO D PROCS

      include "mpif.h"

      INTEGER :: requests(2), num_of_requests
      INTEGER :: mpistatus(MPI_STATUS_SIZE,2)
      integer :: rank,comm,err
      integer :: nx,ny,nz
      integer :: dst,src,ident
      integer :: npts,layer,width

      integer, dimension(MPI_STATUS_SIZE) :: status

      real, dimension(1-width:nx+width,1-width:ny+width,1-width:nz+width) :: b

      real, dimension(1:nx,1:ny,1:layer) :: bto
      real, dimension(1:nx,1:ny,nz+1:nz+layer) :: bfrom

      npts = layer*nx*ny
      num_of_requests=0

      if(dst /= MPI_PROC_NULL) then
        bto = b(1:nx,1:ny,1:layer)
        num_of_requests=num_of_requests+1
        call MPI_ISEND(bto,npts,MPI_REAL,dst,10100+ident,comm,requests(num_of_requests),err)
      end if

      if(src /= MPI_PROC_NULL) then
        num_of_requests=num_of_requests+1
        call MPI_IRECV(bfrom,npts,MPI_REAL,src,10100+ident,comm,requests(num_of_requests),err)
      end if

      CALL MPI_WAITALL(num_of_requests,requests,mpistatus,err)
      if(src /= MPI_PROC_NULL) then
        b(1:nx,1:ny,nz+1:nz+layer) = bfrom
      end if

      return
      end


      subroutine velswp_s(rank,comm,ton,tos,toe,tow,tou,tod,
     +                  send_d2u_1,send_u2d_1,send_d2u_2,send_u2d_2,
     +                  send_s2n_1,send_n2s_1,send_s2n_2,send_n2s_2,
     +                  send_w2e_1,send_e2w_1,send_w2e_2,send_e2w_2)

      use parstat
      integer :: rank,comm,err
      integer :: toe,tow,ton,tos,tou,tod
      integer :: send_d2u_1,send_u2d_1,send_d2u_2,send_u2d_2
      integer :: send_s2n_1,send_n2s_1,send_s2n_2,send_n2s_2
      integer :: send_w2e_1,send_e2w_1,send_w2e_2,send_e2w_2

      if(ton /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(u1(-1,-1,-1),1,send_s2n_1,ton,10001,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(v1(-1,-1,-1),1,send_s2n_2,ton,10002,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(w1(-1,-1,-1),1,send_s2n_1,ton,10003,comm,requests_v(num_of_requests_v),err)
      end if

      if(tos /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(u1(-1,-1,-1),1,send_n2s_2,tos,10021,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(v1(-1,-1,-1),1,send_n2s_1,tos,10022,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(w1(-1,-1,-1),1,send_n2s_2,tos,10023,comm,requests_v(num_of_requests_v),err)
      end if

      if(toe /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(u1(-1,-1,-1),1,send_w2e_1,toe,10041,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(v1(-1,-1,-1),1,send_w2e_2,toe,10042,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(w1(-1,-1,-1),1,send_w2e_2,toe,10043,comm,requests_v(num_of_requests_v),err)
      end if

      if(tow /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(u1(-1,-1,-1),1,send_e2w_2,tow,10061,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(v1(-1,-1,-1),1,send_e2w_1,tow,10062,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(w1(-1,-1,-1),1,send_e2w_1,tow,10063,comm,requests_v(num_of_requests_v),err)
      end if

      if(tou /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(u1(-1,-1,-1),1,send_d2u_1,tou,10081,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(v1(-1,-1,-1),1,send_d2u_1,tou,10082,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(w1(-1,-1,-1),1,send_d2u_2,tou,10083,comm,requests_v(num_of_requests_v),err)
      end if

      if(tod /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(u1(-1,-1,-1),1,send_u2d_2,tod,10101,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(v1(-1,-1,-1),1,send_u2d_2,tod,10102,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(w1(-1,-1,-1),1,send_u2d_1,tod,10103,comm,requests_v(num_of_requests_v),err)
      end if

      return
      end

      subroutine velswp_r(rank,comm,froms,fromn,fromw,frome,fromd,fromu,
     +                  recv_d2u_1,recv_u2d_1,recv_d2u_2,recv_u2d_2,
     +                  recv_s2n_1,recv_n2s_1,recv_s2n_2,recv_n2s_2,
     +                  recv_w2e_1,recv_e2w_1,recv_w2e_2,recv_e2w_2)

      use parstat
      integer :: rank,comm,err
      integer :: fromw,frome,froms,fromn,fromd,fromu
      integer :: recv_d2u_1,recv_u2d_1,recv_d2u_2,recv_u2d_2
      integer :: recv_s2n_1,recv_n2s_1,recv_s2n_2,recv_n2s_2
      integer :: recv_w2e_1,recv_e2w_1,recv_w2e_2,recv_e2w_2

      num_of_requests_v = 0

      if(froms /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(u1(-1,-1,-1),1,recv_s2n_1,froms,10001,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(v1(-1,-1,-1),1,recv_s2n_2,froms,10002,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(w1(-1,-1,-1),1,recv_s2n_1,froms,10003,comm,requests_v(num_of_requests_v),err)
      end if

      if(fromn /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(u1(-1,-1,-1),1,recv_n2s_2,fromn,10021,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(v1(-1,-1,-1),1,recv_n2s_1,fromn,10022,comm,requests_v(num_of_requests_v),err)  
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(w1(-1,-1,-1),1,recv_n2s_2,fromn,10023,comm,requests_v(num_of_requests_v),err)  
      end if

      if(fromw /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(u1(-1,-1,-1),1,recv_w2e_1,fromw,10041,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(v1(-1,-1,-1),1,recv_w2e_2,fromw,10042,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(w1(-1,-1,-1),1,recv_w2e_2,fromw,10043,comm,requests_v(num_of_requests_v),err)
      end if

      if(frome /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(u1(-1,-1,-1),1,recv_e2w_2,frome,10061,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(v1(-1,-1,-1),1,recv_e2w_1,frome,10062,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(w1(-1,-1,-1),1,recv_e2w_1,frome,10063,comm,requests_v(num_of_requests_v),err)
      end if

      if(fromd /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(u1(-1,-1,-1),1,recv_d2u_1,fromd,10081,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(v1(-1,-1,-1),1,recv_d2u_1,fromd,10082,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(w1(-1,-1,-1),1,recv_d2u_2,fromd,10083,comm,requests_v(num_of_requests_v),err)
      end if

      if(fromu /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(u1(-1,-1,-1),1,recv_u2d_2,fromu,10101,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(v1(-1,-1,-1),1,recv_u2d_2,fromu,10102,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(w1(-1,-1,-1),1,recv_u2d_1,fromu,10103,comm,requests_v(num_of_requests_v),err)
      end if

      return
      end

      subroutine wait_all_vel()
      use parstat
      integer :: mpistatus(MPI_STATUS_SIZE,36)

      call MPI_WAITALL(num_of_requests_v,requests_v,mpistatus,err)

      return
      end

      subroutine strswp_r(rank,comm,froms,fromn,fromw,frome,fromd,fromu,
     +                  recv_d2u_1,recv_u2d_1,recv_d2u_2,recv_u2d_2,
     +                  recv_s2n_1,recv_n2s_1,recv_s2n_2,recv_n2s_2,
     +                  recv_w2e_1,recv_e2w_1,recv_w2e_2,recv_e2w_2)

      use parstat
      integer :: rank,comm,err
      integer :: fromw,frome,froms,fromn,fromd,fromu
      integer :: recv_d2u_1,recv_u2d_1,recv_d2u_2,recv_u2d_2
      integer :: recv_s2n_1,recv_n2s_1,recv_s2n_2,recv_n2s_2
      integer :: recv_w2e_1,recv_e2w_1,recv_w2e_2,recv_e2w_2

      num_of_requests_v = 0

      if(froms /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(yy(-1,-1,-1),1,recv_s2n_1,froms,10005,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(yz(-1,-1,-1),1,recv_s2n_2,froms,10008,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(xy(-1,-1,-1),1,recv_s2n_2,froms,10009,comm,requests_v(num_of_requests_v),err)
      end if

      if(fromn /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(yy(-1,-1,-1),1,recv_n2s_2,fromn,10025,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(yz(-1,-1,-1),1,recv_n2s_1,fromn,10028,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(xy(-1,-1,-1),1,recv_n2s_1,fromn,10029,comm,requests_v(num_of_requests_v),err)
      end if

      if(fromw /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(xx(-1,-1,-1),1,recv_w2e_2,fromw,10044,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(xz(-1,-1,-1),1,recv_w2e_1,fromw,10047,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(xy(-1,-1,-1),1,recv_w2e_1,fromw,10049,comm,requests_v(num_of_requests_v),err)
      end if

      if(frome /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(xx(-1,-1,-1),1,recv_e2w_1,frome,10064,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(xz(-1,-1,-1),1,recv_e2w_2,frome,10067,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(xy(-1,-1,-1),1,recv_e2w_2,frome,10069,comm,requests_v(num_of_requests_v),err)
      end if

      if(fromd /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(zz(-1,-1,-1),1,recv_d2u_1,fromd,10086,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(xz(-1,-1,-1),1,recv_d2u_2,fromd,10087,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(yz(-1,-1,-1),1,recv_d2u_2,fromd,10088,comm,requests_v(num_of_requests_v),err)
      end if

      if(fromu /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(zz(-1,-1,-1),1,recv_u2d_2,fromu,10106,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(xz(-1,-1,-1),1,recv_u2d_1,fromu,10107,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_IRECV(yz(-1,-1,-1),1,recv_u2d_1,fromu,10108,comm,requests_v(num_of_requests_v),err)
      end if

      return
      end

      subroutine strswp_s(rank,comm,ton,tos,toe,tow,tou,tod,
     +                  send_d2u_1,send_u2d_1,send_d2u_2,send_u2d_2,
     +                  send_s2n_1,send_n2s_1,send_s2n_2,send_n2s_2,
     +                  send_w2e_1,send_e2w_1,send_w2e_2,send_e2w_2)

      use parstat
      integer :: rank,comm,err
      integer :: toe,tow,ton,tos,tou,tod
      integer :: send_d2u_1,send_u2d_1,send_d2u_2,send_u2d_2
      integer :: send_s2n_1,send_n2s_1,send_s2n_2,send_n2s_2
      integer :: send_w2e_1,send_e2w_1,send_w2e_2,send_e2w_2

      if(ton /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(yy(-1,-1,-1),1,send_s2n_1,ton,10005,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(yz(-1,-1,-1),1,send_s2n_2,ton,10008,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(xy(-1,-1,-1),1,send_s2n_2,ton,10009,comm,requests_v(num_of_requests_v),err)
      end if

      if(tos /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(yy(-1,-1,-1),1,send_n2s_2,tos,10025,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(yz(-1,-1,-1),1,send_n2s_1,tos,10028,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(xy(-1,-1,-1),1,send_n2s_1,tos,10029,comm,requests_v(num_of_requests_v),err)
      end if

      if(toe /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(xx(-1,-1,-1),1,send_w2e_2,toe,10044,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(xz(-1,-1,-1),1,send_w2e_1,toe,10047,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(xy(-1,-1,-1),1,send_w2e_1,toe,10049,comm,requests_v(num_of_requests_v),err)
      end if

      if(tow /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(xx(-1,-1,-1),1,send_e2w_1,tow,10064,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(xz(-1,-1,-1),1,send_e2w_2,tow,10067,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(xy(-1,-1,-1),1,send_e2w_2,tow,10069,comm,requests_v(num_of_requests_v),err)
      end if

      if(tou /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(zz(-1,-1,-1),1,send_d2u_1,tou,10086,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(xz(-1,-1,-1),1,send_d2u_2,tou,10087,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(yz(-1,-1,-1),1,send_d2u_2,tou,10088,comm,requests_v(num_of_requests_v),err)
      end if

      if(tod /= MPI_PROC_NULL) then
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(zz(-1,-1,-1),1,send_u2d_2,tod,10106,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(xz(-1,-1,-1),1,send_u2d_1,tod,10107,comm,requests_v(num_of_requests_v),err)
        num_of_requests_v=num_of_requests_v+1
        call MPI_ISEND(yz(-1,-1,-1),1,send_u2d_1,tod,10108,comm,requests_v(num_of_requests_v),err)
      end if

      return
      end

