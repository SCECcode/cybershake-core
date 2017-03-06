CCC   'sgsndyna.f' COMPUTES VELOCITIES AND STRESSES OF SGSN DYNAMIC FAULT MODEL
CCC    of Dalguer and Day (2006)
CCC   
CCC  The fault is positioned at the XZ plane
CCC
CCC WARNING!!, The fault plane should be located at a distance of at least two grids 
CCC from the borders of the processor (y directions)
CCC           
CCC

      subroutine dynavel(ix_d,fx_d,y_d,iz_d,fz_d,
     +                   nxt,nzt,dh,dt,fsr)

      use parstat

      integer :: ix_d,fx_d,y_d,iz_d,fz_d,ix,iz,fx,fz,fsr      
      fx=0
      fz=0      
      ix=0
      iz=0            
      if(nxt.gt.fx_d) fx=1
      if(nzt.gt.fz_d) fz=1      

      if(ix_d.gt.1) ix=1
      if(iz_d.gt.1) iz=1         
      
c 2th order finite-difference 
       call ux2th(ix_d,fx_d,y_d+1,y_d+1,iz_d,fz_d,dh,dt,y_d)
       call wz2th(ix_d,fx_d-fx,y_d+1,y_d+1,iz_d,fz_d-fz,dh,dt,y_d)       
       call ux2th(ix_d,fx_d,y_d-1,y_d-1,iz_d,fz_d,dh,dt,y_d)       
       call wz2th(ix_d,fx_d-fx,y_d-1,y_d-1,iz_d,fz_d-fz,dh,dt,y_d)   
       call vy2th(ix_d,fx_d-fx,y_d-1,y_d,iz_d,fz_d,dh,dt,y_d)
c west border        
       if(ix.eq.1) then
         call uxedge(ix_d-3,ix_d-1,y_d,y_d,iz_d-iz,fz_d+fz,dh,dt)
         call wzedge(ix_d-3,ix_d-1,y_d,y_d,iz_d-iz,fz_d+fz,dh,dt)       
       endif
c east border
       if(fx.eq.1) then
         call uxedge(fx_d+1,fx_d+3,y_d,y_d,iz_d-iz,fz_d+fz,dh,dt)
         call wzedge(fx_d,fx_d+3,y_d,y_d,iz_d-iz,fz_d+fz,dh,dt)       
       endif
c base border
       if(iz.eq.1) then
         call uxedge(ix_d,fx_d,y_d,y_d,iz_d-3,iz_d-1,dh,dt)
         call wzedge(ix_d,fx_d-fx,y_d,y_d,iz_d-3,iz_d-1,dh,dt)       
       endif

c top border
       if(fz.eq.1) then
         call uxedge(ix_d,fx_d,y_d,y_d,fz_d+1,fz_d+3,dh,dt)
         call wzedge(ix_d,fx_d-fx,y_d,y_d,fz_d,fz_d+3,dh,dt)       
       endif       
       
              
c split nodes velocities

         call ux2thsplit(ix_d,fx_d,y_d,iz_d,fz_d,fsr,dh,dt)
         call wz2thsplit(ix_d,fx_d-fx,y_d,iz_d,fz_d-fz,fsr,dh,dt)         
c free surface rupture (fsr = 1)
c No free surface rupture (fsr = 0)


       
      return
      end

      subroutine dynastress(it,ix_d,fx_d,y_d,iz_d,fz_d,nxt,nzt,
     +                      fsr,dh,dt)
     
c     it     current time step in FD code (int)
      use parstat

      integer :: ix_d,fx_d,y_d,iz_d,fz_d,ix,iz,fx,fz,fsr

      fx=0
      fz=0      
      ix=0
      iz=0            
      if(nxt.gt.fx_d) fx=1
      if(nzt.gt.fz_d) fz=1      

      if(ix_d.gt.1) ix=1
      if(iz_d.gt.1) iz=1      

      
c 4th order finite-difference including split nodes 
       call xy4th_plus(ix_d,fx_d,y_d+1,y_d+1,iz_d,fz_d,dh,dt,y_d)
       
       call yz4th_plus(ix_d,fx_d-fx,y_d+1,y_d+1,iz_d,fz_d-fz,dh,dt,y_d)

c 2th order finite-difference 
       call xyz2th(ix_d,fx_d-fx,y_d+1,y_d+1,iz_d,fz_d,dh,dt,y_d)
       call xyz2th(ix_d,fx_d-fx,y_d-1,y_d-1,iz_d,fz_d,dh,dt,y_d)       

c 2th order finite-difference including split nodes
       call xy2th_plus(ix_d,fx_d,y_d,y_d,iz_d,fz_d,dh,dt,y_d)
       call xy2th_minus(ix_d,fx_d,y_d-1,y_d-1,iz_d,fz_d,dh,dt,y_d)       
       
       call yz2th_plus(ix_d,fx_d-fx,y_d,y_d,iz_d,fz_d-fz,dh,dt,y_d)
       call yz2th_minus(ix_d,fx_d-fx,y_d-1,y_d-1,iz_d,fz_d-fz,dh,dt,y_d)       

c west border        
       if(ix.eq.1) then
         call xzedge(ix_d-1,ix_d-1,y_d,y_d,iz_d,fz_d-fz)
         call xyzedge(ix_d-2,ix_d-1,y_d,y_d,iz_d,fz_d)       
       endif
c east border
       if(fx.eq.1) then
         call xzedge(fx_d+1,fx_d+1,y_d,y_d,iz_d,fz_d-fz)
         call xyzedge(fx_d,fx_d+1,y_d,y_d,iz_d,fz_d)       
       endif
c base border
       if(iz.eq.1) then
         call xyzedge(ix_d,fx_d-fx,y_d,y_d,iz_d-1,iz_d-1)
         call xzedge(ix_d,fx_d,y_d,y_d,iz_d-2,iz_d-1)       
       endif

c top border
       if(fz.eq.1) then
         call xyzedge(ix_d,fx_d-fx,y_d,y_d,fz_d+1,fz_d+1)
         call xzedge(ix_d,fx_d,y_d,y_d,fz_d,fz_d+1)       
       endif
       
c split nodes stresses 
c free surface rupture (fsr = 1)
c No free surface rupture (fsr = 0)

       call xyzsplit(ix_d,fx_d-fx,y_d,iz_d,fz_d,dh,dt)
       call xzsplit(ix_d,fx_d,y_d,iz_d,fz_d-fz,fsr,dh,dt)  

      return
      end

      subroutine dynarestrac(it,ix_d,fx_d,y_d,iz_d,fz_d,nxt,nzt,
     +                      fsr,dh,dt,eta)
C Compute restoring forces and traction on the fault     
     
c     it     current time step in FD code (int)
      use parstat
      real :: eta
      integer :: ix_d,fx_d,y_d,iz_d,fz_d,ix,iz,fx,fz,fsr
      fx=0
      fz=0      
      ix=0
      iz=0            
      if(nxt.gt.fx_d) fx=1
      if(nzt.gt.fz_d) fz=1      

      if(ix_d.gt.1) ix=1
      if(iz_d.gt.1) iz=1      

      
       call resttract_u(ix_d,fx_d,y_d,iz_d,fz_d,fsr,dh,dt,eta)
       call resttract_w(ix_d,fx_d-fx,y_d,iz_d,fz_d-fz,fsr,dh,dt,eta)               

      return
      end

       subroutine ux2th(nxb,nxe,nyb,nye,nzb,nze,dh,dt,jf)

c     2nd order finite-difference of u1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)
c     jf    fault plane location (y direction)  (integer)   (sent)

      use parstat

      dth = dt/dh
c     Find u-displacement fields at time t+1/2
      jj=jf-3
      do 50 k= nzb,nze
      do 50 j= nyb,nye
      do 50 i= nxb,nxe

      d = 0.25*(d1(i,j,k)+d1(i,j-1,k)+
     +          d1(i,j,k-1)+d1(i,j-1,k-1))

       dxx=xx(i,j,k)-xx(i-1,j,k)
       dxy=xy(i,j,k)-xy(i,j-1,k)
       dxz=xz(i,j,k)-xz(i,j,k-1)

       u_f(i,j-jj,k)=u_f(i,j-jj,k)+(dth/d)*(dxx+dxy+dxz)     
       u1(i,j,k)=u_f(i,j-jj,k)     
   50 continue

      return
      end

       subroutine uxedge(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     2nd order finite-difference of u1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)
c     jf    fault plane location (y direction)  (integer)   (sent)

      use parstat

      dth = dt/dh
c     Find u-displacement fields at time t+1/2
      do 50 k= nzb,nze
      do 50 j= nyb,nye
      do 50 i= nxb,nxe
      
       uplus(i,k)=u1(i,j,k)
     
   50 continue

      return
      end


       subroutine wzedge(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     2nd order finite-difference of u1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)
c     jf    fault plane location (y direction)  (integer)   (sent)

      use parstat

      dth = dt/dh
c     Find u-displacement fields at time t+1/2

      do 50 k= nzb,nze
      do 50 j= nyb,nye
      do 50 i= nxb,nxe

       wplus(i,k)=w1(i,j,k)            

   50 continue

      return
      end             

      
      subroutine vy2th(nxb,nxe,nyb,nye,nzb,nze,dh,dt,jf)

c     2nd order finite-difference of u1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)
c     jf    fault plane location (y direction)  (integer)   (sent)

      use parstat

      dth = dt/dh
c     Find u-displacement fields at time t+1/2
      jj=jf-3
      do 50 k= nzb,nze
      do 50 j= nyb,nye
      do 50 i= nxb,nxe

      d = 0.25*(d1(i,j,k)+d1(i+1,j,k)+
     +          d1(i,j,k-1)+d1(i+1,j,k-1))

       dyy=yy(i,j+1,k)-yy(i,j,k)
       dxy=xy(i+1,j,k)-xy(i,j,k)
       dyz=yz(i,j,k)-yz(i,j,k-1)

       v_f(i,j-jj,k)=v_f(i,j-jj,k)+(dth/d)*(dyy+dxy+dyz)     
       v1(i,j,k)=v_f(i,j-jj,k)     
   50 continue

      return
      end       

       subroutine wz2th(nxb,nxe,nyb,nye,nzb,nze,dh,dt,jf)

c     2nd order finite-difference of u1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)
c     jf    fault plane location (y direction)  (integer)   (sent)

      use parstat

      dth = dt/dh
c     Find u-displacement fields at time t+1/2
      jj=jf-3
      do 50 k= nzb,nze
      do 50 j= nyb,nye
      do 50 i= nxb,nxe

      d = 0.25*(d1(i,j,k)+d1(i+1,j,k)+
     +          d1(i,j-1,k)+d1(i+1,j-1,k))

       dzz=zz(i,j,k+1)-zz(i,j,k)
       dxz=xz(i+1,j,k)-xz(i,j,k)
       dyz=yz(i,j,k)-yz(i,j-1,k)

       w_f(i,j-jj,k)=w_f(i,j-jj,k)+(dth/d)*(dzz+dxz+dyz)     
       w1(i,j,k)=w_f(i,j-jj,k)     
   50 continue

      return
      end       
 


       subroutine ux2thsplit(nxb,nxe,nyb,nzb,nze,fsr,dh,dt)     
       use parstat
c     Find u-velocity fields of split nodes on the fault
      integer :: fsr
      jj=nyb-3
      j=nyb
      do 50 k= nzb,nze
      do 50 i= nxb,nxe
       area=dh*dh
       tru1(i,k)=tru2(i,k)       
       temp=(tru1(i,k)-strisx(i,k))*area
c plus side of fault  
       d = 0.5*(d1(i,j,k)+d1(i,j,k-1))
       smas=d*area*dh/2.   
 
       uplus(i,k)=uplus(i,k)+dt*(Rplusu(i,k)-temp)/smas 

c minus side of fault       
       d = 0.5*(d1(i,j-1,k)+d1(i,j-1,k-1))
       smas=d*area*dh/2.
       u_f(i,j-jj,k)=u_f(i,j-jj,k)+
     +               dt*(Rminusu(i,k)+temp)/smas 

c        if(brokenu(i,k).eq.0) then
c         uplus(i,k)=(uplus(i,k)+u_f(i,j-jj,k))/2.
c         u_f(i,j-jj,k)=uplus(i,k)
c        endif

 
         uplus(i,k)=(uplus(i,k)*(brokenu(i,k)+1)+u_f(i,j-jj,k)*(1-brokenu(i,k)))/2.
         u_f(i,j-jj,k)=u_f(i,j-jj,k)*brokenu(i,k) + uplus(i,k)*(1-brokenu(i,k))        
        
        u1(i,j,k)=u_f(i,j-jj,k) 
        
   50 continue              
       return
       end

       subroutine wz2thsplit(nxb,nxe,nyb,nzb,nze,fsr,dh,dt)     
       use parstat

c     Find w-velocity fields of split nodes on the fault
      integer :: fsr

      jj=nyb-3
      j=nyb
      do 50 k= nzb,nze
      do 50 i= nxb,nxe
       area=dh*dh/(1.+fsr*int(k/nze))      
c      if(fsr.eq.1.and.k.eq.nze) area=dh*dh/2.    

c plus side of fault  
       d = 0.5*(d1(i,j,k)+d1(i+1,j,k))
       smas=d*area*dh/2.   
       trw1(i,k)=trw2(i,k)  
       temp=(trw1(i,k)-strisz(i,k))*area       
       wplus(i,k)=wplus(i,k)+dt*(Rplusw(i,k)-temp)/smas 
c
c minus side of fault       
       d = 0.5*(d1(i,j-1,k)+d1(i+1,j-1,k))
       smas=d*area*dh/2.
       w_f(i,j-jj,k)=w_f(i,j-jj,k)+
     +               dt*(Rminusw(i,k)+temp)/smas 

c        if(brokenw(i,k).eq.0) then
c         wplus(i,k)=(wplus(i,k)+w_f(i,j-jj,k))/2.
c         w_f(i,j-jj,k)=wplus(i,k)
c        endif
        
         wplus(i,k)=(wplus(i,k)*(brokenw(i,k)+1)+w_f(i,j-jj,k)*(1-brokenw(i,k)))/2.
         w_f(i,j-jj,k)=w_f(i,j-jj,k)*brokenw(i,k) + wplus(i,k)*(1-brokenw(i,k))

         w1(i,j,k)=w_f(i,j-jj,k)           
        
   50 continue     
   
       return
       end
       
       
       subroutine sourceslip(nxb,nxe,nyb,nzb,nze,nzt,fsr,dh,dt)

       use parstat
       integer :: fsr
          
       
        j=nyb

      do 50 k= nzb,nze
      do 50 i= nxb,nxe
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Compute V components
c
c---V Slip-rate magnitude and slip-rate components X and Z
c
        rateu_u(i,k) = uplus(i,k) - u1(i,j,k)
        rateu_w(i,k) = (wplus(i,k)+wplus(i-1,k) +
     &                  wplus(i,k-1)+wplus(i-1,k-1))/4.0 - 
     &                 (w1(i,j,k)+w1(i-1,j,k) +
     &                  w1(i,j,k-1)+w1(i-1,j,k-1))/4.0
        rateu(i,k) = sqrt(rateu_u(i,k)**2 + rateu_w(i,k)**2)
c
c---V Slip magnitude and slip components X and Z
c
        d_u(i,k) = d_u(i,k) + dt*rateu(i,k)
        slipu_u(i,k) = slipu_u(i,k) + dt*rateu_u(i,k)
        slipu_w(i,k) = slipu_w(i,k) + dt*rateu_w(i,k)
c
c---V Velocity field components X and Z at both sides of the fault
c
        wplusi(i,k) = (wplus(i,k)+wplus(i-1,k) +
     &                 wplus(i,k-1)+wplus(i-1,k-1))/4.0

        umin(i,k) = u1(i,j,k)
        wmin(i,k) = (w1(i,j,k)+w1(i-1,j,k) +
     &               w1(i,j,k-1)+w1(i-1,j,k-1))/4.0
c
c---V Displacement field components X and Z at both sides of the fault
c
        duplus(i,k) = duplus(i,k) + dt*uplus(i,k)
        dwplus(i,k) = dwplus(i,k) + dt*wplusi(i,k)
        dumin(i,k) = dumin(i,k) + dt*umin(i,k)
        dwmin(i,k) = dwmin(i,k) + dt*wmin(i,k)
c        
c +++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Compute W components
c

c safe slip and slip rate in variables originally set in the code
      u2(i,j,k)=d_u(i,k)
      u3(i,j,k)=rateu(i,k)

      ratew_w(i,k)=wplus(i,k)-w1(i,j,k)

      tp1=1.+fsr*int(k/nze)
	  tp2=1.-fsr*int(k/nze)

        ratew_u(i,k)=        
     +  (uplus(i,k)*tp1+uplus(i+1,k)*tp1+
     +   uplus(i,k+1)*tp2+uplus(i+1,k+1)*tp2)/4. - 
     +  (u1(i,j,k)*tp1+u1(i+1,j,k)*tp1+
     +   u1(i,j,k+1)*tp2+u1(i+1,j,k+1)*tp2)/4.
     
c        if(fsr.eq.1.and.k.eq.nze) ratew_u(i,k)=
c     +   (uplus(i,k)+uplus(i+1,k))/2. - (u1(i,j,k)+u1(i+1,j,k))/2.
      
     
        ratew(i,k) = sqrt(ratew_w(i,k)**2 + ratew_u(i,k)**2)

        d_w(i,k) = d_w(i,k) + dt*ratew(i,k)


50      continue

       return
        end	  


        subroutine xy4th_plus(nxb,nxe,nyb,nye,nzb,nze,dh,dt,jf)

c     4nd order finite-difference of xy at t+1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)
c     jf    fault plane location (y direction)  (integer)   (sent)

      use parstat

      jj=jf-3
      dth = dt/dh
      c1 = 9./8.
      c2 = -1./24.

      do 50 k= nzb,nze
      do 50 j= nyb,nye
      do 50 i= nxb,nxe

      xmu = 2./(mu(i,j,k)+mu(i,j,k-1))

      xy_f(i,j-jj,k) = xy_f(i,j-jj,k) + dth*xmu*

     +(c1*( u1(i,j+1,k)  - u1(i,j,k)   )  +
     +c2*( u1(i,j+2,k)  - uplus(i,k) )  +

     +c1*( v1(i,j,k)  - v1(i-1,j,k)   )  +
     +c2*( v1(i+1,j,k)  - v1(i-2,j,k) ))

       xy(i,j,k) = xy_f(i,j-jj,k)
   50 continue

      return
      end
        
      subroutine yz4th_plus(nxb,nxe,nyb,nye,nzb,nze,dh,dt,jf)

c     4nd order finite-difference of yz at t+1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)
c     jf    fault plane location (y direction)  (integer)   (sent)

      use parstat

      jj=jf-3
      dth = dt/dh
      c1 = 9./8.
      c2 = -1./24.

      do 50 k= nzb,nze
      do 50 j= nyb,nye
      do 50 i= nxb,nxe

      xmu = 2./(mu(i,j,k) + mu(i+1,j,k))

      yz_f(i,j-jj,k) = yz_f(i,j-jj,k) + dth*xmu*

     +(c1*( v1(i,j,k+1)  - v1(i,j,k)   )  +
     +c2*( v1(i,j,k+2)  - v1(i,j,k-1) )  +

     +c1*( w1(i,j+1,k)    - w1(i,j,k) )  +
     +c2*( w1(i,j+2,k)  - wplus(i,k) ))
     
      yz(i,j,k) = yz_f(i,j-jj,k)
   50 continue

      return
      end

       subroutine xyz2th(nxb,nxe,nyb,nye,nzb,nze,dh,dt,jf)

c     2th order finite-difference of normal stresses at t+1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)
c     jf    fault plane location (y direction)  (integer)   (sent)

      use parstat
      jj=jf-3
      dth = dt/dh
      do 50 k= nzb,nze
      do 50 j= nyb,nye
      do 50 i= nxb,nxe

      xl=8./(lam(i,j,k)+lam(i+1,j,k) +
     +       lam(i,j-1,k)+lam(i+1,j-1,k) +
     +       lam(i,j,k-1)+lam(i+1,j,k-1) +
     +       lam(i,j-1,k-1)+lam(i+1,j-1,k-1))     
     
      xm=8./(mu(i,j,k)+mu(i+1,j,k) +
     +       mu(i,j-1,k)+mu(i+1,j-1,k) +
     +       mu(i,j,k-1)+mu(i+1,j,k-1) +
     +       mu(i,j-1,k-1)+mu(i+1,j-1,k-1))      


      a = xl + 2.*xm
      b = xl

c     find xx stress

      xx_f(i,j-jj,k)=xx_f(i,j-jj,k) + dth*a*(u1(i+1,j,k) - u1(i,j,k))+
     +      dth*b*(v1(i,j,k) - v1(i,j-1,k) + w1(i,j,k) - w1(i,j,k-1))

c    find yy stress

      yy_f(i,j-jj,k)=yy_f(i,j-jj,k) + dth*a*(v1(i,j,k) - v1(i,j-1,k))+
     +      dth*b*(u1(i+1,j,k) - u1(i,j,k) + w1(i,j,k) - w1(i,j,k-1))

c    find zz stress

      zz_f(i,j-jj,k)=zz_f(i,j-jj,k) + dth*a*(w1(i,j,k) - w1(i,j,k-1))+
     +      dth*b*(u1(i+1,j,k) - u1(i,j,k) + v1(i,j,k) - v1(i,j-1,k))
      
       xx(i,j,k)=xx_f(i,j-jj,k)
       yy(i,j,k)=yy_f(i,j-jj,k)
       zz(i,j,k)=zz_f(i,j-jj,k)
   50 continue

      return
      end
      
      
      subroutine xy2th_plus(nxb,nxe,nyb,nye,nzb,nze,dh,dt,jf)

c     2nd order finite-difference of xy at t+1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)
c     jf    fault plane location (y direction)  (integer)   (sent)

      use parstat
      jj=jf-3
      
      dth = dt/dh
      do 50 k= nzb,nze
      do 50 j= nyb,nye
      do 50 i= nxb,nxe

      xmu = 2./(mu(i,j,k)+mu(i,j,k-1))
      
      xy_f(i,j-jj,k) = xy_f(i,j-jj,k) + dth*xmu*
     +        (u1(i,j+1,k) - uplus(i,k)+ v1(i,j,k) - v1(i-1,j,k))
      xy(i,j,k) = xy_f(i,j-jj,k)
   50 continue

      return
      end
      
      subroutine xy2th_minus(nxb,nxe,nyb,nye,nzb,nze,dh,dt,jf)

c     2nd order finite-difference of xy at t+1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)
c     jf    fault plane location (y direction)  (integer)   (sent)

      use parstat
      jj=jf-3
      
      dth = dt/dh
      do 50 k= nzb,nze
      do 50 j= nyb,nye
      do 50 i= nxb,nxe

      xmu = 2./(mu(i,j,k)+mu(i,j,k-1))
      
      xy_f(i,j-jj,k) = xy_f(i,j-jj,k) + dth*xmu*
     +        (u1(i,j+1,k) - u1(i,j,k)+ v1(i,j,k) - v1(i-1,j,k))
      xy(i,j,k) = xy_f(i,j-jj,k)
   50 continue

      return
      end

      subroutine yz2th_plus(nxb,nxe,nyb,nye,nzb,nze,dh,dt,jf)

c     2nd order finite-difference of yz at t+1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)
c     jf    fault plane location (y direction)  (integer)   (sent)

      use parstat

      jj=jf-3
      dth = dt/dh

      do 50 k= nzb,nze
      do 50 j= nyb,nye
      do 50 i= nxb,nxe

      xmu = 2./(mu(i,j,k) + mu(i+1,j,k))

      yz_f(i,j-jj,k) = yz_f(i,j-jj,k) + dth*xmu*
     +   (v1(i,j,k+1) - v1(i,j,k) + w1(i,j+1,k)- wplus(i,k))
     
      yz(i,j,k) = yz_f(i,j-jj,k)
   50 continue

      return
      end


      subroutine yz2th_minus(nxb,nxe,nyb,nye,nzb,nze,dh,dt,jf)

c     2nd order finite-difference of yz at t+1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)
c     jf    fault plane location (y direction)  (integer)   (sent)

      use parstat

      jj=jf-3
      dth = dt/dh

      do 50 k= nzb,nze
      do 50 j= nyb,nye
      do 50 i= nxb,nxe

      xmu = 2./(mu(i,j,k) + mu(i+1,j,k))

      yz_f(i,j-jj,k) = yz_f(i,j-jj,k) + dth*xmu*
     +   (v1(i,j,k+1) - v1(i,j,k) + w1(i,j+1,k)- w1(i,j,k))
     
      yz(i,j,k) = yz_f(i,j-jj,k)
   50 continue

      return
      end
      
      subroutine xzedge(nxb,nxe,nyb,nye,nzb,nze)

c     4nd order finite-difference of xz at t+1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)
c     jf    fault plane location (y direction)  (integer)   (sent)

      use parstat

      do 50 k= nzb,nze
      do 50 j= nyb,nye
      do 50 i= nxb,nxe

      xzplus(i,k) = xz(i,j,k)                 
      
   50 continue

      return
      end
      
       subroutine xyzedge(nxb,nxe,nyb,nye,nzb,nze)

c     2th order finite-difference of normal stresses at t+1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)
c     jf    fault plane location (y direction)  (integer)   (sent)

      use parstat
      do 50 k= nzb,nze
      do 50 j= nyb,nye
      do 50 i= nxb,nxe

       xxplus(i,k)=xx(i,j,k)
       zzplus(i,k)=zz(i,j,k)
       
   50 continue

      return
      end
      

      
        subroutine xyzsplit(nxb,nxe,nyb,nzb,nze,dh,dt)     

c     2th order finite-difference of split normal stresses at t+1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)
c     nyb    fault plane location (y direction)  (integer)   (sent)

      use parstat
      jj=nyb-3
      dth = dt/dh
      j= nyb      
      do 50 k= nzb,nze
      do 50 i= nxb,nxe

      xl_p=4./(lam(i,j,k)+lam(i+1,j,k) +
     +       lam(i,j,k-1)+lam(i+1,j,k-1))     
     
      xm_p=4./(mu(i,j,k)+mu(i+1,j,k) +
     +       mu(i,j,k-1)+mu(i+1,j,k-1))      


      xl_m=4./(lam(i,j-1,k)+lam(i+1,j-1,k) +
     +       lam(i,j-1,k-1)+lam(i+1,j-1,k-1))     
     
      xm_m=4./(mu(i,j-1,k)+mu(i+1,j-1,k) +
     +       mu(i,j-1,k-1)+mu(i+1,j-1,k-1))      
      
       du=u1(i+1,j,k) - u1(i,j,k)
       dw=w1(i,j,k) - w1(i,j,k-1)
              
       dupl=uplus(i+1,k) - uplus(i,k)
       dwpl=wplus(i,k) - wplus(i,k-1)
      
       a_m=2*(xl_m + 2.*xm_m)/dh
       a_p=2*(xl_p + 2.*xm_p)/dh       
       
       b_m=xl_m*(du+dw)/dh
       b_p=xl_p*(dupl+dwpl)/dh          

       v11=(a_p*v1(i,j,k)+a_m*v1(i,j-1,k)+b_p-b_m)/(a_p+a_m)       
c
c---V Particle velocity and displacement component Y (No Fault Opening)
c
       vplus(i,k) = v11
       vmin(i,k) = v11
       dvplus(i,k) = dvplus(i,k) + dt*vplus(i,k)
       dvmin(i,k) = dvmin(i,k) + dt*vmin(i,k)
c
       dv=v11-v1(i,j-1,k)
       dvpl=v1(i,j,k)-v11
       
      a_m = xl_m + 2.*xm_m
      b_m = xl_m

      a_p = xl_p + 2.*xm_p
      b_p = xl_p


        
c  (+) side             
      xxplus(i,k)=xxplus(i,k) + dth*a_p*(dupl)+
     +      dth*b_p*(2*dvpl + dwpl)

      zzplus(i,k)=zzplus(i,k) + dth*a_p*(dwpl)+
     +      dth*b_p*(2*dvpl + dupl)

      yy_p=yy_f(i,j-jj,k) + 2*dth*a_p*(dvpl)+
     +      dth*b_p*(dupl + dwpl)

c  (-) side             
      xx_f(i,j-jj,k)=xx_f(i,j-jj,k) + dth*a_m*(du)+
     +      dth*b_m*(2*dv + dw)

      zz_f(i,j-jj,k)=zz_f(i,j-jj,k) + dth*a_m*(dw)+
     +      dth*b_m*(2*dv + du)

      yy_m=yy_f(i,j-jj,k) + 2*dth*a_m*(dv)+
     +      dth*b_m*(du + dw)

       yy_f(i,j-jj,k)=(yy_m+yy_p)/2.
       
       xx(i,j,k)=xx_f(i,j-jj,k)
       yy(i,j,k)=yy_f(i,j-jj,k)
       zz(i,j,k)=zz_f(i,j-jj,k)
   50 continue
        
       return
       end
       
        subroutine xzsplit(nxb,nxe,nyb,nzb,nze,fsr,dh,dt)       


c     4nd order finite-difference of xz split on the fault at t+1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)
c     nyb    fault plane location (y direction)  (integer)   (sent)

      use parstat
      integer :: fsr      
      jj=nyb-3
      dth = dt/dh
      j= nyb  

      c1 = 9./8.
      c2 = -1./24.

      do 50 k= nzb,nze
      do 50 i= nxb,nxe

      xmu_p = mu(i,j,k)
      xmu_m = mu(i,j-1,k)      
      
c       if(fsr.eq.1.and.k.eq.nze-1) then
c        xz_f(i,j-jj,k) = 0.0
c        xzplus(i,k)=0.0
c      xzplus(i,k) = xzplus(i,k) + dth*xmu_p*
c     +(uplus(i,k+1) - uplus(i,k) + wplus(i,k) - wplus(i-1,k))
c
c      xz_f(i,j-jj,k) = xz_f(i,j-jj,k) + dth*xmu_m*
c     +(u1(i,j,k+1) - u1(i,j,k) + w1(i,j,k) - w1(i-1,j,k))
c
c       else

      xzplus(i,k) = xzplus(i,k) + dth*xmu_p*

     +(c1*( uplus(i,k+1)  - uplus(i,k)   )  +
     +c2*( uplus(i,k+2)  - uplus(i,k-1) )  +

     +c1*( wplus(i,k)    - wplus(i-1,k) )  +
     +c2*( wplus(i+1,k)  - wplus(i-2,k) ))

      xz_f(i,j-jj,k) = xz_f(i,j-jj,k) + dth*xmu_m*

     +(c1*( u1(i,j,k+1)  - u1(i,j,k)   )  +
     +c2*( u1(i,j,k+2)  - u1(i,j,k-1) )  +

     +c1*( w1(i,j,k)    - w1(i-1,j,k) )  +
     +c2*( w1(i+1,j,k)  - w1(i-2,j,k) ))

c       endif
      xz(i,j,k) = xz_f(i,j-jj,k)
      
   50 continue

      return
      end


       subroutine resttract_u(nxb,nxe,nyb,nzb,nze,fsr,dh,dt,eta)
     
      use parstat
      real :: eta      
      integer :: fsr      
      jj=nyb-3
      dth = dt/dh
      j= nyb  

c u component        

      do 50 k= nzb,nze
      do 50 i= nxb,nxe
       area=dh*dh
  
c (-) side                        
       d = 0.5*(d1(i,j-1,k)+d1(i,j-1,k-1))
       smas_m=d*area*dh/2.

       dxx=xx(i,j,k)-xx(i-1,j,k)       
       dxy=-xy(i,j-1,k)
       dxz=xz(i,j,k)-xz(i,j,k-1)

        Rrate= (area*(dxy + dxx/2. + dxz/2.)
     +      - Rminusue(i,k))/dt
       
        Rminusue(i,k)=area*(dxy + dxx/2. + dxz/2.)
                
        Rminusu(i,k)= Rminusue(i,k) + dt*eta*Rrate

c (+) side        
       d = 0.5*(d1(i,j,k)+d1(i,j,k-1))
       smas_p=d*area*dh/2.   

       dxx=xxplus(i,k)-xxplus(i-1,k)       
       dxy=xy(i,j,k)
       dxz=xzplus(i,k)-xzplus(i,k-1)
        
        Rrate= (area*(dxy + dxx/2. + dxz/2.)
     +      - Rplusue(i,k))/dt
       
        Rplusue(i,k)=area*(dxy + dxx/2. + dxz/2.)
                
        Rplusu(i,k)= Rplusue(i,k) + dt*eta*Rrate

c traction        
         temp=area*(smas_p+smas_m)
          deltau=uplus(i,k)-u1(i,j,k)
          
        tru1(i,k)=strisx(i,k)+(smas_p*smas_m*deltau/dt+
     +           smas_m*Rplusu(i,k)-smas_p*Rminusu(i,k))/temp
        tru2(i,k)=tru1(i,k)
        
   50 continue       
       return
       end

       subroutine resttract_w(nxb,nxe,nyb,nzb,nze,fsr,dh,dt,eta)
     
      use parstat
      real :: eta      
      integer :: fsr      
      
      jj=nyb-3
      dth = dt/dh
      j= nyb  

c w component        

      do 50 k= nzb,nze
      do 50 i= nxb,nxe
       area=dh*dh/(1.+fsr*int(k/nze))             
c       if(fsr.eq.1.and.k.eq.nze) area=dh*dh/2.   
c (-) side
       d = 0.5*(d1(i,j-1,k)+d1(i+1,j-1,k))
       smas_m=d*area*dh/2.

       dzz=zz(i,j,k+1)-zz(i,j,k)       
       dyz=-yz(i,j-1,k)
       dxz=xz(i+1,j,k)-xz(i,j,k)

        Rrate= (area*(dyz + dzz/2. + dxz/2.)
     +      - Rminuswe(i,k))/dt
       
        Rminuswe(i,k)=area*(dyz + dzz/2. + dxz/2.)
                
        Rminusw(i,k)= Rminuswe(i,k) + dt*eta*Rrate

c (+) side        
       d = 0.5*(d1(i,j,k)+d1(i+1,j,k))
       smas_p=d*area*dh/2.   

       dzz=zzplus(i,k+1)-zzplus(i,k)       
       dyz=yz(i,j,k)
       dxz=xzplus(i+1,k)-xzplus(i,k)

        Rrate= (area*(dyz + dzz/2. + dxz/2.)
     +      - Rpluswe(i,k))/dt
       
        Rpluswe(i,k)=area*(dyz + dzz/2. + dxz/2.)
                
        Rplusw(i,k)= Rpluswe(i,k) + dt*eta*Rrate

c traction        
         temp=area*(smas_p+smas_m)
          deltaw=wplus(i,k)-w1(i,j,k)
          
        trw1(i,k)=strisz(i,k) + (smas_p*smas_m*deltaw/dt+
     +           smas_m*Rplusw(i,k)-smas_p*Rminusw(i,k))/temp
        trw2(i,k)=trw1(i,k)
        
   50 continue       
       return
       end
       
     
       
       subroutine sourcestres(it,nxb,nxe,nyb,nzb,nze,nzt,fsr,dh,dt)
       
c     it     current time step in FD code (int)

      use parstat
      real    :: mus,mud,iyy      
      integer :: fz,fsr
      iyy=1.0
      fz=0      
      if(nzt.gt.fz_d) fz=1            
      
      j= nyb       
      do 50 k= nzb,nze
      do 50 i= nxb,nxe
c u componet (strike)              
      itmpbroken= brokenu(i,k)    
	  tx=tru1(i,k)
	  tz=(trw1(i,k)+trw1(i-1,k)+trw1(i,k-1)+trw1(i-1,k-1))/4.	  

	  tot1=sqrt(tx**2 + tz**2) 

c friction model (slip weakening model)
          dd = d_u(i,k)
	  tyy =-striny(i,k)+yy(i,j,k)/2.+yy(i-1,j,k)/2.
	  if(tyy.gt.0) tyy=0.0
c
c---V Fault shear traction Tyz
c
          trv1(i,k) = tyy
c
    	  tot= -tyy*mu_s(i,k)
      if(dd.le.d0(i,k).and.d0(i,k).gt.0.0) 
     +   tot= -tyy*((1-dd/d0(i,k))*mu_s(i,k) + mu_d(i,k)*dd/d0(i,k))
	  if(dd.gt.d0(i,k)) tot= -tyy*mu_d(i,k)

c checking if fault break

	  if(tot1.gt.tot) then
	    brokenu(i,k)=1
        trupt(i,k)=it*dt*(1-itmpbroken) + trupt(i,k)*itmpbroken        
c=  safe in variable originally set in the code
         rupt(i,j,k)=trupt(i,k)
    	 broken(i,j,k)=brokenu(i,k)         
c=        
    	tru2(i,k) = tot*tru2(i,k)/tot1
	  end if   
	  
c w componet (dip)              
	  
	  tp1=1.+fsr*int(k/nze)
	  tp2=1.-fsr*int(k/nze)

      mus=(mu_s(i,k)*tp1+mu_s(i+1,k)*tp1+mu_s(i,k+1)*tp2+mu_s(i+1,k+1)*tp2)/4.
      mud=(mu_d(i,k)*tp1+mu_d(i+1,k)*tp1+mu_d(i,k+1)*tp2+mu_d(i+1,k+1)*tp2)/4.      
      d0_w=(d0(i,k)*tp1+d0(i+1,k)*tp1+d0(i,k+1)*tp2+d0(i+1,k+1)*tp2)/4.            
	  tx=(tru1(i,k)*tp1+tru1(i+1,k)*tp1+tru1(i,k+1)*tp2+tru1(i+1,k+1)*tp2)/4.
	  
c	  if(fsr.eq.1.and.k.eq.nze) then
c       mus=(mu_s(i,k)+mu_s(i+1,k))/2.
c       mud=(mu_d(i,k)+mu_d(i+1,k))/2.      
c       d0_w=(d0(i,k)+d0(i+1,k))/2.            
c	   tx=(tru1(i,k)+tru1(i+1,k))/2.
c	   iyy=0.0
c	  endif	

       
	  tz=trw1(i,k)

	  tot1=sqrt(tx**2 + tz**2) 

c friction model (slip weakening model)

          dd = d_w(i,k)
	  tyy =-striny(i,k)+yy(i,j,k)*tp1/2.+yy(i,j,k+1)*tp2/2.
	  if(tyy.gt.0) tyy=0.0	  
    	  tot= -tyy*mus
      if(dd.le.d0_w.and.d0_w.gt.0.0) 
     +   tot= -tyy*((1-dd/d0_w)*mus + mud*dd/d0_w)
	  if(dd.gt.d0_w) tot= -tyy*mud

c checking if fault break

	  if(tot1.gt.tot) then
	    brokenw(i,k)=1
    	trw2(i,k) = tot*trw2(i,k)/tot1
	  end if   	  
	  
50        continue	  

      return
      end


      subroutine fvelxysplit(nxb,nxe,nzt)

c     free-surface B.C. for u split velocities of SGSN dynamic model

      use parstat

      do 10 i=nxb-1,nxe+1

      uplus(i,nzt+1) = uplus(i,nzt) -
     +(wplus(i,nzt) - wplus(i-1,nzt))

   10 continue

      return
      end

      subroutine fvelzsplit(nxb,nxe,nyb,nzt)

c     free-surface B.C. for w split velocities


      use parstat

      j=nyb
      do 30 i=nxb-1,nxe+1
      xl = 2.*(lam(i,j,nzt)+lam(i+1,j,nzt))
      xl2m = xl + 2.*2.*(mu(i,j,nzt)+mu(i+1,j,nzt))
     
      wplus(i,nzt+1) = wplus(i,nzt-1) - (xl/xl2m)*(

     +(uplus(i+1,nzt+1) - uplus(i,nzt+1)) +
     +(uplus(i+1,nzt) - uplus(i,nzt)) +

     +(v1(i,j,nzt+1) - v1(i,j-1,nzt+1)) +
     +(v1(i,j,nzt) - v1(i,j-1,nzt)))

      xl = 2./(lam(i,j-1,nzt)+lam(i+1,j-1,nzt))
      xl2m = xl + 2.*2./(mu(i,j-1,nzt)+mu(i+1,j-1,nzt))
     
      w1(i,j,nzt+1) = w1(i,j,nzt-1) - (xl/xl2m)*(

     +(u1(i+1,j,nzt+1) - u1(i,j,nzt+1)) +
     +(u1(i+1,j,nzt) - u1(i,j,nzt)) +

     +(v1(i,j,nzt+1) - v1(i,j-1,nzt+1)) +
     +(v1(i,j,nzt) - v1(i,j-1,nzt)))

   30 continue

      return
      end

      subroutine fstrsplit(nxb,nxe,nyb,nzt)  

                                                               
c     free-surface B.C. for splits stresses
                                                           
      use parstat
                                                           
      j=nyb
      do 10 i=nxb-2,nxe+2
                                                      
c     asymmetry reflection above free surface
c (+) side                                                             
      zzplus(i,nzt+1) = -zzplus(i,nzt)
      zzplus(i,nzt+2) = -zzplus(i,nzt-1)
                                                            
      xzplus(i,nzt+1) = -xzplus(i,nzt-1)
      xzplus(i,nzt+2) = -xzplus(i,nzt-2)                                                                 
                                                                
      xzplus(i,nzt) = 0.
      
c (-) side                                                             
      zz(i,j,nzt+1) = -zz(i,j,nzt)
      zz(i,j,nzt+2) = -zz(i,j,nzt-1)
                                                            
      xz(i,j,nzt+1) = -xz(i,j,nzt-1)
      xz(i,j,nzt+2) = -xz(i,j,nzt-2)                                                                 
                                                                
      xz(i,j,nzt) = 0.
                                                                
   10 continue
                                                         
      return
      end

      

      subroutine fsrcprm(nxb,nxe,nzt)  
                                                               
c     source parameters above free-surface
                                                           
      use parstat
                                                           
      do 10 i=nxb-2,nxe+2
       d0(i,nzt+1) = d0(i,nzt)
       mu_s(i,nzt+1) = mu_s(i,nzt)
       mu_d(i,nzt+1) = mu_d(i,nzt)      
   10 continue                                                         
      return
      end
