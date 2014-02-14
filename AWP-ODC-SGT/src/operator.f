CCC   'operator.f' COMPUTES INTERNAL REGION VELOCITIES AND STRESSES


      subroutine dvel(dh,dt)

c     4th order finite-difference of velocity components at t+1/2

c     nxt   nodal points in x dir  (integer)(sent)
c     nyt   nodal points in y dir  (integer)(sent)
c     nzt   nodal points in z dir  (integer)(sent)
c     dh    spatial discretization (real)   (sent)
c     dt    temporal discretization(real)   (sent)
c     nd    damping zone width     (integer)(sent) 

      use parstat

      call uxx1(xi(18),xf(18),yi(18),yf(18),zi(18),zf(18), dh,dt)
      call vyy1(xi(18),xf(18),yi(18),yf(18),zi(18),zf(18), dh,dt)
      call wzz1(xi(18),xf(18),yi(18),yf(18),zi(18),zf(18), dh,dt)

      return
      end



      subroutine dstr(dh,dt)

c     4th order finite-difference of stress components

c     nxt   nodal points in x dir          (integer)(sent)
c     nyt   nodal points in y dir          (integer)(sent)
c     nzt   nodal points in z dir          (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      use parstat

      call xyz1(xi(18),xf(18),yi(18),yf(18),zi(18),zf(18), dh,dt)
      call  xz1(xi(18),xf(18),yi(18),yf(18),zi(18),zf(18), dh,dt)
      call  yz1(xi(18),xf(18),yi(18),yf(18),zi(18),zf(18), dh,dt)
      call  xy1(xi(18),xf(18),yi(18),yf(18),zi(18),zf(18), dh,dt)

      return
      end



      subroutine dvelc(dh,dt,nxt,nyt,nzt)

c     4th order finite-difference of velocity components at t+1/2

c     nxt   nodal points in x dir  (integer)(sent)
c     nyt   nodal points in y dir  (integer)(sent)
c     nzt   nodal points in z dir  (integer)(sent)
c     dh    spatial discretization (real)   (sent)
c     dt    temporal discretization(real)   (sent)
c     nd    damping zone width     (integer)(sent)

      use parstat

      call uxx1(1,nxt,1,nyt,1,nzt, dh,dt)
      call vyy1(1,nxt,1,nyt,1,nzt, dh,dt)
      call wzz1(1,nxt,1,nyt,1,nzt, dh,dt)

      return
      end



      subroutine dstrc(dh,dt,nxt,nyt,nzt)

c     4th order finite-difference of stress components

c     nxt   nodal points in x dir          (integer)(sent)
c     nyt   nodal points in y dir          (integer)(sent)
c     nzt   nodal points in z dir          (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      use parstat

      call xyz1(1,nxt,1,nyt,1,nzt, dh,dt)
      call  xz1(1,nxt,1,nyt,1,nzt, dh,dt)
      call  yz1(1,nxt,1,nyt,1,nzt, dh,dt)
      call  xy1(1,nxt,1,nyt,1,nzt, dh,dt)

      return
      end



      subroutine uxx1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     4nd order finite-difference of u1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      use parstat
      
      dth = dt/dh
      c1 = 9./8.
      c2 = -1./24.

c     Find u-displacement fields at time t+1/2

      do kk= nzb,nze,kblock
      do jj= nyb,nye,jblock
      do 50 k=kk,MIN(kk+kblock-1,nze)
      do 50 j=jj,MIN(jj+jblock-1,nye)
      do 50 i= nxb,nxe

C == FOR SGSN DYNAMIC FAULT MODEL
      d = 0.25*((d1(i,j,k)+d1(i,j-1,k))+
     +          (d1(i,j,k-1)+d1(i,j-1,k-1)))
C ==
     
      u1(i,j,k)=u1(i,j,k)+(dth/d)*(

     +c1*(xx(i,j,k)-xx(i-1,j,k))+ 
     +c2*(xx(i+1,j,k)-xx(i-2,j,k))+

     +c1*(xy(i,j,k)-xy(i,j-1,k))+
     +c2*(xy(i,j+1,k)-xy(i,j-2,k))+

     +c1*(xz(i,j,k)-xz(i,j,k-1))+
     +c2*(xz(i,j,k+1)-xz(i,j,k-2)))

   50 continue
      enddo
      enddo
      return
      end


 
      subroutine vyy1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     4nd order finite-difference of v1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      use parstat

      dth = dt/dh
      c1 = 9./8.
      c2 = -1./24.

c     find v-displacement fields at time t+1/2
       
      do kk= nzb,nze,kblock
      do jj= nyb,nye,jblock
      do 50 k=kk,MIN(kk+kblock-1,nze)
      do 50 j=jj,MIN(jj+jblock-1,nye)
      do 50 i= nxb,nxe

C == FOR SGSN DYNAMIC FAULT MODEL
      d = 0.25*((d1(i,j,k)+d1(i+1,j,k))+
     +          (d1(i,j,k-1)+d1(i+1,j,k-1)))

C ==
      v1(i,j,k)=v1(i,j,k)+(dth/d)* (

     +c1*(xy(i+1,j,k)-xy(i,j,k))+
     +c2*(xy(i+2,j,k)-xy(i-1,j,k))+

     +c1*(yy(i,j+1,k)-yy(i,j,k))+
     +c2*(yy(i,j+2,k)-yy(i,j-1,k))+  

     +c1*(yz(i,j,k)-yz(i,j,k-1))+
     +c2*(yz(i,j,k+1)-yz(i,j,k-2)))
C ==

   50 continue
      enddo
      enddo
      return
      end


 
      subroutine wzz1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     4nd order finite-difference of w1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      use parstat

      dth = dt/dh
      c1 = 9./8.
      c2 = -1./24.


c     find w-displacement fields at time t+1/2

      do kk= nzb,nze,kblock
      do jj= nyb,nye,jblock
      do 50 k=kk,MIN(kk+kblock-1,nze)
      do 50 j=jj,MIN(jj+jblock-1,nye)
      do 50 i= nxb,nxe

C == FOR SGSN DYNAMIC FAULT MODEL
      d = 0.25*((d1(i,j,k)+d1(i+1,j,k))+
     +          (d1(i,j-1,k)+d1(i+1,j-1,k)))

C ==
      w1(i,j,k)=w1(i,j,k)+(dth/d)*(

     +c1*(xz(i+1,j,k)-xz(i,j,k))+
     +c2*(xz(i+2,j,k)-xz(i-1,j,k))+

     +c1*(yz(i,j,k)-yz(i,j-1,k))+
     +c2*(yz(i,j+1,k)-yz(i,j-2,k))+

     +c1*(zz(i,j,k+1)-zz(i,j,k))+
     +c2*(zz(i,j,k+2)-zz(i,j,k-1)))
C ==

   50 continue
      enddo
      enddo
      return
      end

      subroutine fvelxy(nxt,nyt,nzt)

c     free-surface B.C. for velocities

c     nxt   nodal points in x dir (integer)(sent)
c     nyt   nodal points in y dir (integer)(sent)
c     nzt   nodal points in z dir (integer)(sent)
      use parstat

      do 10 j=1,nyt
      do 10 i=1,nxt

      u1(i,j,nzt+1) = u1(i,j,nzt) -
     +(w1(i,j,nzt) - w1(i-1,j,nzt))

   10 continue

      do 20 j=1,nyt
      do 20 i=1,nxt

      v1(i,j,nzt+1) = v1(i,j,nzt) -
     +(w1(i,j+1,nzt) - w1(i,j,nzt))

   20 continue

      return
      end

      subroutine fvelz(nxt,nyt,nzt)

c     free-surface B.C. for velocities

c     nxt   nodal points in x dir (integer)(sent)
c     nyt   nodal points in y dir (integer)(sent)
c     nzt   nodal points in z dir (integer)(sent)

      use parstat

      do jj= 1,nyt,kblock
      do ii= 1,nxt,jblock
      do 30 j=jj,MIN(jj+kblock-1,nyt)
      do 30 i=ii,MIN(ii+jblock-1,nxt)

      xl = 1./lam(i,j,nzt)
      xl2m = 2./mu(i,j,nzt)+ xl
C == FOR SGSN DYNAMIC FAULT MODEL

C ==      
      w1(i,j,nzt+1) = w1(i,j,nzt-1) - (xl/xl2m)*(

     +(u1(i+1,j,nzt+1) - u1(i,j,nzt+1)) +
     +(u1(i+1,j,nzt) - u1(i,j,nzt)) +

     +(v1(i,j,nzt+1) - v1(i,j-1,nzt+1)) +
     +(v1(i,j,nzt) - v1(i,j-1,nzt)))

   30 continue
      enddo
      enddo
      return
      end



      subroutine fvelzma(nxt,nyt,nzt)

CCC   DAY< BRADLEY< OLSEN VERSION OF FS WITH MEDIA AVERAGING

      use parstat      

      do jj= 1,nyt,kblock
      do ii= 1,nxt,jblock
      do 30 j=jj,MIN(jj+kblock-1,nyt)
      do 30 i=ii,MIN(ii+jblock-1,nxt)

      al=0.5*(1./lam(i,j,nzt)+1./lam(i,j,nzt+1) + 2.*(1./mu(i,j,nzt)+1./mu(i,j,nzt+1)))
      amu=0.5*(1./mu(i,j,nzt)+1./mu(i,j,nzt+1))

      a=al-2*amu
      b=al

      al=0.5*(1./lam(i,j,nzt)+1./lam(i,j,nzt-1) + 2.*(1./mu(i,j,nzt)+1./mu(i,j,nzt-1)))
      amu=0.5*(1./mu(i,j,nzt)+1./mu(i,j,nzt-1))

      a1=al-2*amu
      b1=al

      w1(i,j,nzt+1)=w1(i,j,nzt)-
     +              (a/b)*
     +              (u1(i+1,j,nzt+1)-u1(i,j,nzt+1)+
     +              v1(i,j,nzt+1)-v1(i,j-1,nzt+1))-
     +              (w1(i,j,nzt)-w1(i,j,nzt-1))-
     +              (a1/b1)*
     +              (u1(i+1,j,nzt)-u1(i,j,nzt)+
     +              v1(i,j,nzt)-v1(i,j-1,nzt))

   30 continue
      enddo
      enddo
      return
      end

      subroutine xyz1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     4th order finite-difference of normal stresses at t+1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      use parstat

      dth = dt/dh
      c1 = 9./8.
      c2 = -1./24.

      do kk= nzb,nze,kblock
      do jj= nyb,nye,jblock
      do 50 k=kk,MIN(kk+kblock-1,nze)
      do 50 j=jj,MIN(jj+jblock-1,nye)
      do 50 i= nxb,nxe

C == FOR SGSN DYNAMIC FAULT MODEL
      xl=8./(lam(i,j,k)+lam(i+1,j,k) +
     +       lam(i,j-1,k)+lam(i+1,j-1,k) +
     +       lam(i,j,k-1)+lam(i+1,j,k-1) +
     +       lam(i,j-1,k-1)+lam(i+1,j-1,k-1))     
     
      xm=8./(mu(i,j,k)+mu(i+1,j,k) +
     +       mu(i,j-1,k)+mu(i+1,j-1,k) +
     +       mu(i,j,k-1)+mu(i+1,j,k-1) +
     +       mu(i,j-1,k-1)+mu(i+1,j-1,k-1))      
C==      
      a = xl + 2.*xm
      b = xl

c     find xx stress

      xx(i,j,k) = xx(i,j,k) + dth*a*

     +(c1*( u1(i+1,j,k) - u1(i,j,k)   )  +
     +c2*( u1(i+2,j,k) - u1(i-1,j,k) )     )

     ++dth*b*( c1*( v1(i,j,k)   - v1(i,j-1,k) )  +
     +c2*( v1(i,j+1,k) - v1(i,j-2,k) )  +

     +c1*( w1(i,j,k)   - w1(i,j,k-1) )  +
     +c2*( w1(i,j,k+1) - w1(i,j,k-2) ))


c    find yy stress

      yy(i,j,k) = yy(i,j,k) + dth*a*

     +(c1*( v1(i,j,k) - v1(i,j-1,k)   )  +
     +c2*( v1(i,j+1,k) - v1(i,j-2,k) )     )

     ++dth*b*( c1*( u1(i+1,j,k)   - u1(i,j,k) )  +
     +c2*( u1(i+2,j,k) - u1(i-1,j,k) )  +

     +c1*( w1(i,j,k)   - w1(i,j,k-1) )  +
     +c2*( w1(i,j,k+1) - w1(i,j,k-2) ))


c    find zz stress

      zz(i,j,k) = zz(i,j,k) + dth*a*

     +(c1*( w1(i,j,k) - w1(i,j,k-1)   )  +
     +c2*( w1(i,j,k+1) - w1(i,j,k-2) )     )

     ++dth*b*( c1*( u1(i+1,j,k)   - u1(i,j,k) )  +
     +c2*( u1(i+2,j,k) - u1(i-1,j,k) )  +

     +c1*( v1(i,j,k)   - v1(i,j-1,k) )  +
     +c2*( v1(i,j+1,k) - v1(i,j-2,k) ))

   50 continue
      enddo
      enddo
      return
      end


      subroutine xy1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     4nd order finite-difference of xy at t+1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      use parstat

      dth = dt/dh
      c1 = 9./8.
      c2 = -1./24.

      do kk= nzb,nze,kblock
      do jj= nyb,nye,jblock
      do 50 k=kk,MIN(kk+kblock-1,nze)
      do 50 j=jj,MIN(jj+jblock-1,nye)
      do 50 i= nxb,nxe

C == FOR SGSN DYNAMIC FAULT MODEL
      xmu = 2./(mu(i,j,k)+mu(i,j,k-1))
C ==

      xy(i,j,k) = xy(i,j,k) + dth*xmu*

     +(c1*( u1(i,j+1,k)  - u1(i,j,k)   )  +
     +c2*( u1(i,j+2,k)  - u1(i,j-1,k) )  +

     +c1*( v1(i,j,k)  - v1(i-1,j,k)   )  +
     +c2*( v1(i+1,j,k)  - v1(i-2,j,k) ))

   50 continue
      enddo
      enddo
      return
      end



      subroutine xz1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     4nd order finite-difference of xz at t+1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      use parstat

      dth = dt/dh
      c1 = 9./8.
      c2 = -1./24.

      do kk= nzb,nze,kblock
      do jj= nyb,nye,jblock
      do 50 k=kk,MIN(kk+kblock-1,nze)
      do 50 j=jj,MIN(jj+jblock-1,nye)
      do 50 i= nxb,nxe

C == FOR SGSN DYNAMIC FAULT MODEL  
      xmu = 2./(mu(i,j,k) + mu(i,j-1,k))
C ==

      xz(i,j,k) = xz(i,j,k) + dth*xmu*

     +(c1*( u1(i,j,k+1)  - u1(i,j,k)   )  +
     +c2*( u1(i,j,k+2)  - u1(i,j,k-1) )  +

     +c1*( w1(i,j,k)    - w1(i-1,j,k) )  +
     +c2*( w1(i+1,j,k)  - w1(i-2,j,k) ))

   50 continue
      enddo
      enddo
      return
      end



      subroutine yz1(nxb,nxe,nyb,nye,nzb,nze,dh,dt)

c     4nd order finite-difference of yz at t+1

c     nxb   starting point for FD in x dir (integer)(sent)
c     nxe   ending point for FD in x dir   (integer)(sent)
c     nyb   starting point for FD in y dir (integer)(sent)
c     nye   ending point for FD in y dir   (integer)(sent)
c     nzb   starting point for FD in z dir (integer)(sent)
c     nze   ending point for FD in z dir   (integer)(sent)
c     dh    spatial discretization         (real)   (sent)
c     dt    temporal discretization        (real)   (sent)

      use parstat

      dth = dt/dh
      c1 = 9./8.
      c2 = -1./24.

      do kk= nzb,nze,kblock
      do jj= nyb,nye,jblock
      do 50 k=kk,MIN(kk+kblock-1,nze)
      do 50 j=jj,MIN(jj+jblock-1,nye)
      do 50 i= nxb,nxe

c      xmu = 0.5*(mu(i,j,k) + mu(i,j+1,k))
c
c      xmu = mu(i,j,k)

C == FOR SGSN DYNAMIC FAULT MODEL  
      xmu = 2./(mu(i,j,k) + mu(i+1,j,k))
C ==


      yz(i,j,k) = yz(i,j,k) + dth*xmu*

     +(c1*( v1(i,j,k+1)  - v1(i,j,k)   )  +
     +c2*( v1(i,j,k+2)  - v1(i,j,k-1) )  +

     +c1*( w1(i,j+1,k)    - w1(i,j,k) )  +
     +c2*( w1(i,j+2,k)  - w1(i,j-1,k) ))

   50 continue
      enddo
      enddo
      return
      end



      subroutine fstr(nxt,nyt,nzt)
                                                               
c     free-surface B.C. for stresses
                                                           
c     nxt   nodal points in x dir  (integer)(sent)
c     nyt   nodal points in y dir  (integer)(sent)
c     nzt   nodal points in z dir  (integer)(sent)
c     dh    spatial discretization (real)   (sent)
c     dt    temporal discretization(real)   (sent)
                                                                  
      use parstat
                                                           
      do jj= 1,nyt,kblock
      do ii= 1,nxt,jblock
      do 10 j=jj,MIN(jj+kblock-1,nyt)
      do 10 i=ii,MIN(ii+jblock-1,nxt)
                                                      
c     asymmetry reflection above free surface
                                                             
      zz(i,j,nzt+1) = -zz(i,j,nzt)
      zz(i,j,nzt+2) = -zz(i,j,nzt-1)
                                                            
      xz(i,j,nzt+1) = -xz(i,j,nzt-1)
      xz(i,j,nzt+2) = -xz(i,j,nzt-2)
                                                               
      yz(i,j,nzt+1) = -yz(i,j,nzt-1)
      yz(i,j,nzt+2) = -yz(i,j,nzt-2)
                                                                 
c     xz & yz on free surface --> 0
                                                                
      xz(i,j,nzt) = 0.
      yz(i,j,nzt) = 0.
                                                                
   10 continue
      enddo
      enddo                                                   
      return
      end

