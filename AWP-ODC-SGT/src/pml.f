CCC   'pml.f' PERFECTLY MATCHED LAYERS EXTERNAL BOUNDARY CONDITIONS 


      subroutine inipml(dh,vel,ndt,arbc)

CCC   ROUTINE TO INITIALIZE DAMPING SERIES

C     DH	SPATIAL INCREMENT	   (REAL)  (SENT)	
C     VEL       EXAMPLE MODEL VELOCITY     (REAL)  (SENT)
C     DAMP	PML REFLECTION COEFFICIENT (REAL)  (SENT)
C     RH	HALF-NODE DAMPING ARRAY
C     RF	FULL_NODE DAMPING ARRAY

      use parstat

      dimension vel(2)
                     
C     REFLECTION COEFFICIENT - PML WIDTH RELATION (PML 1-20 ONLY!)
                                                  
      parameter(c0=0., c1=8./15., c2=-3./100., c3=1./1500.)
      if (rank==master) write(*,*) 'in the beginning of inipml, pht=', pht

      f = c0 + ndt*c1 + ndt*ndt*c2 + ndt*ndt*ndt*c3
      v = 1./(0.5*(1./vel(1) + 1./vel(2)))

      coef = arbc*f*v/(ndt*dh)

C     PML DAMPING ARRAYS

      do 50 i=1,ndt
        rf(ndt+1-i) = coef*i*i/(ndt*ndt)
   50 continue

C     ARBITRARY CALCULATION OF FINAL EDGE VALUE

      rh(1) = rf(1) + 0.50*(rf(1) - rf(2))

C     HALF NODE VALUES = AVERAGE OF FULL NODE VALUES
 
      do 60 i=2,ndt
        rh(i) = 0.5*(rf(i-1) + rf(i))
   60 continue

      return
      end
 

 
      subroutine pmlvel(nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE PML MODEL IN BOUNDARY ZONES

C     U(X), V(X) and W(X) ARE THE VELOCITY FIELDS
C     N(Y), S(Y), E(X), W(X) CODE INDICATES GEOGRAPHIC DIRECTION
C     B(Z) CODE INDICATES MODEL SPACE BOTTOM
C     SINGLE CODE IS A PLANE NORMAL TO THE DIRECTION AXIS (5 TOTAL)
C     COMBINATION OF 2 INDICATES AN EDGE (8 TOTAL)
C     COMBINATION OF 3 INDICATES A CORNER (4 TOTAL)
C     ROUTINE CALLS ORDERED PLANES, EDGES, CORNERS
 
C     N* 	FINAL COORDINATE ADDRESSES
C     DT	TIME INCREMENT
C     DH        SPATIAL INCREMENT

C     PLANES: N, S, W, E AND B

      use parstat

      call m_us(xi(1),xf(1),yi(1),yf(1),zi(1),zf(1), nx,ny,nz,dh,dt)
      call m_vs(xi(1),xf(1),yi(1),yf(1),zi(1),zf(1), nx,ny,nz,dh,dt)
      call m_ws(xi(1),xf(1),yi(1),yf(1),zi(1),zf(1), nx,ny,nz,dh,dt)

      call m_uw(xi(2),xf(2),yi(2),yf(2),zi(2),zf(2), nx,ny,nz,dh,dt)
      call m_vw(xi(2),xf(2),yi(2),yf(2),zi(2),zf(2), nx,ny,nz,dh,dt)
      call m_ww(xi(2),xf(2),yi(2),yf(2),zi(2),zf(2), nx,ny,nz,dh,dt)

      call m_un(xi(3),xf(3),yi(3),yf(3),zi(3),zf(3), nx,ny,nz,dh,dt)
      call m_vn(xi(3),xf(3),yi(3),yf(3),zi(3),zf(3), nx,ny,nz,dh,dt)
      call m_wn(xi(3),xf(3),yi(3),yf(3),zi(3),zf(3), nx,ny,nz,dh,dt)

      call m_ue(xi(4),xf(4),yi(4),yf(4),zi(4),zf(4), nx,ny,nz,dh,dt)
      call m_ve(xi(4),xf(4),yi(4),yf(4),zi(4),zf(4), nx,ny,nz,dh,dt)
      call m_we(xi(4),xf(4),yi(4),yf(4),zi(4),zf(4), nx,ny,nz,dh,dt)

      call m_ub(xi(5),xf(5),yi(5),yf(5),zi(5),zf(5), nx,ny,nz,dh,dt)
      call m_vb(xi(5),xf(5),yi(5),yf(5),zi(5),zf(5), nx,ny,nz,dh,dt)
      call m_wb(xi(5),xf(5),yi(5),yf(5),zi(5),zf(5), nx,ny,nz,dh,dt)

C     VERTICAL EDGES: SW, SE, NW AND NE

      call m_usw(xi(6),xf(6),yi(6),yf(6),zi(6),zf(6), nx,ny,nz,dh,dt)
      call m_vsw(xi(6),xf(6),yi(6),yf(6),zi(6),zf(6), nx,ny,nz,dh,dt)
      call m_wsw(xi(6),xf(6),yi(6),yf(6),zi(6),zf(6), nx,ny,nz,dh,dt)

      call m_use(xi(8),xf(8),yi(8),yf(8),zi(8),zf(8), nx,ny,nz,dh,dt)
      call m_vse(xi(8),xf(8),yi(8),yf(8),zi(8),zf(8), nx,ny,nz,dh,dt)
      call m_wse(xi(8),xf(8),yi(8),yf(8),zi(8),zf(8), nx,ny,nz,dh,dt)

      call m_unw(xi(7),xf(7),yi(7),yf(7),zi(7),zf(7), nx,ny,nz,dh,dt)
      call m_vnw(xi(7),xf(7),yi(7),yf(7),zi(7),zf(7), nx,ny,nz,dh,dt)
      call m_wnw(xi(7),xf(7),yi(7),yf(7),zi(7),zf(7), nx,ny,nz,dh,dt)

      call m_une(xi(9),xf(9),yi(9),yf(9),zi(9),zf(9), nx,ny,nz,dh,dt)
      call m_vne(xi(9),xf(9),yi(9),yf(9),zi(9),zf(9), nx,ny,nz,dh,dt)
      call m_wne(xi(9),xf(9),yi(9),yf(9),zi(9),zf(9), nx,ny,nz,dh,dt)

C     HORIZONTAL EDGES: SB, NB, WB AND EB

      call m_usb(xi(10),xf(10),yi(10),yf(10),zi(10),zf(10), nx,ny,nz,dh,dt)
      call m_vsb(xi(10),xf(10),yi(10),yf(10),zi(10),zf(10), nx,ny,nz,dh,dt)
      call m_wsb(xi(10),xf(10),yi(10),yf(10),zi(10),zf(10), nx,ny,nz,dh,dt)

      call m_unb(xi(12),xf(12),yi(12),yf(12),zi(12),zf(12), nx,ny,nz,dh,dt)
      call m_vnb(xi(12),xf(12),yi(12),yf(12),zi(12),zf(12), nx,ny,nz,dh,dt)
      call m_wnb(xi(12),xf(12),yi(12),yf(12),zi(12),zf(12), nx,ny,nz,dh,dt)

      call m_uwb(xi(11),xf(11),yi(11),yf(11),zi(11),zf(11), nx,ny,nz,dh,dt)
      call m_vwb(xi(11),xf(11),yi(11),yf(11),zi(11),zf(11), nx,ny,nz,dh,dt)
      call m_wwb(xi(11),xf(11),yi(11),yf(11),zi(11),zf(11), nx,ny,nz,dh,dt)

      call m_ueb(xi(13),xf(13),yi(13),yf(13),zi(13),zf(13), nx,ny,nz,dh,dt)
      call m_veb(xi(13),xf(13),yi(13),yf(13),zi(13),zf(13), nx,ny,nz,dh,dt)
      call m_web(xi(13),xf(13),yi(13),yf(13),zi(13),zf(13), nx,ny,nz,dh,dt)

C     BOTTOM CORNERS: SWB, NWB, SEB AND NEB
      
      call m_uswb(xi(14),xf(14),yi(14),yf(14),zi(14),zf(14), nx,ny,nz,dh,dt)
      call m_vswb(xi(14),xf(14),yi(14),yf(14),zi(14),zf(14), nx,ny,nz,dh,dt)
      call m_wswb(xi(14),xf(14),yi(14),yf(14),zi(14),zf(14), nx,ny,nz,dh,dt)

      call m_unwb(xi(15),xf(15),yi(15),yf(15),zi(15),zf(15), nx,ny,nz,dh,dt)
      call m_vnwb(xi(15),xf(15),yi(15),yf(15),zi(15),zf(15), nx,ny,nz,dh,dt)
      call m_wnwb(xi(15),xf(15),yi(15),yf(15),zi(15),zf(15), nx,ny,nz,dh,dt)

      call m_useb(xi(16),xf(16),yi(16),yf(16),zi(16),zf(16), nx,ny,nz,dh,dt)
      call m_vseb(xi(16),xf(16),yi(16),yf(16),zi(16),zf(16), nx,ny,nz,dh,dt)
      call m_wseb(xi(16),xf(16),yi(16),yf(16),zi(16),zf(16), nx,ny,nz,dh,dt)

      call m_uneb(xi(17),xf(17),yi(17),yf(17),zi(17),zf(17), nx,ny,nz,dh,dt)
      call m_vneb(xi(17),xf(17),yi(17),yf(17),zi(17),zf(17), nx,ny,nz,dh,dt)
      call m_wneb(xi(17),xf(17),yi(17),yf(17),zi(17),zf(17), nx,ny,nz,dh,dt)

      return
      end



      subroutine pmlstr(nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE PML MODEL IN BOUNDARY ZONES

C     XZ, YZ AND XY ARE THE SHEAR STRESS FIELDS
C     XYZ CONTAINS THE NORMAL STRESS FIELDS XX, YY AND ZZ
C     N(Y), S(Y), E(X), W(X) CODE INDICATES GEOGRAPHIC DIRECTION
C     B(Z) CODE INDICATES MODEL SPACE BOTTOM
C     SINGLE CODE IS A PLANE NORMAL TO THE DIRECTION AXIS (5 TOTAL)
C     COMBINATION OF 2 INDICATES AN EDGE (8 TOTAL)
C     COMBINATION OF 3 INDICATES A CORNER (4 TOTAL)
C     ROUTINE CALLS ORDERED PLANES, EDGES, CORNERS

C     N*        FINAL COORDINATE ADDRESSES
C     DT        TIME INCREMENT
C     DH        SPATIAL INCREMENT

C     PLANES: N, S, W, E AND B

      use parstat

      call m_xyzn(xi(3),xf(3),yi(3),yf(3),zi(3),zf(3), nx,ny,nz,dh,dt)
      call m_xzn(xi(3),xf(3),yi(3),yf(3),zi(3),zf(3), nx,ny,nz,dh,dt)
      call m_yzn(xi(3),xf(3),yi(3),yf(3),zi(3),zf(3), nx,ny,nz,dh,dt)
      call m_xyn(xi(3),xf(3),yi(3),yf(3),zi(3),zf(3), nx,ny,nz,dh,dt)

      call m_xyzs(xi(1),xf(1),yi(1),yf(1),zi(1),zf(1), nx,ny,nz,dh,dt)
      call m_xzs(xi(1),xf(1),yi(1),yf(1),zi(1),zf(1), nx,ny,nz,dh,dt)
      call m_yzs(xi(1),xf(1),yi(1),yf(1),zi(1),zf(1), nx,ny,nz,dh,dt)
      call m_xys(xi(1),xf(1),yi(1),yf(1),zi(1),zf(1), nx,ny,nz,dh,dt)

      call m_xyzw(xi(2),xf(2),yi(2),yf(2),zi(2),zf(2), nx,ny,nz,dh,dt)
      call m_xzw(xi(2),xf(2),yi(2),yf(2),zi(2),zf(2), nx,ny,nz,dh,dt)
      call m_yzw(xi(2),xf(2),yi(2),yf(2),zi(2),zf(2), nx,ny,nz,dh,dt)
      call m_xyw(xi(2),xf(2),yi(2),yf(2),zi(2),zf(2), nx,ny,nz,dh,dt)

      call m_xyze(xi(4),xf(4),yi(4),yf(4),zi(4),zf(4), nx,ny,nz,dh,dt)
      call m_xze(xi(4),xf(4),yi(4),yf(4),zi(4),zf(4), nx,ny,nz,dh,dt)
      call m_yze(xi(4),xf(4),yi(4),yf(4),zi(4),zf(4), nx,ny,nz,dh,dt)
      call m_xye(xi(4),xf(4),yi(4),yf(4),zi(4),zf(4), nx,ny,nz,dh,dt)

      call m_xyzb(xi(5),xf(5),yi(5),yf(5),zi(5),zf(5), nx,ny,nz,dh,dt)
      call m_xzb(xi(5),xf(5),yi(5),yf(5),zi(5),zf(5), nx,ny,nz,dh,dt)
      call m_yzb(xi(5),xf(5),yi(5),yf(5),zi(5),zf(5), nx,ny,nz,dh,dt)
      call m_xyb(xi(5),xf(5),yi(5),yf(5),zi(5),zf(5), nx,ny,nz,dh,dt)

C     VERTICAL EDGES: SW, SE, NW AND NE

      call m_xyzsw(xi(6),xf(6),yi(6),yf(6),zi(6),zf(6), nx,ny,nz,dh,dt)
      call m_xzsw(xi(6),xf(6),yi(6),yf(6),zi(6),zf(6), nx,ny,nz,dh,dt)
      call m_yzsw(xi(6),xf(6),yi(6),yf(6),zi(6),zf(6), nx,ny,nz,dh,dt)
      call m_xysw(xi(6),xf(6),yi(6),yf(6),zi(6),zf(6), nx,ny,nz,dh,dt)

      call m_xyzse(xi(8),xf(8),yi(8),yf(8),zi(8),zf(8), nx,ny,nz,dh,dt)
      call m_xzse(xi(8),xf(8),yi(8),yf(8),zi(8),zf(8), nx,ny,nz,dh,dt)
      call m_yzse(xi(8),xf(8),yi(8),yf(8),zi(8),zf(8), nx,ny,nz,dh,dt)
      call m_xyse(xi(8),xf(8),yi(8),yf(8),zi(8),zf(8), nx,ny,nz,dh,dt)

      call m_xyznw(xi(7),xf(7),yi(7),yf(7),zi(7),zf(7), nx,ny,nz,dh,dt)
      call m_xznw(xi(7),xf(7),yi(7),yf(7),zi(7),zf(7), nx,ny,nz,dh,dt)
      call m_yznw(xi(7),xf(7),yi(7),yf(7),zi(7),zf(7), nx,ny,nz,dh,dt)
      call m_xynw(xi(7),xf(7),yi(7),yf(7),zi(7),zf(7), nx,ny,nz,dh,dt)

      call m_xyzne(xi(9),xf(9),yi(9),yf(9),zi(9),zf(9), nx,ny,nz,dh,dt)
      call m_xzne(xi(9),xf(9),yi(9),yf(9),zi(9),zf(9), nx,ny,nz,dh,dt)
      call m_yzne(xi(9),xf(9),yi(9),yf(9),zi(9),zf(9), nx,ny,nz,dh,dt)
      call m_xyne(xi(9),xf(9),yi(9),yf(9),zi(9),zf(9), nx,ny,nz,dh,dt)

C     HORIZONTAL EDGES: SB, NB, WB AND EB

      call m_xyzsb(xi(10),xf(10),yi(10),yf(10),zi(10),zf(10), nx,ny,nz,dh,dt)
      call m_xzsb(xi(10),xf(10),yi(10),yf(10),zi(10),zf(10), nx,ny,nz,dh,dt)
      call m_yzsb(xi(10),xf(10),yi(10),yf(10),zi(10),zf(10), nx,ny,nz,dh,dt)
      call m_xysb(xi(10),xf(10),yi(10),yf(10),zi(10),zf(10), nx,ny,nz,dh,dt)

      call m_xyznb(xi(12),xf(12),yi(12),yf(12),zi(12),zf(12), nx,ny,nz,dh,dt)
      call m_xznb(xi(12),xf(12),yi(12),yf(12),zi(12),zf(12), nx,ny,nz,dh,dt)
      call m_yznb(xi(12),xf(12),yi(12),yf(12),zi(12),zf(12), nx,ny,nz,dh,dt)
      call m_xynb(xi(12),xf(12),yi(12),yf(12),zi(12),zf(12), nx,ny,nz,dh,dt)

      call m_xyzwb(xi(11),xf(11),yi(11),yf(11),zi(11),zf(11), nx,ny,nz,dh,dt)
      call m_xzwb(xi(11),xf(11),yi(11),yf(11),zi(11),zf(11), nx,ny,nz,dh,dt)
      call m_yzwb(xi(11),xf(11),yi(11),yf(11),zi(11),zf(11), nx,ny,nz,dh,dt)
      call m_xywb(xi(11),xf(11),yi(11),yf(11),zi(11),zf(11), nx,ny,nz,dh,dt)

      call m_xyzeb(xi(13),xf(13),yi(13),yf(13),zi(13),zf(13), nx,ny,nz,dh,dt)
      call m_xzeb(xi(13),xf(13),yi(13),yf(13),zi(13),zf(13), nx,ny,nz,dh,dt)
      call m_yzeb(xi(13),xf(13),yi(13),yf(13),zi(13),zf(13), nx,ny,nz,dh,dt)
      call m_xyeb(xi(13),xf(13),yi(13),yf(13),zi(13),zf(13), nx,ny,nz,dh,dt)

C     BOTTOM CORNERS: SWB, NWB, SEB AND NEB

      call m_xyzswb(xi(14),xf(14),yi(14),yf(14),zi(14),zf(14), nx,ny,nz,dh,dt)
      call m_xzswb(xi(14),xf(14),yi(14),yf(14),zi(14),zf(14), nx,ny,nz,dh,dt)
      call m_yzswb(xi(14),xf(14),yi(14),yf(14),zi(14),zf(14), nx,ny,nz,dh,dt)
      call m_xyswb(xi(14),xf(14),yi(14),yf(14),zi(14),zf(14), nx,ny,nz,dh,dt)

      call m_xyznwb(xi(15),xf(15),yi(15),yf(15),zi(15),zf(15), nx,ny,nz,dh,dt)
      call m_xznwb(xi(15),xf(15),yi(15),yf(15),zi(15),zf(15), nx,ny,nz,dh,dt)
      call m_yznwb(xi(15),xf(15),yi(15),yf(15),zi(15),zf(15), nx,ny,nz,dh,dt)
      call m_xynwb(xi(15),xf(15),yi(15),yf(15),zi(15),zf(15), nx,ny,nz,dh,dt)

      call m_xyzseb(xi(16),xf(16),yi(16),yf(16),zi(16),zf(16), nx,ny,nz,dh,dt)
      call m_xzseb(xi(16),xf(16),yi(16),yf(16),zi(16),zf(16), nx,ny,nz,dh,dt)
      call m_yzseb(xi(16),xf(16),yi(16),yf(16),zi(16),zf(16), nx,ny,nz,dh,dt)
      call m_xyseb(xi(16),xf(16),yi(16),yf(16),zi(16),zf(16), nx,ny,nz,dh,dt)

      call m_xyzneb(xi(17),xf(17),yi(17),yf(17),zi(17),zf(17), nx,ny,nz,dh,dt)
      call m_xzneb(xi(17),xf(17),yi(17),yf(17),zi(17),zf(17), nx,ny,nz,dh,dt)
      call m_yzneb(xi(17),xf(17),yi(17),yf(17),zi(17),zf(17), nx,ny,nz,dh,dt)
      call m_xyneb(xi(17),xf(17),yi(17),yf(17),zi(17),zf(17), nx,ny,nz,dh,dt)

      return
      end


C VELOCITY PLANES =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

                                                                                 
      subroutine m_un(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE U IN NORTH PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      unT(i,ny-j+1,k) = (2.*dth*d*(xy(i,j,k) - xy(i,j-1,k)) +
     +                unT(i,ny-j+1,k)*
     +                (2. - dt*rf(ny-j+1)))/(2. + dt*rf(ny-j+1))

      
      unH(i,ny-j+1,k) = (2.*dth*d* ((xx(i,j,k) - xx(i-1,j,k)) + 
     +                  (xz(i,j,k) - xz(i,j,k-1))) +
     +                  unH(i,ny-j+1,k)* 
     +                  (2. - dt * pht * rf(ny-j+1)))/
     +                  (2. + dt * pht * rf(ny-j+1))

      u1(i,j,k) = unT(i,ny-j+1,k) + unH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_vn(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE V IN NORTH PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      
      vnT(i,ny-j+1,k) = (2.*dth*d*(yy(i,j+1,k) - yy(i,j,k)) +
     +                vnT(i,ny-j+1,k)*
     +                (2. - dt*rh(ny-j+1)))/(2. + dt*rh(ny-j+1))

      vnH(i,ny-j+1,k) = (2.*dth*d*((xy(i+1,j,k) - xy(i,j,k)) +
     +                (yz(i,j,k) - yz(i,j,k-1))) +
     +                vnH(i,ny-j+1,k) *
     +                (2. - dt*pht*rh(ny-j+1)))/
     +                (2. + dt*pht*rh(ny-j+1))

      v1(i,j,k) = vnT(i,ny-j+1,k) + vnH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_wn(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W IN NORTH PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)
      

      wnT(i,ny-j+1,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j-1,k)) +
     +                wnT(i,ny-j+1,k)*
     +                (2. - dt*rf(ny-j+1)))/(2. + dt*rf(ny-j+1))

      wnH(i,ny-j+1,k) = (2.*dth*d*((xz(i+1,j,k) - xz(i,j,k)) +
     +                (zz(i,j,k+1) - zz(i,j,k))) +
     +                wnH(i,ny-j+1,k) * 
     +                (2. - dt * pht*rf(ny-j+1))) / 
     +                (2. + dt * pht*rf(ny-j+1)) 

      w1(i,j,k) = wnT(i,ny-j+1,k) + wnH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_us(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE U IN SOUTH PLANE

      use parstat
 
      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      usT(i,j,k) = (2.*dth*d*(xy(i,j,k) - xy(i,j-1,k)) +
     +             usT(i,j,k)*
     +             (2. - dt*rh(j)))/(2. + dt*rh(j))

      usH(i,j,k) = (2.*d*dth*((xx(i,j,k) - xx(i-1,j,k)) +
     +             (xz(i,j,k) - xz(i,j,k-1)))+
     +             usH(i,j,k)*
     +             (2. - dt*pht*rh(j)))/(2. + dt*pht*rh(j))
                 
      u1(i,j,k) = usT(i,j,k) + usH(i,j,k)
 
   50 continue
      enddo
      enddo
      return
      end



      subroutine m_vs(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE V IN SOUTH PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      vsT(i,j,k) = (2.*dth*d*(yy(i,j+1,k) - yy(i,j,k)) +
     +             vsT(i,j,k)*
     +             (2. - dt*rf(j)))/(2. + dt*rf(j))

      vsH(i,j,k) = (2.*d*dth*((xy(i+1,j,k) - xy(i,j,k)) +
     +             (yz(i,j,k) - yz(i,j,k-1))) +
     +             vsH(i,j,k)*
     +             (2. - dt*pht*rf(j)))/(2. + dt*pht*rf(j))

      v1(i,j,k) = vsT(i,j,k) + vsH(i,j,k)

   50 continue 
      enddo
      enddo
      return
      end



      subroutine m_ws(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W IN SOUTH PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      wsT(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j-1,k)) +
     +             wsT(i,j,k)*
     +             (2. - dt*rh(j)))/(2. + dt*rh(j))

      wsH(i,j,k) = (2.*d*dth*((xz(i+1,j,k) - xz(i,j,k)) +
     +             (zz(i,j,k+1) - zz(i,j,k)))+
     +             wsH(i,j,k) *
     +             (2. - dt*pht*rh(j)))/(2. + dt*pht*rh(j))

      w1(i,j,k) = wsT(i,j,k) + wsH(i,j,k)
 
   50 continue
      enddo
      enddo
      return
      end



      subroutine m_uw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE U IN WEST PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      uwT(i,j,k) = (2.*dth*d*(xx(i,j,k) - xx(i-1,j,k)) +
     +             uwT(i,j,k)*
     +             (2. - dt*rh(i)))/(2. + dt*rh(i))

      uwH(i,j,k) = (2.*d*dth*((xy(i,j,k) - xy(i,j-1,k)) +
     +             (xz(i,j,k) - xz(i,j,k-1)))+
     +             uwH(i,j,k)*
     +             (2. - dt*pht*rh(i)))/(2. + dt*pht*rh(i))

      u1(i,j,k) = uwT(i,j,k) + uwH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_vw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE V IN WEST PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      vwT(i,j,k) = (2.*dth*d*(xy(i+1,j,k) - xy(i,j,k)) +
     +             vwT(i,j,k)*
     +             (2. - dt*rf(i)))/(2. + dt*rf(i))

      vwH(i,j,k) = (2.*d*dth*((yy(i,j+1,k) - yy(i,j,k)) +
     +             (yz(i,j,k) - yz(i,j,k-1))) +
     +             vwH(i,j,k) * 
     +             (2. - dt*pht*rf(i)))/(2. + dt*pht*rf(i))

      v1(i,j,k) = vwT(i,j,k) + vwH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_ww(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W IN WEST PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      wwT(i,j,k) = (2.*dth*d*(xz(i+1,j,k) - xz(i,j,k)) +
     +             wwT(i,j,k)*
     +             (2. - dt*rf(i)))/(2. + dt*rf(i))

      wwH(i,j,k) = (2.*d*dth*((yz(i,j,k) - yz(i,j-1,k)) +
     +             (zz(i,j,k+1) - zz(i,j,k))) + 
     +             wwH(i,j,k)*
     +             (2. - dt*pht*rf(i)))/(2. + dt*pht*rf(i))

      w1(i,j,k) = wwT(i,j,k) + wwH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_ue(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE U IN EAST PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      ueT(nx-i+1,j,k) = (2.*dth*d*(xx(i,j,k) - xx(i-1,j,k)) +
     +                ueT(nx-i+1,j,k)*
     +                (2. - dt*rf(nx-i+1)))/(2. + dt*rf(nx-i+1))

      ueH(nx-i+1,j,k) = (2.*d*dth*((xy(i,j,k) - xy(i,j-1,k)) +
     +                (xz(i,j,k) - xz(i,j,k-1))) +
     +                ueH(nx-i+1,j,k)*
     +                (2. - dt*pht*rf(nx-i+1)))/(2. + dt*pht*rf(nx-i+1))

      u1(i,j,k) = ueT(nx-i+1,j,k) + ueH(nx-i+1,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_ve(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE V IN EAST PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      veT(nx-i+1,j,k) = (2.*dth*d*(xy(i+1,j,k) - xy(i,j,k)) +
     +                veT(nx-i+1,j,k)*
     +                (2. - dt*rh(nx-i+1)))/(2. + dt*rh(nx-i+1))

      veH(nx-i+1,j,k) = (2.*d*dth*(
     +                (yy(i,j+1,k) - yy(i,j,k)) +
     +                (yz(i,j,k) - yz(i,j,k-1))) +
     +                veH(nx-i+1,j,k) *
     +                (2. - dt*pht*rh(nx-i+1)))/(2. + dt*pht*rh(nx-i+1))

      v1(i,j,k) = veT(nx-i+1,j,k) + veH(nx-i+1,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_we(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W IN EAST PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      weT(nx-i+1,j,k) = (2.*dth*d*(xz(i+1,j,k) - xz(i,j,k)) +
     +                weT(nx-i+1,j,k)*
     +                (2. - dt*rh(nx-i+1)))/(2. + dt*rh(nx-i+1))

      weH(nx-i+1,j,k) = (2.*d*dth*(
     +                (yz(i,j,k) - yz(i,j-1,k)) +
     +                (zz(i,j,k+1) - zz(i,j,k))) +
     +                weH(nx-i+1,j,k) *
     +                (2. - dt*pht*rh(nx-i+1)))/(2. + dt*pht*rh(nx-i+1))

      w1(i,j,k) = weT(nx-i+1,j,k) + weH(nx-i+1,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_ub(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE U IN BOTTOM PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      ubT(i,j,k) = (2.*dth*d*(xz(i,j,k) - xz(i,j,k-1)) +
     +             ubT(i,j,k)*
     +             (2. - dt*rh(k)))/(2. + dt*rh(k))

      ubH(i,j,k) = (2.*d*dth*(
     +             (xy(i,j,k) - xy(i,j-1,k)) +
     +             (xx(i,j,k) - xx(i-1,j,k))) +
     +             ubH(i,j,k)*
     +             (2. - dt*pht*rh(k)))/(2. + dt*pht*rh(k))

      u1(i,j,k) = ubT(i,j,k) + ubH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_vb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE V IN BOTTOM PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      vbT(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j,k-1)) +
     +             vbT(i,j,k)*
     +             (2. - dt*rh(k)))/(2. + dt*rh(k))

      vbH(i,j,k) = (2.*d*dth*(
     +             (xy(i+1,j,k) - xy(i,j,k)) +
     +             (yy(i,j+1,k) - yy(i,j,k))) +
     +             vbH(i,j,k)*
     +             (2. - dt*pht*rh(k)))/(2. + dt*pht*rh(k))

      v1(i,j,k) = vbT(i,j,k) + vbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_wb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W IN BOTTOM PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      wbT(i,j,k) = (2.*dth*d*(zz(i,j,k+1) - zz(i,j,k)) +
     +             wbT(i,j,k)*
     +             (2. - dt*rf(k)))/(2. + dt*rf(k))

      wbH(i,j,k) =  (2.*d*dth*(
     +             (yz(i,j,k) - yz(i,j-1,k)) +
     +             (xz(i+1,j,k) - xz(i,j,k))) +
     +             wbH(i,j,k)*
     +             (2. - dt*pht*rf(k)))/(2. + dt*pht*rf(k))

      w1(i,j,k) = wbT(i,j,k) + wbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end


C VELOCITY PLANES =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
C VERTICAL VELOCITY EDGES =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


      subroutine m_usw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE U AT THE SOUTH-WEST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rh(j) + pht*rh(i)
      usT(i,j,k) = (2.*dth*d*(xy(i,j,k) - xy(i,j-1,k)) +
     +             usT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rh(i)
      uwT(i,j,k) = (2.*dth*d*(xx(i,j,k) - xx(i-1,j,k)) +
     +             uwT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rh(i))
      usH(i,j,k) = (2.*dth*d*(xz(i,j,k) - xz(i,j,k-1)) +
     +             usH(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      u1(i,j,k) = usT(i,j,k) + uwT(i,j,k) + usH(i,j,k) 

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_vsw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE V AT THE SOUTH-WEST EDGE  

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rf(j) + pht*rf(i)
      vsT(i,j,k) = (2.*dth*d*(yy(i,j+1,k) - yy(i,j,k)) +
     +             vsT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + rf(i)
      vwT(i,j,k) = (2.*dth*d*(xy(i+1,j,k) - xy(i,j,k)) +
     +             vwT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(j) + rf(i))
      vsH(i,j,k) =  (2.*dth*d*(yz(i,j,k) - yz(i,j,k-1)) +
     +              vsH(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      v1(i,j,k) = vsT(i,j,k) + vwT(i,j,k) + vsH(i,j,k)

   50 continue
      enddo
      enddo
      return  
      end     



      subroutine m_wsw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W AT THE SOUTH-WEST EDGE   

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rh(j) + pht*rf(i)
      wsT(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j-1,k)) +
     +             wsT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rf(i)
      wwT(i,j,k) = (2.*dth*d*(xz(i+1,j,k) - xz(i,j,k)) +
     +             wwT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rf(i))
      wsH(i,j,k) = (2.*dth*d*(zz(i,j,k+1) - zz(i,j,k)) +
     +             wsH(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      w1(i,j,k) = wsT(i,j,k) + wwT(i,j,k) + wsH(i,j,k)

   50 continue
      enddo
      enddo
      return  
      end     



      subroutine m_use(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE U AT THE SOUTH-EAST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rh(j) + pht*rf(nx-i+1)
      usT(i,j,k) = (2.*dth*d*(xy(i,j,k) - xy(i,j-1,k)) +
     +             usT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rf(nx-i+1)
      ueT(nx-i+1,j,k) = (2.*dth*d*(xx(i,j,k) - xx(i-1,j,k)) +
     +                ueT(nx-i+1,j,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rf(nx-i+1))
      usH(i,j,k) = (2.*dth*d*(xz(i,j,k) - xz(i,j,k-1)) +
     +             usH(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      u1(i,j,k) = usT(i,j,k) + ueT(nx-i+1,j,k) + usH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_vse(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE V AT THE SOUTH-EAST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rf(j) + pht*rh(nx-i+1)
      vsT(i,j,k) = (2.*dth*d*(yy(i,j+1,k) - yy(i,j,k)) +
     +             vsT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + rh(nx-i+1)
      veT(nx-i+1,j,k) = (2.*dth*d*(xy(i+1,j,k) - xy(i,j,k)) +
     +                veT(nx-i+1,j,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(j) + rh(nx-i+1))
      vsH(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j,k-1))+
     +             vsH(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      v1(i,j,k) = vsT(i,j,k) + veT(nx-i+1,j,k) + vsH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_wse(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W AT THE SOUTH-EAST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rh(j) + pht * rh(nx-i+1)
      wsT(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j-1,k)) +
     +             wsT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rh(nx-i+1)
      weT(nx-i+1,j,k) = (2.*dth*d*(xz(i+1,j,k) - xz(i,j,k)) +
     +                weT(nx-i+1,j,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rh(nx-i+1))
      wsH(i,j,k) = (2.*dth*d*(zz(i,j,k+1) - zz(i,j,k)) +
     +             wsH(i,j,k) *
     +             (2. - dt*dc))/(2. + dt*dc)

      w1(i,j,k) = wsT(i,j,k) + weT(nx-i+1,j,k) + wsH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_unw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE U AT THE NORTH-WEST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc = rf(ny-j+1) + pht*rh(i)
      unT(i,ny-j+1,k) = (2.*dth*d*(xy(i,j,k) - xy(i,j-1,k)) +
     +                unT(i,ny-j+1,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc = pht*rf(ny-j+1) + rh(i)
      uwT(i,j,k) = (2.*dth*d*(xx(i,j,k) - xx(i-1,j,k)) +
     +             uwT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc = pht*(rf(ny-j+1) + rh(i))
      unH(i,ny-j+1,k) =(2.*dth*d*(xz(i,j,k) - xz(i,j,k-1)) +
     +                unH(i,ny-j+1,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      u1(i,j,k) = unT(i,ny-j+1,k) + uwT(i,j,k) + unH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_vnw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE V AT THE NORTH-WEST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc = rh(ny-j+1) + pht * rf(i)
      vnT(i,ny-j+1,k) = (2.*dth*d*(yy(i,j+1,k) - yy(i,j,k)) +
     +                vnT(i,ny-j+1,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc = pht * rh(ny-j+1) + rf(i)
      vwT(i,j,k) = (2.*dth*d*(xy(i+1,j,k) - xy(i,j,k)) +
     +             vwT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc = pht * (rh(ny-j+1) + rf(i))
      vnH(i,ny-j+1,k) = (2.* dth*d*(yz(i,j,k) - yz(i,j,k-1)) +
     +                vnH(i,ny-j+1,k) *
     +             (2. - dt*dc))/(2. + dt*dc)

      v1(i,j,k) = vnT(i,ny-j+1,k) + vwT(i,j,k) + vnH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_wnw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W AT THE NORTH-WEST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rf(ny-j+1) + pht * rf(i)
      wnT(i,ny-j+1,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j-1,k)) +
     +                wnT(i,ny-j+1,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht * rf(ny-j+1) + rf(i)
      wwT(i,j,k) = (2.*dth*d*(xz(i+1,j,k) - xz(i,j,k)) +
     +             wwT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc = pht * (rf(ny-j+1) + rf(i))
      wnH(i,ny-j+1,k) = (2.*dth*d*(zz(i,j,k+1) - zz(i,j,k)) +
     +              wnH(i,ny-j+1,k) *
     +             (2. - dt*dc))/(2. + dt*dc)

      w1(i,j,k) = wnT(i,ny-j+1,k) + wwT(i,j,k) + wnH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_une(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE U AT THE NORTH-EAST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rf(ny-j+1) + pht*rf(nx-i+1)
      unT(i,ny-j+1,k) = (2.*dth*d*(xy(i,j,k) - xy(i,j-1,k)) +
     +                unT(i,ny-j+1,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rf(nx-i+1)
      ueT(nx-i+1,j,k) = (2.*dth*d*(xx(i,j,k) - xx(i-1,j,k)) +
     +                ueT(nx-i+1,j,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(ny-j+1) + rf(nx-i+1))
      unH(i,ny-j+1,k) = (2.*dth*d*(xz(i,j,k) - xz(i,j,k-1))+
     +                unH(i,ny-j+1,k)* 
     +                (2. - dt*dc))/(2. + dt*dc)

      u1(i,j,k) = unT(i,ny-j+1,k) + ueT(nx-i+1,j,k) + unH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_vne(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE V AT THE NORTH-EAST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rh(ny-j+1) + pht*rh(nx-i+1)
      vnT(i,ny-j+1,k) = (2.*dth*d*(yy(i,j+1,k) - yy(i,j,k)) +
     +                vnT(i,ny-j+1,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + rh(nx-i+1)
      veT(nx-i+1,j,k) = (2.*dth*d*(xy(i+1,j,k) - xy(i,j,k)) +
     +                veT(nx-i+1,j,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(ny-j+1) + rh(nx-i+1))
      vnH(i,ny-j+1,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j,k-1)) +
     +                vnH(i,ny-j+1,k) *
     +                (2. - dt*dc))/(2. + dt*dc)

      v1(i,j,k) = vnT(i,ny-j+1,k) + veT(nx-i+1,j,k) + vnH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_wne(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W AT THE NORTH-EAST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rf(ny-j+1) + pht * rh(nx-i+1) 
      wnT(i,ny-j+1,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j-1,k)) +
     +                wnT(i,ny-j+1,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rh(nx-i+1) 
      weT(nx-i+1,j,k) = (2.*dth*d*(xz(i+1,j,k) - xz(i,j,k)) +
     +                weT(nx-i+1,j,k) *
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(ny-j+1) + rh(nx-i+1)) 
      wnH(i,ny-j+1,k) = (2.*dth*d*(zz(i,j,k+1) - zz(i,j,k)) +
     +                wnH(i,ny-j+1,k) *
     +                (2. - dt*dc))/(2. + dt*dc)

      w1(i,j,k) = wnT(i,ny-j+1,k) + weT(nx-i+1,j,k) + wnH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



C VERTICAL VELOCITY EDGES =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
C HORIZONTAL VELOCITY EDGES =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


      subroutine m_usb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE U AT THE SOUTH-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rh(j) + pht*rh(k)
      usT(i,j,k) = (2.*dth*d*(xy(i,j,k) - xy(i,j-1,k)) +
     +             usT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rh(k)
      ubT(i,j,k) = (2.*dth*d*(xz(i,j,k) - xz(i,j,k-1)) +
     +             ubT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rh(k))
      ubH(i,j,k) = (2.*dth*d*(xx(i,j,k) - xx(i-1,j,k)) +
     +             ubH(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      u1(i,j,k) = usT(i,j,k) + ubT(i,j,k) + ubH(i,j,k)

   50 continue
      enddo
      enddo
      return  
      end     



      subroutine m_vsb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE V AT THE SOUTH-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rf(j) + pht*rh(k)
      vsT(i,j,k) = (2.*dth*d*(yy(i,j+1,k) - yy(i,j,k)) +
     +             vsT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + rh(k)
      vbT(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j,k-1)) +
     +             vbT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(j) + rh(k))
      vbH(i,j,k) = (2.*dth*d*(xy(i+1,j,k) - xy(i,j,k)) +
     +             vbH(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      v1(i,j,k) = vsT(i,j,k) + vbT(i,j,k) + vbH(i,j,k)

   50 continue
      enddo
      enddo
      return  
      end     



      subroutine m_wsb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W AT THE SOUTH-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rh(j) + pht*rf(k)
      wsT(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j-1,k)) +
     +             wsT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rf(k)
      wbT(i,j,k) = (2.*dth*d*(zz(i,j,k+1) - zz(i,j,k)) +
     +             wbT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rf(k))
      wbH(i,j,k) = (2.*dth*d*(xz(i+1,j,k) - xz(i,j,k)) +
     +             wbH(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      w1(i,j,k) = wsT(i,j,k) + wbT(i,j,k) + wbH(i,j,k)

   50 continue
      enddo
      enddo
      return  
      end     



      subroutine m_unb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE U AT THE NORTH-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rf(ny-j+1) + pht*rh(k)
      unT(i,ny-j+1,k) = (2.*dth*d*(xy(i,j,k) - xy(i,j-1,k)) +
     +                unT(i,ny-j+1,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rh(k)
      ubT(i,j,k) = (2.*dth*d*(xz(i,j,k) - xz(i,j,k-1)) +
     +             ubT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(ny-j+1) + rh(k))
      ubH(i,j,k) = (2.*dth*d*(xx(i,j,k) - xx(i-1,j,k)) +
     +             ubH(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      u1(i,j,k) = unT(i,ny-j+1,k) + ubT(i,j,k) + ubH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_vnb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE V AT THE NORTH-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rh(ny-j+1) + pht*rh(k)
      vnT(i,ny-j+1,k) = (2.*dth*d*(yy(i,j+1,k) - yy(i,j,k)) +
     +                vnT(i,ny-j+1,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + rh(k)
      vbT(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j,k-1)) +
     +             vbT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(ny-j+1) + rh(k))
      vbH(i,j,k) = (2.*dth*d*(xy(i+1,j,k) - xy(i,j,k)) +
     +             vbH(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      v1(i,j,k) = vnT(i,ny-j+1,k) + vbT(i,j,k) + vbH(i,j,k)

   50 continue
      enddo
      enddo
      return 
      end    



      subroutine m_wnb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W AT THE NORTH-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rf(ny-j+1) + pht*rf(k)
      wnT(i,ny-j+1,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j-1,k)) +
     +                wnT(i,ny-j+1,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rf(k)
      wbT(i,j,k) = (2.*dth*d*(zz(i,j,k+1) - zz(i,j,k)) +
     +             wbT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(ny-j+1) + rf(k))
      wbH(i,j,k) = (2.*dth*d*(xz(i+1,j,k) - xz(i,j,k)) +
     +             wbH(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      w1(i,j,k) = wnT(i,ny-j+1,k) + wbT(i,j,k) + wbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_uwb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W AT THE WEST-BOTTOM EDGE

      use parstat

      dth = dt/dh 

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)
      
      dc=rh(i) + pht*rh(k)
      uwT(i,j,k) = (2.*dth*d*(xx(i,j,k) - xx(i-1,j,k)) +
     +             uwT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(i) + rh(k)
      ubT(i,j,k) = (2.*dth*d*(xz(i,j,k) - xz(i,j,k-1)) +
     +             ubT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)
 
      dc=pht*(rh(i) + rh(k))
      ubH(i,j,k) = (2.*dth*d*(xy(i,j,k) - xy(i,j-1,k)) +
     +             ubH(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      u1(i,j,k) = uwT(i,j,k) + ubT(i,j,k) + ubH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_vwb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W AT THE WEST-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)
     
      dc=rf(i) + pht*rh(k)
      vwT(i,j,k) = (2.*dth*d*(xy(i+1,j,k) - xy(i,j,k)) +
     +             vwT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(i) + rh(k)
      vbT(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j,k-1)) +
     +             vbT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(i) + rh(k))
      vbH(i,j,k) = (2.*dth*d*(yy(i,j+1,k) - yy(i,j,k)) +
     +             vbH(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      v1(i,j,k) = vwT(i,j,k) + vbT(i,j,k) + vbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_wwb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W AT THE WEST-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rf(i) + pht*rf(k)
      wwT(i,j,k) = (2.*dth*d*(xz(i+1,j,k) - xz(i,j,k)) +
     +             wwT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(i) + rf(k)
      wbT(i,j,k) = (2.*dth*d*(zz(i,j,k+1) - zz(i,j,k)) +
     +             wbT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(i) + rf(k))
      wbH(i,j,k) = (2.* dth*d*(yz(i,j,k) - yz(i,j-1,k)) +
     +              wbH(i,j,k) *
     +             (2. - dt*dc))/(2. + dt*dc)

      w1(i,j,k) = wwT(i,j,k) + wbT(i,j,k) + wbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_ueb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W AT THE EAST-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rf(nx-i+1) + pht*rh(k)
      ueT(nx-i+1,j,k) = (2.*dth*d*(xx(i,j,k) - xx(i-1,j,k)) +
     +                ueT(nx-i+1,j,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(nx-i+1) + rh(k)
      ubT(i,j,k) = (2.*dth*d*(xz(i,j,k) - xz(i,j,k-1)) +
     +             ubT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(nx-i+1) + rh(k))
      ubH(i,j,k) = (2.*dth*d*(xy(i,j,k) - xy(i,j-1,k))+
     +             ubH(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      u1(i,j,k) = ueT(nx-i+1,j,k) + ubT(i,j,k) + ubH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_veb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W AT THE EAST-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rh(nx-i+1) + pht*rh(k)
      veT(nx-i+1,j,k) = (2.*dth*d*(xy(i+1,j,k) - xy(i,j,k)) +
     +                veT(nx-i+1,j,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(nx-i+1) + rh(k)
      vbT(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j,k-1)) +
     +             vbT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(nx-i+1) + rh(k))
      vbH(i,j,k) = (2.*dth*d*(yy(i,j+1,k) - yy(i,j,k)) +
     +             vbH(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      v1(i,j,k) = veT(nx-i+1,j,k) + vbT(i,j,k) + vbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_web(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W AT THE EAST-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rh(nx-i+1) + pht*rf(k)
      weT(nx-i+1,j,k) = (2.*dth*d*(xz(i+1,j,k) - xz(i,j,k)) +
     +                weT(nx-i+1,j,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(nx-i+1) + rf(k)
      wbT(i,j,k) = (2.*dth*d*(zz(i,j,k+1) - zz(i,j,k)) +
     +             wbT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(nx-i+1) + rf(k))
      wbH(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j-1,k)) +
     +             wbH(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      w1(i,j,k) = weT(nx-i+1,j,k) + wbT(i,j,k) + wbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end


C HORIZONTAL VELOCITY EDGES -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
C VELOCITY CORNERS =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


      subroutine m_uswb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE U AT THE SOUTH-WEST-BOTTOM CORNER 

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)
      
      dc=rh(j) + pht*rh(i) + pht*rh(k)
      usT(i,j,k) = (2.*dth*d*(xy(i,j,k) - xy(i,j-1,k)) +
     +             usT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rh(i) + pht*rh(k)
      uwT(i,j,k) = (2.*dth*d*(xx(i,j,k) - xx(i-1,j,k)) +
     +             uwT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + pht*rh(i) + rh(k)
      ubT(i,j,k) = (2.*dth*d*(xz(i,j,k) - xz(i,j,k-1)) +
     +             ubT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      u1(i,j,k) = usT(i,j,k) + uwT(i,j,k) + ubT(i,j,k)

   50 continue
      enddo
      enddo
      return  
      end     



      subroutine m_vswb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE V AT THE SOUTH-WEST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rf(j) + pht*rf(i) + pht*rh(k)
      vsT(i,j,k) = (2.*dth*d*(yy(i,j+1,k) - yy(i,j,k)) +
     +             vsT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + rf(i) + pht*rh(k)
      vwT(i,j,k) = (2.*dth*d*(xy(i+1,j,k) - xy(i,j,k)) +
     +             vwT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + pht*rf(i) + rh(k)
      vbT(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j,k-1)) +
     +             vbT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      v1(i,j,k) = vsT(i,j,k) + vwT(i,j,k) + vbT(i,j,k)

   50 continue
      enddo
      enddo
      return  
      end     



      subroutine m_wswb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W AT THE SOUTH-WEST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)
  
      dc=rh(j) + pht*rf(i) + pht*rf(k)
      wsT(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j-1,k)) +
     +             wsT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rf(i) + pht*rf(k)
      wwT(i,j,k) = (2.*dth*d*(xz(i+1,j,k) - xz(i,j,k)) +
     +             wwT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + pht*rf(i) + rf(k)
      wbT(i,j,k) = (2.*dth*d*(zz(i,j,k+1) - zz(i,j,k)) +
     +             wbT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      w1(i,j,k) = wsT(i,j,k) + wwT(i,j,k) + wbT(i,j,k)

   50 continue
      enddo
      enddo
      return  
      end     



      subroutine m_unwb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE U AT THE NORTH-WEST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rf(ny-j+1) + pht*rh(i) + pht*rh(k)
      unT(i,ny-j+1,k) = (2.*dth*d*(xy(i,j,k) - xy(i,j-1,k)) +
     +                unT(i,ny-j+1,k)*
     +                (2. - dt*dc))/(2. + dt*dc)
 
      dc=pht*rf(ny-j+1) + rh(i) + pht*rh(k)
      uwT(i,j,k) = (2.*dth*d*(xx(i,j,k) - xx(i-1,j,k)) +
     +             uwT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + pht*rh(i) + rh(k)
      ubT(i,j,k) = (2.*dth*d*(xz(i,j,k) - xz(i,j,k-1)) +
     +             ubT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      u1(i,j,k) = unT(i,ny-j+1,k) + uwT(i,j,k) + ubT(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_vnwb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE V AT THE NORTH-WEST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rh(ny-j+1) + pht*rf(i) + pht*rh(k)
      vnT(i,ny-j+1,k) = (2.*dth*d*(yy(i,j+1,k) - yy(i,j,k)) +
     +                vnT(i,ny-j+1,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + rf(i) + pht*rh(k)
      vwT(i,j,k) = (2.*dth*d*(xy(i+1,j,k) - xy(i,j,k)) +
     +             vwT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)
      
      dc=pht*rh(ny-j+1) + pht*rf(i) + rh(k)
      vbT(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j,k-1)) +
     +             vbT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      v1(i,j,k) = vnT(i,ny-j+1,k) + vwT(i,j,k) + vbT(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_wnwb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W AT THE NORTH-WEST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rf(ny-j+1) + pht*rf(i) + pht*rf(k)
      wnT(i,ny-j+1,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j-1,k)) +
     +                wnT(i,ny-j+1,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rf(i) + pht*rf(k)
      wwT(i,j,k) = (2.*dth*d*(xz(i+1,j,k) - xz(i,j,k)) +
     +             wwT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + pht*rf(i) + rf(k)
      wbT(i,j,k) = (2.*dth*d*(zz(i,j,k+1) - zz(i,j,k)) +
     +             wbT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      w1(i,j,k) = wnT(i,ny-j+1,k) + wwT(i,j,k) + wbT(i,j,k) 

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_useb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE U AT THE SOUTH-EAST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rh(j) + pht*rf(nx-i+1) + pht*rh(k)
      usT(i,j,k) = (2.*dth*d*(xy(i,j,k) - xy(i,j-1,k)) +
     +             usT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rf(nx-i+1) + pht*rh(k)
      ueT(nx-i+1,j,k) = (2.*dth*d*(xx(i,j,k) - xx(i-1,j,k)) +
     +                ueT(nx-i+1,j,k)*
     +                (2. - dt*dc))/(2. + dt*dc)
    
      dc=pht*rh(j) + pht*rf(nx-i+1) + rh(k)
      ubT(i,j,k) = (2.*dth*d*(xz(i,j,k) - xz(i,j,k-1)) +
     +             ubT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)
 
      u1(i,j,k) = usT(i,j,k) + ueT(nx-i+1,j,k) + ubT(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_vseb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE V AT THE SOUTH-EAST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rf(j) + pht*rh(nx-i+1) + pht*rh(k)
      vsT(i,j,k) = (2.*dth*d*(yy(i,j+1,k) - yy(i,j,k)) +
     +             vsT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + rh(nx-i+1) + pht*rh(k)
      veT(nx-i+1,j,k) = (2.*dth*d*(xy(i+1,j,k) - xy(i,j,k)) +
     +                veT(nx-i+1,j,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + pht*rh(nx-i+1) + rh(k)
      vbT(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j,k-1)) +
     +             vbT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)
 
      v1(i,j,k) = vsT(i,j,k) + veT(nx-i+1,j,k) + vbT(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_wseb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W AT THE SOUTH-EAST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rh(j) + pht*rh(nx-i+1) + pht*rf(k)
      wsT(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j-1,k)) +
     +             wsT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rh(nx-i+1) + pht*rf(k)
      weT(nx-i+1,j,k) = (2.*dth*d*(xz(i+1,j,k) - xz(i,j,k)) +
     +                weT(nx-i+1,j,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + pht*rh(nx-i+1) + rf(k)
      wbT(i,j,k) = (2.*dth*d*(zz(i,j,k+1) - zz(i,j,k)) +
     +             wbT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)
           
      w1(i,j,k) = wsT(i,j,k) + weT(nx-i+1,j,k) + wbT(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_uneb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE U AT THE NORTH-EAST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rf(ny-j+1) + pht*rf(nx-i+1) + pht*rh(k)
      unT(i,ny-j+1,k) = (2.*dth*d*(xy(i,j,k) - xy(i,j-1,k)) +
     +                unT(i,ny-j+1,k)*
     +                (2. - dt*dc))/(2. + dt*dc)
      
      dc=pht*rf(ny-j+1) + rf(nx-i+1) + pht*rh(k)
      ueT(nx-i+1,j,k) = (2.*dth*d*(xx(i,j,k) - xx(i-1,j,k)) +
     +                ueT(nx-i+1,j,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + pht*rf(nx-i+1) + rh(k)
      ubT(i,j,k) = (2.*dth*d*(xz(i,j,k) - xz(i,j,k-1)) +
     +             ubT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      u1(i,j,k) = unT(i,ny-j+1,k) + ueT(nx-i+1,j,k) + ubT(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_vneb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE V AT THE NORTH-EAST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rh(ny-j+1) + pht*rh(nx-i+1) + pht*rh(k)
      vnT(i,ny-j+1,k) = (2.*dth*d*(yy(i,j+1,k) - yy(i,j,k)) +
     +                vnT(i,ny-j+1,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + rh(nx-i+1) + pht*rh(k)
      veT(nx-i+1,j,k) = (2.*dth*d*(xy(i+1,j,k) - xy(i,j,k)) +
     +                veT(nx-i+1,j,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + pht*rh(nx-i+1) + rh(k)
      vbT(i,j,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j,k-1)) +
     +             vbT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      v1(i,j,k) = vnT(i,ny-j+1,k) + veT(nx-i+1,j,k) + vbT(i,j,k) 

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_wneb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE W AT THE NORTH-EAST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      d = 1./d1(i,j,k)

      dc=rf(ny-j+1) + pht*rh(nx-i+1) + pht*rf(k)
      wnT(i,ny-j+1,k) = (2.*dth*d*(yz(i,j,k) - yz(i,j-1,k)) +
     +                wnT(i,ny-j+1,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rh(nx-i+1) + pht*rf(k)
      weT(nx-i+1,j,k) = (2.*dth*d*(xz(i+1,j,k) - xz(i,j,k)) +
     +                weT(nx-i+1,j,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + pht*rh(nx-i+1) + rf(k)
      wbT(i,j,k) = (2.*dth*d*(zz(i,j,k+1) - zz(i,j,k)) +
     +             wbT(i,j,k)*
     +             (2. - dt*dc))/(2. + dt*dc)

      w1(i,j,k) = wnT(i,ny-j+1,k) + weT(nx-i+1,j,k) + wbT(i,j,k)

   50 continue
      enddo
      enddo
      return
      end


C VELOCITY CORNERS =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C STRESS PLANES -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


      subroutine m_xyzn(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ IN NORTH PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl
      

      xxnT(i,ny-j+1,k) = (2.*dth*xl*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 xxnT(i,ny-j+1,k)*
     +                 (2. - dt*rf(ny-j+1)))/(2. + dt*rf(ny-j+1))

      xxnH(i,ny-j+1,k) = (2.*dth* (xl2m*(u1(i+1,j,k) - u1(i,j,k)) +
     +                 xl*(w1(i,j,k) - w1(i,j,k-1))) +
     +                 xxnH(i,ny-j+1,k) * 
     +                 (2. - dt * pht * rf(ny-j+1)))/
     +                 (2. + dt * pht * rf(ny-j+1))

      xx(i,j,k) = xxnT(i,ny-j+1,k) + xxnH(i,ny-j+1,k)


      yynT(i,ny-j+1,k) = (2.*dth*xl2m*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 yynT(i,ny-j+1,k)*
     +                 (2. - dt*rf(ny-j+1)))/(2. + dt*rf(ny-j+1))

      yynH(i,ny-j+1,k) = (2.*dth*xl*((u1(i+1,j,k) - u1(i,j,k)) +
     +                 (w1(i,j,k) - w1(i,j,k-1))) +
     +                 yynH(i,ny-j+1,k)*
     +                 (2. - dt * pht * rf(ny-j+1)))/
     +                 (2. + dt * pht * rf(ny-j+1))

      yy(i,j,k) = yynT(i,ny-j+1,k) + yynH(i,ny-j+1,k)


      zznT(i,ny-j+1,k) = (2.*dth*xl*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 zznT(i,ny-j+1,k)*
     +                 (2. - dt*rf(ny-j+1)))/(2. + dt*rf(ny-j+1))

      zznH(i,ny-j+1,k) = (2.*dth*(xl2m*(w1(i,j,k) - w1(i,j,k-1)) +
     +                 xl*(u1(i+1,j,k) - u1(i,j,k))) +
     +                 zznH(i,ny-j+1,k)*
     +                 (2. - dt * pht * rf(ny-j+1)))/
     +                 (2. + dt * pht * rf(ny-j+1))


      zz(i,j,k) = zznT(i,ny-j+1,k) + zznH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xzn(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ IN NORTH PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)
      

      xznT(i,ny-j+1,k) = xznT(i,ny-j+1,k)*
     +                 (2. - dt*rf(ny-j+1))/(2. + dt*rf(ny-j+1))

      xznH(i,ny-j+1,k) = (2.*xm*dth*((u1(i,j,k+1) - u1(i,j,k)) +
     +                 (w1(i,j,k) - w1(i-1,j,k))) +
     +                  xznH(i,ny-j+1,k) *
     +                  (2. - dt * pht * rf(ny-j+1)))/
     +                  (2. + dt * pht * rf(ny-j+1))

      xz(i,j,k) = xznT(i,ny-j+1,k) + xznH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yzn(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ IN NORTH PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)
      

      yznT(i,ny-j+1,k) = (2.*dth*xm*
     +                 (w1(i,j+1,k) - w1(i,j,k)) +
     +                 yznT(i,ny-j+1,k)*
     +                 (2. - dt*rh(ny-j+1)))/(2. + dt*rh(ny-j+1))

      yznH(i,ny-j+1,k) = (2. * dth * xm*(v1(i,j,k+1) - v1(i,j,k)) +
     +                  yznH(i,ny-j+1,k) * 
     +                  (2. - dt * pht * rh(ny-j+1)))/
     +                  (2. + dt * pht * rh(ny-j+1))

      yz(i,j,k) = yznT(i,ny-j+1,k) + yznH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyn(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY IN NORTH PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)
      

      xynT(i,ny-j+1,k) = (2.*dth*xm*
     +                 (u1(i,j+1,k) - u1(i,j,k)) +
     +                 xynT(i,ny-j+1,k)*
     +                 (2. - dt*rh(ny-j+1)))/(2. + dt*rh(ny-j+1))

      xynH(i,ny-j+1,k) = (2.*xm*dth*(v1(i,j,k) - v1(i-1,j,k)) +
     +                  xynH(i,ny-j+1,k) *
     +                  (2. - dt * pht * rh(ny-j+1)))/
     +                  (2. + dt * pht * rh(ny-j+1))

      xy(i,j,k) = xynT(i,ny-j+1,k) + xynH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyzs(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ IN SOUTH PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl


      xxsT(i,j,k) = (2.*dth*xl*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              xxsT(i,j,k)*
     +              (2. - dt*rh(j)))/(2. + dt*rh(j))
      
      xxsH(i,j,k) = (2.*dth*(xl2m*(u1(i+1,j,k) - u1(i,j,k)) +
     +              xl*(w1(i,j,k) - w1(i,j,k-1))) +
     +              xxsH(i,j,k)*
     +              (2. - dt*pht*rh(j)))/(2. + dt*pht*rh(j))

      xx(i,j,k) = xxsT(i,j,k) + xxsH(i,j,k)


      yysT(i,j,k) = (2.*dth*xl2m*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              yysT(i,j,k)*
     +              (2. - dt*rh(j)))/(2. + dt*rh(j))

      yysH(i,j,k) = (2.*dth*xl*((u1(i+1,j,k) - u1(i,j,k)) +
     +              (w1(i,j,k) - w1(i,j,k-1))) +
     +              yysH(i,j,k)*
     +              (2. - dt*pht*rh(j)))/(2. + dt*pht*rh(j))

      yy(i,j,k) = yysT(i,j,k) + yysH(i,j,k)


      zzsT(i,j,k) = (2.*dth*xl*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              zzsT(i,j,k)*
     +              (2. - dt*rh(j)))/(2. + dt*rh(j))

      zzsH(i,j,k) = (2.*dth*(xl2m*(w1(i,j,k) - w1(i,j,k-1)) +
     +              xl*(u1(i+1,j,k) - u1(i,j,k)))+
     +              zzsH(i,j,k)*
     +              (2. - dt*pht*rh(j)))/(2. + dt*pht*rh(j))

      zz(i,j,k) = zzsT(i,j,k) + zzsH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xzs(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ IN SOUTH PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      xzsT(i,j,k) = xzsT(i,j,k)*
     +              (2. - dt*rh(j))/(2. + dt*rh(j)) 

      xzsH(i,j,k) = (2.*xm*dth*((u1(i,j,k+1) - u1(i,j,k)) +
     +              (w1(i,j,k) - w1(i-1,j,k))) +
     +              xzsH(i,j,k)*
     +              (2. - dt*pht*rh(j)))/(2. + dt*pht*rh(j)) 

      xz(i,j,k) = xzsT(i,j,k) + xzsH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yzs(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ IN SOUTH PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      yzsT(i,j,k) = (2.*dth*xm*
     +              (w1(i,j+1,k) - w1(i,j,k)) + 
     +              yzsT(i,j,k)*
     +              (2. - dt*rf(j)))/(2. + dt*rf(j))

      yzsH(i,j,k) = (2.*xm*dth*(v1(i,j,k+1) - v1(i,j,k)) +
     +              yzsH(i,j,k)*
     +              (2. - dt*pht*rf(j)))/(2. + dt*pht*rf(j))

      yz(i,j,k) = yzsT(i,j,k) + yzsH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xys(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY IN SOUTH PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      xysT(i,j,k) = (2*dth*xm*
     +              (u1(i,j+1,k) - u1(i,j,k)) +
     +              xysT(i,j,k)*
     +              (2. - dt*rf(j)))/(2. + dt*rf(j))

      xysH(i,j,k) = (2.* xm*dth*(v1(i,j,k) - v1(i-1,j,k))  +
     +              xysH(i,j,k) *
     +              (2. - dt*pht*rf(j)))/(2. + dt*pht*rf(j))

      xy(i,j,k) = xysT(i,j,k) + xysH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyzw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ IN WEST PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl


      xxwT(i,j,k) = (2.*dth*xl2m*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              xxwT(i,j,k)*
     +              (2. - dt*rf(i)))/(2. + dt*rf(i))

      xxwH(i,j,k) = (2.*dth*xl*((v1(i,j,k) - v1(i,j-1,k)) +
     +              (w1(i,j,k) - w1(i,j,k-1))) +
     +             xxwH(i,j,k) *
     +             (2. - dt*pht*rf(i)))/(2. + dt*pht*rf(i))

      xx(i,j,k) = xxwT(i,j,k) + xxwH(i,j,k)


      yywT(i,j,k) = (2.*dth*xl*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              yywT(i,j,k)*
     +              (2. - dt*rf(i)))/(2. + dt*rf(i))

      yywH(i,j,k) = (2.*dth*(xl2m*(v1(i,j,k) - v1(i,j-1,k)) +
     +              xl*(w1(i,j,k) - w1(i,j,k-1))) +
     +             yywH(i,j,k) *
     +             (2. - dt*pht*rf(i)))/(2. + dt*pht*rf(i))

      yy(i,j,k) = yywT(i,j,k) + yywH(i,j,k)


      zzwT(i,j,k) = (2.*dth*xl*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              zzwT(i,j,k)*
     +              (2. - dt*rf(i)))/(2. + dt*rf(i))

      zzwH(i,j,k) = (2.*dth*(xl*(v1(i,j,k) - v1(i,j-1,k)) +
     +              xl2m*(w1(i,j,k) - w1(i,j,k-1))) +
     +              zzwH(i,j,k) *
     +              (2. - dt*pht*rf(i)))/(2. + dt*pht*rf(i))

      zz(i,j,k) = zzwT(i,j,k) + zzwH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xzw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ IN WEST PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      xzwT(i,j,k) = (2.*dth*xm*
     +              (w1(i,j,k) - w1(i-1,j,k)) +
     +              xzwT(i,j,k)*
     +              (2. - dt*rh(i)))/(2. + dt*rh(i))

      xzwH(i,j,k) = (2.*xm*dth*(u1(i,j,k+1) - u1(i,j,k)) +
     +              xzwH(i,j,k) *
     +              (2. - dt*pht*rh(i)))/(2. + dt*pht*rh(i))

      xz(i,j,k) = xzwT(i,j,k) + xzwH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yzw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ IN WEST PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      yzwT(i,j,k) = yzwT(i,j,k)*
     +              (2. - dt*rf(i))/(2. + dt*rf(i))

      yzwH(i,j,k) = (2.*xm*dth*((v1(i,j,k+1) - v1(i,j,k)) +
     +              (w1(i,j+1,k) - w1(i,j,k)))+
     +              yzwH(i,j,k) *
     +              (2. - dt*pht*rf(i)))/(2. + dt*pht*rf(i))

      yz(i,j,k) = yzwT(i,j,k) + yzwH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY IN WEST PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      xywT(i,j,k) = (2*dth*xm*
     +              (v1(i,j,k) - v1(i-1,j,k)) +
     +              xywT(i,j,k)*
     +              (2. - dt*rh(i)))/(2. + dt*rh(i))

      xywH(i,j,k) = (2.*xm*dth*(u1(i,j+1,k) - u1(i,j,k)) +
     +              xywH(i,j,k)*
     +              (2. - dt*pht*rh(i)))/(2. + dt*pht*rh(i))

      xy(i,j,k) = xywT(i,j,k) + xywH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyze(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ IN EAST PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl


      xxeT(nx-i+1,j,k) = (2.*dth*xl2m*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 xxeT(nx-i+1,j,k)*
     +                 (2. - dt*rh(nx-i+1)))/(2. + dt*rh(nx-i+1))

      xxeH(nx-i+1,j,k) = (2.*dth*xl*((v1(i,j,k) - v1(i,j-1,k)) +
     +                 (w1(i,j,k) - w1(i,j,k-1))) +
     +                 xxeH(nx-i+1,j,k)*
     +                 (2. - dt*pht*rh(nx-i+1)))/
     +                 (2. + dt*pht*rh(nx-i+1))

      xx(i,j,k) = xxeT(nx-i+1,j,k) + xxeH(nx-i+1,j,k)


      yyeT(nx-i+1,j,k) = (2.*dth*xl*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 yyeT(nx-i+1,j,k)*
     +                 (2. - dt*rh(nx-i+1)))/(2. + dt*rh(nx-i+1))

      yyeH(nx-i+1,j,k) = (2.*dth*(xl2m*(v1(i,j,k) - v1(i,j-1,k)) +
     +                 xl*(w1(i,j,k) - w1(i,j,k-1)))+
     +                 yyeH(nx-i+1,j,k)*
     +                 (2. - dt*pht*rh(nx-i+1)))/
     +                 (2. + dt*pht*rh(nx-i+1))

      yy(i,j,k) = yyeT(nx-i+1,j,k) + yyeH(nx-i+1,j,k)


      zzeT(nx-i+1,j,k) = (2.*dth*xl*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 zzeT(nx-i+1,j,k)*
     +                 (2. - dt*rh(nx-i+1)))/(2. + dt*rh(nx-i+1))

      zzeH(nx-i+1,j,k) = (2.*dth*(xl*(v1(i,j,k) - v1(i,j-1,k)) +
     +                 xl2m*(w1(i,j,k) - w1(i,j,k-1))) +
     +                 zzeH(nx-i+1,j,k)*
     +                 (2. - dt*pht*rh(nx-i+1)))/
     +                 (2. + dt*pht*rh(nx-i+1))

      zz(i,j,k) = zzeT(nx-i+1,j,k) + zzeH(nx-i+1,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xze(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ IN EAST PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      xzeT(nx-i+1,j,k) = (2.*dth*xm*
     +                 (w1(i,j,k) - w1(i-1,j,k)) +
     +                 xzeT(nx-i+1,j,k)*
     +                 (2. - dt*rf(nx-i+1)))/(2. + dt*rf(nx-i+1))

      xzeH(nx-i+1,j,k) = (2.*xm*dth*(u1(i,j,k+1) - u1(i,j,k)) +
     +                 xzeH(nx-i+1,j,k)*
     +                 (2. - dt*pht*rf(nx-i+1)))/
     +                 (2. + dt*pht*rf(nx-i+1))

      xz(i,j,k) = xzeT(nx-i+1,j,k) + xzeH(nx-i+1,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yze(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ IN EAST PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      yzeT(nx-i+1,j,k) = yzeT(nx-i+1,j,k)*
     +                 (2. - dt*rh(nx-i+1))/(2. + dt*rh(nx-i+1))

      yzeH(nx-i+1,j,k) =(2.*xm*dth*((v1(i,j,k+1) - v1(i,j,k)) +
     +                 (w1(i,j+1,k) - w1(i,j,k)))+
     +                 yzeH(nx-i+1,j,k)*
     +                 (2. - dt*pht*rh(nx-i+1)))/
     +                 (2. + dt*pht*rh(nx-i+1))

      yz(i,j,k) = yzeT(nx-i+1,j,k) + yzeH(nx-i+1,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xye(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY IN EAST PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      xyeT(nx-i+1,j,k) = (2*dth*xm*
     +                 (v1(i,j,k) - v1(i-1,j,k)) +
     +                 xyeT(nx-i+1,j,k)*
     +                 (2. - dt*rf(nx-i+1)))/(2. + dt*rf(nx-i+1))

      xyeH(nx-i+1,j,k) = (2.*xm*dth*(u1(i,j+1,k) - u1(i,j,k)) +
     +                 xyeH(nx-i+1,j,k)*
     +                 (2. - dt*pht*rf(nx-i+1)))/
     +                 (2. + dt*pht*rf(nx-i+1))

      xy(i,j,k) = xyeT(nx-i+1,j,k) + xyeH(nx-i+1,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyzb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ IN BOTTOM PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl


      xxbT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              xxbT(i,j,k)*
     +              (2. - dt*rh(k)))/(2. + dt*rh(k))

      xxbH(i,j,k) = (2.* dth*(xl2m*(u1(i+1,j,k) - u1(i,j,k)) +
     +              xl*(v1(i,j,k) - v1(i,j-1,k))) +
     +              xxbH(i,j,k)*
     +              (2. - dt*pht*rh(k)))/(2. + dt*pht*rh(k))

      xx(i,j,k) = xxbT(i,j,k) + xxbH(i,j,k)


      yybT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              yybT(i,j,k)*
     +              (2. - dt*rh(k)))/(2. + dt*rh(k))

      yybH(i,j,k) = (2.* dth*(xl*(u1(i+1,j,k) - u1(i,j,k)) +
     +              xl2m*(v1(i,j,k) - v1(i,j-1,k))) +
     +              yybH(i,j,k)*
     +              (2. - dt*pht*rh(k)))/(2. + dt*pht*rh(k))

      yy(i,j,k) = yybT(i,j,k) + yybH(i,j,k)


      zzbT(i,j,k) = (2.*dth*xl2m*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              zzbT(i,j,k)*
     +              (2. - dt*rh(k)))/(2. + dt*rh(k))

      zzbH(i,j,k) = (2.*dth*(xl*(u1(i+1,j,k) - u1(i,j,k)) +
     +              xl*(v1(i,j,k) - v1(i,j-1,k))) +
     +              zzbH(i,j,k)*
     +              (2. - dt*pht*rh(k)))/(2. + dt*pht*rh(k))

      zz(i,j,k) = zzbT(i,j,k) + zzbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xzb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ IN BOTTOM PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      xzbT(i,j,k) = (2.*dth*xm*
     +              (u1(i,j,k+1) - u1(i,j,k)) +
     +              xzbT(i,j,k)*
     +              (2. - dt*rf(k)))/(2. + dt*rf(k))

      xzbH(i,j,k) = (2.*xm*dth*(w1(i,j,k) - w1(i-1,j,k)) +
     +              xzbH(i,j,k)*
     +              (2. - dt*pht*rf(k)))/(2. + dt*pht*rf(k))

      xz(i,j,k) = xzbT(i,j,k) + xzbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yzb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ IN BOTTOM PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      yzbT(i,j,k) = (2.*dth*xm*
     +              (v1(i,j,k+1) - v1(i,j,k)) +
     +              yzbT(i,j,k)*
     +              (2. - dt*rf(k)))/(2. + dt*rf(k))

      yzbH(i,j,k) = (2.* xm*dth*(w1(i,j+1,k) - w1(i,j,k))+
     +              yzbH(i,j,k)*
     +              (2. - dt*pht*rf(k)))/(2. + dt*pht*rf(k))

      yz(i,j,k) = yzbT(i,j,k) + yzbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY IN BOTTOM PLANE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      xybT(i,j,k) = xybT(i,j,k)*
     +              (2. - dt*rh(k))/(2. + dt*rh(k))

      xybH(i,j,k) = (2.* xm*dth*((u1(i,j+1,k) - u1(i,j,k)) +
     +              (v1(i,j,k) - v1(i-1,j,k))) +
     +              xybH(i,j,k)*
     +              (2. - dt*pht*rh(k)))/(2. + dt*pht*rh(k))

      xy(i,j,k) = xybT(i,j,k) + xybH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end


C STRESS PLANES =-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C VERTICAL STRESS EDGES =--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


      subroutine m_xyzsw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ ON THE SOUTH-WEST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl
   
      dc=rh(j) + pht*rf(i)
      xxsT(i,j,k) = (2.*dth*xl*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              xxsT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rf(i)
      xxwT(i,j,k) = (2.*dth*xl2m*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              xxwT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)
 
      dc=pht*(rh(j) + rf(i))
      xxsH(i,j,k) = (2.*dth*xl*(w1(i,j,k) - w1(i,j,k-1)) +
     +              xxsH(i,j,k) *
     +              (2. - dt*dc))/(2. + dt*dc)

      xx(i,j,k) = xxsT(i,j,k) + xxwT(i,j,k) + xxsH(i,j,k) 


      dc=rh(j) + pht*rf(i)
      yysT(i,j,k) = (2.*dth*xl2m*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              yysT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rf(i)
      yywT(i,j,k) = (2.*dth*xl*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              yywT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rf(i))
      yysH(i,j,k) = (2.*dth*xl*(w1(i,j,k) - w1(i,j,k-1)) +
     +              yysH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yy(i,j,k) = yysT(i,j,k) + yywT(i,j,k) + yysH(i,j,k)


      dc=rh(j) + pht*rf(i)
      zzsT(i,j,k) = (2.*dth*xl*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              zzsT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rf(i)
      zzwT(i,j,k) = (2.*dth*xl*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              zzwT(i,j,k)* 
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rf(i))
      zzsH(i,j,k) = (2.* dth*xl2m*(w1(i,j,k) - w1(i,j,k-1))+
     +              zzsH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      zz(i,j,k) = zzsT(i,j,k) + zzwT(i,j,k) + zzsH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xzsw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ ON THE SOUTH-WEST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(j) + pht*rh(i)
      xzsT(i,j,k) = xzsT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      dc=pht*rh(j) + rh(i)
      xzwT(i,j,k) = (2.*dth*xm*
     +              (w1(i,j,k) - w1(i-1,j,k)) +
     +              xzwT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rh(i))
      xzsH(i,j,k) = (2.*dth*xm*(u1(i,j,k+1) - u1(i,j,k)) +
     +              xzsH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xz(i,j,k) = xzsT(i,j,k) + xzwT(i,j,k) + xzsH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yzsw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ ON THE SOUTH-WEST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(j) + pht*rf(i)
      yzsT(i,j,k) = (2.*dth*xm*
     +              (w1(i,j+1,k) - w1(i,j,k)) +
     +              yzsT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + rf(i)
      yzwT(i,j,k) = yzwT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      dc=pht*(rf(j) + rf(i))
      yzsH(i,j,k) = (2.*dth*xm*(v1(i,j,k+1) - v1(i,j,k)) +
     +              yzsH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yz(i,j,k) = yzsT(i,j,k) + yzwT(i,j,k) + yzsH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xysw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY ON THE SOUTH-WEST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(j) + pht*rh(i)
      xysT(i,j,k) = (2*dth*xm*
     +              (u1(i,j+1,k) - u1(i,j,k)) +
     +              xysT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + rh(i)
      xywT(i,j,k) = (2*dth*xm*
     +              (v1(i,j,k) - v1(i-1,j,k)) +
     +              xywT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(j) + rh(i))
      xysH(i,j,k) = xysH(i,j,k) *
     +              (2. - dt*dc)/(2. + dt*dc)

      xy(i,j,k) = xysT(i,j,k) + xywT(i,j,k) + xysH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyzse(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ ON THE SOUTH-EAST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl

      dc=rh(j) + pht*rh(nx-i+1)
      xxsT(i,j,k) = (2.*dth*xl*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              xxsT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rh(nx-i+1)
      xxeT(nx-i+1,j,k) = (2.*dth*xl2m*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 xxeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rh(nx-i+1))
      xxsH(i,j,k) = (2.*dth*xl*(w1(i,j,k) - w1(i,j,k-1))+
     +              xxsH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xx(i,j,k) = xxsT(i,j,k) + xxeT(nx-i+1,j,k) + xxsH(i,j,k)


      dc=rh(j) + pht*rh(nx-i+1)
      yysT(i,j,k) = (2.*dth*xl2m*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              yysT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rh(nx-i+1)
      yyeT(nx-i+1,j,k) = (2.*dth*xl* 
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 yyeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rh(nx-i+1))
      yysH(i,j,k) = (2.*dth*xl*(w1(i,j,k) - w1(i,j,k-1)) +
     +              yysH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yy(i,j,k) = yysT(i,j,k) + yyeT(nx-i+1,j,k) + yysH(i,j,k)


      dc=rh(j) + pht*rh(nx-i+1)
      zzsT(i,j,k) = (2.*dth*xl*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              zzsT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rh(nx-i+1)
      zzeT(nx-i+1,j,k) = (2.*dth*xl*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 zzeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rh(nx-i+1))
      zzsH(i,j,k) = (2.*dth*xl2m*(w1(i,j,k) - w1(i,j,k-1)) +
     +              zzsH(i,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      zz(i,j,k) = zzsT(i,j,k) + zzeT(nx-i+1,j,k) + zzsH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xzse(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ ON THE SOUTH-EAST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(j) + pht*rf(nx-i+1)
      xzsT(i,j,k) = xzsT(i,j,k)* 
     +              (2. - dt*dc)/(2. + dt*dc)

      dc=pht*rh(j) + rf(nx-i+1)
      xzeT(nx-i+1,j,k) = (2.*dth*xm*
     +                 (w1(i,j,k) - w1(i-1,j,k)) +
     +                 xzeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rf(nx-i+1))
      xzsH(i,j,k) = (2.*dth*xm*(u1(i,j,k+1) - u1(i,j,k)) +
     +              xzsH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xz(i,j,k) = xzsT(i,j,k) + xzeT(nx-i+1,j,k) + xzsH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yzse(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ ON THE SOUTH-EAST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(j) + pht*rh(nx-i+1)
      yzsT(i,j,k) = (2.*dth*xm*
     +              (w1(i,j+1,k) - w1(i,j,k)) +
     +              yzsT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + rh(nx-i+1)
      yzeT(nx-i+1,j,k) = yzeT(nx-i+1,j,k)*
     +                 (2. - dt*dc)/(2. + dt*dc)

      dc=pht*(rf(j) + rh(nx-i+1))
      yzsH(i,j,k) = (2.* dth*xm*(v1(i,j,k+1) - v1(i,j,k))+
     +              yzsH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yz(i,j,k) = yzsT(i,j,k) + yzeT(nx-i+1,j,k) + yzsH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyse(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY ON THE SOUTH-EAST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(j) + pht*rf(nx-i+1)
      xysT(i,j,k) = (2*dth*xm*
     +              (u1(i,j+1,k) - u1(i,j,k)) +
     +              xysT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + rf(nx-i+1)
      xyeT(nx-i+1,j,k) = (2*dth*xm*
     +                 (v1(i,j,k) - v1(i-1,j,k)) +
     +                 xyeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(j) + rf(nx-i+1))
      xysH(i,j,k) = 2.*xysH(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      xy(i,j,k) = xysT(i,j,k) + xyeT(nx-i+1,j,k) + xysH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyznw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ ON THE NORTH-WEST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl


      dc=rf(ny-j+1) + rf(i) * pht 
      xxnT(i,ny-j+1,k) = (2.*dth*xl*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 xxnT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht * rf(ny-j+1) + rf(i)
      xxwT(i,j,k) = (2.*dth*xl2m*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              xxwT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)
 
      dc=pht*(rf(ny-j+1)+rf(i))
      xxnH(i,ny-j+1,k) = (2.* dth*xl*(w1(i,j,k) - w1(i,j,k-1)) +
     +              xxnH(i,ny-j+1,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xx(i,j,k) = xxnT(i,ny-j+1,k) + xxwT(i,j,k) + xxnH(i,ny-j+1,k)


      dc=rf(ny-j+1) + pht*rf(i)
      yynT(i,ny-j+1,k) = (2.*dth*xl2m*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 yynT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rf(i) 
      yywT(i,j,k) = (2.*dth*xl*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              yywT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(ny-j+1)+rf(i))
      yynH(i,ny-j+1,k) =(2.* dth*xl*(w1(i,j,k) - w1(i,j,k-1)) +
     +              yynH(i,ny-j+1,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yy(i,j,k) = yynT(i,ny-j+1,k) + yywT(i,j,k) + yynH(i,ny-j+1,k)

      
      dc=rf(ny-j+1) + pht*rf(i)
      zznT(i,ny-j+1,k) = (2.*dth*xl*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 zznT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht * rf(ny-j+1) + rf(i)
      zzwT(i,j,k) = (2.*dth*xl*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              zzwT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht * (rf(ny-j+1) + rf(i))
      zznH(i,ny-j+1,k) = (2.*dth*xl2m*(w1(i,j,k) - w1(i,j,k-1)) +
     +                   zznH(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      zz(i,j,k) = zznT(i,ny-j+1,k) + zzwT(i,j,k) + zznH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xznw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ ON THE NORTH-WEST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(ny-j+1) + pht*rh(i)
      xznT(i,ny-j+1,k) = xznT(i,ny-j+1,k)*
     +                 (2. - dt*dc)/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rh(i)
      xzwT(i,j,k) = (2.*dth*xm*
     +              (w1(i,j,k) - w1(i-1,j,k)) +
     +              xzwT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(ny-j+1) + rh(i))
      xznH(i,ny-j+1,k) = (2.*dth*xm*(u1(i,j,k+1) - u1(i,j,k)) +
     +              xznH(i,ny-j+1,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xz(i,j,k) = xznT(i,ny-j+1,k) + xzwT(i,j,k) + xznH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yznw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ ON THE NORTH-WEST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(ny-j+1) + pht * rf(i)
      yznT(i,ny-j+1,k) = (2.*dth*xm*
     +                 (w1(i,j+1,k) - w1(i,j,k)) +
     +                 yznT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)     
  
      dc=pht*rh(ny-j+1) + rf(i)
      yzwT(i,j,k) = yzwT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      dc=pht*(rh(ny-j+1) + rf(i))
      yznH(i,ny-j+1,k) = (2.*dth*xm*(v1(i,j,k+1) - v1(i,j,k))+
     +                 yznH(i,ny-j+1,k)*
     +                (2. - dt*dc))/(2. + dt*dc)

      yz(i,j,k) = yznT(i,ny-j+1,k) + yzwT(i,j,k) + yznH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xynw(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY ON THE NORTH-WEST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(ny-j+1) + pht*rh(i)
      xynT(i,ny-j+1,k) = (2.*dth*xm*
     +                 (u1(i,j+1,k) - u1(i,j,k)) +
     +                 xynT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + rh(i)
      xywT(i,j,k) = (2.*dth*xm*
     +              (v1(i,j,k) - v1(i-1,j,k)) +
     +              xywT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(ny-j+1) + rh(i))
      xynH(i,ny-j+1,k) = xynH(i,ny-j+1,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      xy(i,j,k) = xynT(i,ny-j+1,k) + xywT(i,j,k) + xynH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyzne(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ ON THE NORTH-EAST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl

      dc=rf(ny-j+1) + pht * rh(nx-i+1)
      xxnT(i,ny-j+1,k) = (2.*dth*xl*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 xxnT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rh(nx-i+1)
      xxeT(nx-i+1,j,k) = (2.*dth*xl2m*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 xxeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(ny-j+1) + rh(nx-i+1))
      xxnH(i,ny-j+1,k) = (2.*dth*xl*(w1(i,j,k) - w1(i,j,k-1))+
     +                 xxnH(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      xx(i,j,k) = xxnT(i,ny-j+1,k) + xxeT(nx-i+1,j,k) + xxnH(i,ny-j+1,k)


      dc=rf(ny-j+1) + pht * rh(nx-i+1)
      yynT(i,ny-j+1,k) = (2.*dth*xl2m*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 yynT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rh(nx-i+1)
      yyeT(nx-i+1,j,k) = (2.*dth*xl*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 yyeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(ny-j+1) + rh(nx-i+1))
      yynH(i,ny-j+1,k) = (2.*dth*xl*(w1(i,j,k) - w1(i,j,k-1))+
     +                 yynH(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      yy(i,j,k) = yynT(i,ny-j+1,k) + yyeT(nx-i+1,j,k) +yynH(i,ny-j+1,k)


      dc=rf(ny-j+1) + pht * rh(nx-i+1)
      zznT(i,ny-j+1,k) = (2.*dth*xl*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 zznT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rh(nx-i+1)
      zzeT(nx-i+1,j,k) = (2.*dth*xl*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 zzeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(ny-j+1) + rh(nx-i+1))
      zznH(i,ny-j+1,k) = (2.*dth*xl2m*(w1(i,j,k) - w1(i,j,k-1))+
     +                 zznH(i,ny-j+1,k) *
     +                 (2. - dt*dc))/(2. + dt*dc)

      zz(i,j,k) = zznT(i,ny-j+1,k) + zzeT(nx-i+1,j,k) + zznH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xzne(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ ON THE NORTH-EAST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(ny-j+1) + pht * rf(nx-i+1)
      xznT(i,ny-j+1,k) = xznT(i,ny-j+1,k)*
     +                 (2. - dt*dc)/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rf(nx-i+1)
      xzeT(nx-i+1,j,k) = (2.*dth*xm*
     +                 (w1(i,j,k) - w1(i-1,j,k)) +
     +                 xzeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(ny-j+1) + rf(nx-i+1))
      xznH(i,ny-j+1,k) = (2.*dth*xm*(u1(i,j,k+1) - u1(i,j,k)) + 
     +                 xznH(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      xz(i,j,k) = xznT(i,ny-j+1,k) + xzeT(nx-i+1,j,k) + xznH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yzne(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ ON THE NORTH-EAST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)
    
      dc=rh(ny-j+1) + pht*rh(nx-i+1)
      yznT(i,ny-j+1,k) = (2.*dth*xm*
     +                 (w1(i,j+1,k) - w1(i,j,k)) +
     +                 yznT(i,ny-j+1,k)*
     +                 (2. - dt*rh(ny-j+1)))/(2. + dt*rh(ny-j+1))

      dc=pht*rh(ny-j+1) + rh(nx-i+1)
      yzeT(nx-i+1,j,k) = yzeT(nx-i+1,j,k)*
     +                 (2. - dt*dc)/(2. + dt*dc)

      dc=pht*(rh(ny-j+1) + rh(nx-i+1))
      yznH(i,ny-j+1,k) = (2.* dth*xm*(v1(i,j,k+1) - v1(i,j,k))+
     +                  yznH(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      yz(i,j,k) = yznT(i,ny-j+1,k) + yzeT(nx-i+1,j,k) + yznH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyne(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY ON THE NORTH-EAST EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(ny-j+1) + pht*rf(nx-i+1)
      xynT(i,ny-j+1,k) = (2*dth*xm*
     +                 (u1(i,j+1,k) - u1(i,j,k)) +
     +                 xynT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + rf(nx-i+1)
      xyeT(nx-i+1,j,k) = (2*dth*xm*
     +                 (v1(i,j,k) - v1(i-1,j,k)) +
     +                 xyeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      xynH(i,ny-j+1,k) = 2.*xynH(i,ny-j+1,k)*
     +                 (2. - dt*dc)/(2. + dt*dc)

      xy(i,j,k) = xynT(i,ny-j+1,k) + xyeT(nx-i+1,j,k) + xynH(i,ny-j+1,k)

   50 continue
      enddo
      enddo
      return
      end


C VERTICAL STRESS EDGES =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
C HORIZONTAL STRESS EDGES =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


      subroutine m_xyzsb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ ON THE SOUTH-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl

      dc=rh(j) + pht*rh(k)
      xxsT(i,j,k) = (2.*dth*xl*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              xxsT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rh(k)
      xxbT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              xxbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rh(k))
      xxbH(i,j,k) = (2.* dth*xl2m*(u1(i+1,j,k) - u1(i,j,k)) +
     +              xxbH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xx(i,j,k) = xxsT(i,j,k) + xxbT(i,j,k) + xxbH(i,j,k)


      dc=rh(j) + pht*rh(k)
      yysT(i,j,k) = (2.*dth*xl2m*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              yysT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rh(k)
      yybT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              yybT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rh(k))
      yybH(i,j,k) = (2.* dth*xl*(u1(i+1,j,k) - u1(i,j,k))+
     +              yybH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yy(i,j,k) = yysT(i,j,k) + yybT(i,j,k) + yybH(i,j,k)


      dc=rh(j) + pht*rh(k)
      zzsT(i,j,k) = (2.*dth*xl*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              zzsT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rh(k)
      zzbT(i,j,k) = (2.*dth*xl2m*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              zzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rh(k))
      zzbH(i,j,k) = (2.*dth*xl*(u1(i+1,j,k) - u1(i,j,k)) +
     +              zzbH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      zz(i,j,k) = zzsT(i,j,k) + zzbT(i,j,k) + zzbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xzsb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ ON THE SOUTH-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(j) + pht*rf(k)
      xzsT(i,j,k) = xzsT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      dc=pht*rh(j) + rf(k)
      xzbT(i,j,k) = (2.*dth*xm*
     +              (u1(i,j,k+1) - u1(i,j,k)) +
     +              xzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(j) + rf(k))
      xzbH(i,j,k) = (2.* dth*xm*(w1(i,j,k) - w1(i-1,j,k))+
     +              xzbH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xz(i,j,k) = xzsT(i,j,k) + xzbT(i,j,k) + xzbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yzsb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ ON THE SOUTH-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(j) + pht*rf(k)
      yzsT(i,j,k) = (2.*dth*xm*
     +              (w1(i,j+1,k) - w1(i,j,k)) +
     +              yzsT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + rf(k)
      yzbT(i,j,k) = (2.*dth*xm*
     +              (v1(i,j,k+1) - v1(i,j,k)) +
     +              yzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)
    
      dc=pht*(rf(j) + rf(k))
      yzbH(i,j,k) = yzbH(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      yz(i,j,k) = yzsT(i,j,k) + yzbT(i,j,k) + yzbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xysb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY ON THE SOUTH-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(j) + pht*rh(k)
      xysT(i,j,k) = (2*dth*xm*
     +              (u1(i,j+1,k) - u1(i,j,k)) +
     +              xysT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + rh(k)
      xybT(i,j,k) = xybT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      dc=pht*(rf(j) + rh(k))
      xybH(i,j,k) = (2.*dth*xm*(v1(i,j,k) - v1(i-1,j,k))+
     +              xybH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xy(i,j,k) = xysT(i,j,k) + xybT(i,j,k) + xybH(i,j,k)
      
   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyznb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ ON THE NORTH-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl

      dc=rf(ny-j+1) + pht*rh(k)
      xxnT(i,ny-j+1,k) = (2.*dth*xl*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 xxnT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rh(k)
      xxbT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              xxbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(ny-j+1) + rh(k))
      xxbH(i,j,k) = (2.*dth*xl2m*(u1(i+1,j,k) - u1(i,j,k)) +
     +              xxbH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xx(i,j,k) = xxnT(i,ny-j+1,k) + xxbT(i,j,k) + xxbH(i,j,k)


      dc=rf(ny-j+1) + pht*rh(k)
      yynT(i,ny-j+1,k) = (2.*dth*xl2m*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 yynT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rh(k)
      yybT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              yybT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(ny-j+1) + rh(k))
      yybH(i,j,k) = (2.*dth*xl*(u1(i+1,j,k) - u1(i,j,k)) +
     +              yybH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yy(i,j,k) = yynT(i,ny-j+1,k) + yybT(i,j,k) + yybH(i,j,k)


      dc=rf(ny-j+1) + pht*rh(k)
      zznT(i,ny-j+1,k) = (2.*dth*xl*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 zznT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rh(k)
      zzbT(i,j,k) = (2.*dth*xl2m*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              zzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(ny-j+1) + rh(k))
      zzbH(i,j,k) = (2.* dth*xl*(u1(i+1,j,k) - u1(i,j,k)) +
     +              zzbH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      zz(i,j,k) = zznT(i,ny-j+1,k) + zzbT(i,j,k) + zzbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xznb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ ON THE NORTH-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(ny-j+1) + pht*rf(k)
      xznT(i,ny-j+1,k) = xznT(i,ny-j+1,k)*
     +                 (2. - dt*dc)/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rf(k)
      xzbT(i,j,k) = (2.*dth*xm*
     +              (u1(i,j,k+1) - u1(i,j,k)) +
     +              xzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(ny-j+1) + rf(k))
      xzbH(i,j,k) = (2.*dth*xm*(w1(i,j,k) - w1(i-1,j,k)) +
     +              xzbH(i,j,k) *
     +              (2. - dt*dc))/(2. + dt*dc)

      xz(i,j,k) = xznT(i,ny-j+1,k) + xzbT(i,j,k) + xzbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yznb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ ON THE NORTH-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(ny-j+1) + pht*rf(k)
      yznT(i,ny-j+1,k) = (2.*dth*xm*
     +                 (w1(i,j+1,k) - w1(i,j,k)) +
     +                 yznT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + rf(k)
      yzbT(i,j,k) = (2.*dth*xm*
     +              (v1(i,j,k+1) - v1(i,j,k)) +
     +              yzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(ny-j+1) + rf(k))
      yzbH(i,j,k) = yzbH(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      yz(i,j,k) = yznT(i,ny-j+1,k) + yzbT(i,j,k) + yzbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xynb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY ON THE NORTH-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(ny-j+1) + pht*rh(k)
      xynT(i,ny-j+1,k) = (2*dth*xm*
     +                 (u1(i,j+1,k) - u1(i,j,k)) +
     +                 xynT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + rh(k)
      xybT(i,j,k) = xybT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      dc=pht*(rh(ny-j+1) + rh(k))
      xybH(i,j,k) = (2.* dth*xm*(v1(i,j,k) - v1(i-1,j,k)) +
     +              xybH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xy(i,j,k) = xynT(i,ny-j+1,k) + xybT(i,j,k) + xybH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyzwb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ ON THE WEST-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl

      dc=rf(i) + pht*rh(k)
      xxwT(i,j,k) = (2.*dth*xl2m*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              xxwT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(i) + rh(k)
      xxbT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              xxbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(i) + rh(k))
      xxbH(i,j,k) = (2.*dth*xl*(v1(i,j,k) - v1(i,j-1,k)) +
     +              xxbH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xx(i,j,k) = xxwT(i,j,k) + xxbT(i,j,k) + xxbH(i,j,k)


      dc=rf(i) + pht*rh(k)
      yywT(i,j,k) = (2.*dth*xl*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              yywT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(i) + rh(k)
      yybT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              yybT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(i) + rh(k))
      yybH(i,j,k) = (2.*dth*xl2m*(v1(i,j,k) - v1(i,j-1,k))+
     +              yybH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yy(i,j,k) = yywT(i,j,k) + yybT(i,j,k) + yybH(i,j,k)


      dc=rf(i) + pht*rh(k)
      zzwT(i,j,k) = (2.*dth*xl*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              zzwT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(i) + rh(k)
      zzbT(i,j,k) = (2.*dth*xl2m*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              zzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(i) + rh(k))
      zzbH(i,j,k) = (2.*dth*xl*(v1(i,j,k) - v1(i,j-1,k))+
     +              zzbH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      zz(i,j,k) = zzwT(i,j,k) + zzbT(i,j,k) + zzbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xzwb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ ON THE WEST-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(i) + pht*rf(k)
      xzwT(i,j,k) = (2.*dth*xm*
     +              (w1(i,j,k) - w1(i-1,j,k)) +
     +              xzwT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(i) + rf(k)
      xzbT(i,j,k) = (2.*dth*xm*
     +              (u1(i,j,k+1) - u1(i,j,k)) +
     +              xzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(i) + rf(k))
      xzbH(i,j,k) = xzbH(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      xz(i,j,k) = xzwT(i,j,k) + xzbT(i,j,k) + xzbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yzwb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ ON THE WEST-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(i) + pht*rf(k)
      yzwT(i,j,k) = yzwT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      dc=pht*rf(i) + rf(k)
      yzbT(i,j,k) = (2.*dth*xm*
     +              (v1(i,j,k+1) - v1(i,j,k)) +
     +              yzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(i) + rf(k))
      yzbH(i,j,k) = (2.*dth*xm*(w1(i,j+1,k) - w1(i,j,k))+
     +              yzbH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yz(i,j,k) = yzwT(i,j,k) + yzbT(i,j,k) + yzbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xywb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY ON THE WEST-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(i) + pht*rh(k)
      xywT(i,j,k) = (2*dth*xm*
     +              (v1(i,j,k) - v1(i-1,j,k)) +
     +              xywT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(i) + rh(k)
      xybT(i,j,k) = xybT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      dc=pht*(rh(i) + rh(k))
      xybH(i,j,k) = (2.*dth*xm*(u1(i,j+1,k) - u1(i,j,k))+
     +              xybH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xy(i,j,k) = xywT(i,j,k) + xybT(i,j,k) + xybH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyzeb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ ON THE EAST-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl

      dc=rh(nx-i+1) + pht*rh(k)
      xxeT(nx-i+1,j,k) = (2.*dth*xl2m*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 xxeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(nx-i+1) + rh(k)
      xxbT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              xxbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(nx-i+1) + rh(k))
      xxbH(i,j,k) = (2.*dth*xl*(v1(i,j,k) - v1(i,j-1,k))+
     +              xxbH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xx(i,j,k) = xxeT(nx-i+1,j,k) + xxbT(i,j,k) + xxbH(i,j,k)


      dc=rh(nx-i+1) + pht*rh(k)
      yyeT(nx-i+1,j,k) = (2.*dth*xl*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 yyeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(nx-i+1) + rh(k)
      yybT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              yybT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(nx-i+1) + rh(k))
      yybH(i,j,k) = (2.*dth*xl2m*(v1(i,j,k) - v1(i,j-1,k))+
     +              yybH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yy(i,j,k) = yyeT(nx-i+1,j,k) + yybT(i,j,k) + yybH(i,j,k)


      dc=rh(nx-i+1) + pht*rh(k)
      zzeT(nx-i+1,j,k) = (2.*dth*xl*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 zzeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(nx-i+1) + rh(k)
      zzbT(i,j,k) = (2.*dth*xl2m*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              zzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(nx-i+1) + rh(k))
      zzbH(i,j,k) = (2.*dth*xl*(v1(i,j,k) - v1(i,j-1,k))+
     +              zzbH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      zz(i,j,k) = zzeT(nx-i+1,j,k) + zzbT(i,j,k) + zzbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xzeb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ ON THE EAST-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(nx-i+1) + pht*rf(k)
      xzeT(nx-i+1,j,k) = (2.*dth*xm*
     +                 (w1(i,j,k) - w1(i-1,j,k)) +
     +                 xzeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(nx-i+1) + rf(k)
      xzbT(i,j,k) = (2.*dth*xm*
     +              (u1(i,j,k+1) - u1(i,j,k)) +
     +              xzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rf(nx-i+1) + rf(k))
      xzbH(i,j,k) = xzbH(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      xz(i,j,k) = xzeT(nx-i+1,j,k) + xzbT(i,j,k) + xzbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yzeb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ ON THE EAST-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(nx-i+1) + pht*rf(k)
      yzeT(nx-i+1,j,k) = yzeT(nx-i+1,j,k)*
     +                 (2. - dt*dc)/(2. + dt*dc)

      dc=pht*rh(nx-i+1) + rf(k)
      yzbT(i,j,k) = (2.*dth*xm*
     +              (v1(i,j,k+1) - v1(i,j,k)) +
     +              yzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*(rh(nx-i+1) + rf(k))
      yzbH(i,j,k) = (2.* dth*xm*(w1(i,j+1,k) - w1(i,j,k)) +
     +              yzbH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yz(i,j,k) = yzeT(nx-i+1,j,k) + yzbT(i,j,k) + yzbH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyeb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY ON THE EAST-BOTTOM EDGE

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(nx-i+1) + pht*rh(k)
      xyeT(nx-i+1,j,k) = (2*dth*xm*
     +                 (v1(i,j,k) - v1(i-1,j,k)) +
     +                 xyeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(nx-i+1) + rh(k)
      xybT(i,j,k) = xybT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      dc=pht*(rf(nx-i+1) + rh(k))
      xybH(i,j,k) = (2.*dth*xm*(u1(i,j+1,k) - u1(i,j,k)) +
     +              xybH(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xy(i,j,k) = xyeT(nx-i+1,j,k) + xybT(i,j,k) + xybH(i,j,k)

   50 continue
      enddo
      enddo
      return
      end


C HORIZONTAL STRESS EDGES =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
C STRESS CORNERS =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


      subroutine m_xyzswb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ ON THE SOUTH-WEST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl

      dc=rh(j) + pht*rf(i) + pht*rh(k)
      xxsT(i,j,k) = (2.*dth*xl*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              xxsT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rf(i) + pht*rh(k)
      xxwT(i,j,k) = (2.*dth*xl2m*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              xxwT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + pht*rf(i) + rh(k)
      xxbT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              xxbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xx(i,j,k) = xxsT(i,j,k) + xxwT(i,j,k) + xxbT(i,j,k)


      dc=rh(j) + pht*rf(i) + pht*rh(k)
      yysT(i,j,k) = (2.*dth*xl2m*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              yysT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rf(i) + pht*rh(k)
      yywT(i,j,k) = (2.*dth*xl*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              yywT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + pht*rf(i) + rh(k)
      yybT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              yybT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yy(i,j,k) = yysT(i,j,k) + yywT(i,j,k) + yybT(i,j,k)


      dc=rh(j) + pht*rf(i) + pht*rh(k)
      zzsT(i,j,k) = (2.*dth*xl*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              zzsT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rf(i) + pht*rh(k)
      zzwT(i,j,k) = (2.*dth*xl*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              zzwT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + pht*rf(i) + rh(k)
      zzbT(i,j,k) = (2.*dth*xl2m*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              zzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      zz(i,j,k) = zzsT(i,j,k) + zzwT(i,j,k) + zzbT(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xzswb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ ON THE SOUTH-WEST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(j) + pht*rh(i) + pht*rf(k)
      xzsT(i,j,k) = xzsT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      dc=pht*rh(j) + rh(i) + pht*rf(k)
      xzwT(i,j,k) = (2.*dth*xm*
     +              (w1(i,j,k) - w1(i-1,j,k)) +
     +              xzwT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + pht*rh(i) + rf(k)
      xzbT(i,j,k) = (2.*dth*xm*
     +              (u1(i,j,k+1) - u1(i,j,k)) +
     +              xzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xz(i,j,k) = xzsT(i,j,k) + xzwT(i,j,k) + xzbT(i,j,k) 

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yzswb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ ON THE SOUTH-WEST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(j) + pht*rf(i) + pht*rf(k)
      yzsT(i,j,k) = (2.*dth*xm*
     +              (w1(i,j+1,k) - w1(i,j,k)) +
     +              yzsT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + rf(i) + pht*rf(k)
      yzwT(i,j,k) = yzwT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      dc=pht*rf(j) + pht*rf(i) + rf(k)
      yzbT(i,j,k) = (2.*dth*xm*
     +              (v1(i,j,k+1) - v1(i,j,k)) +
     +              yzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yz(i,j,k) = yzsT(i,j,k) + yzwT(i,j,k) + yzbT(i,j,k) 

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyswb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY ON THE SOUTH-WEST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)
     
      dc=rf(j) + pht*rh(i) + pht*rh(k)
      xysT(i,j,k) = (2*dth*xm*
     +              (u1(i,j+1,k) - u1(i,j,k)) +
     +              xysT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)
 
      dc=pht*rf(j) + rh(i) + pht*rh(k)
      xywT(i,j,k) = (2*dth*xm*
     +              (v1(i,j,k) - v1(i-1,j,k)) +
     +              xywT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + pht*rh(i) + rh(k)
      xybT(i,j,k) = xybT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      xy(i,j,k) = xysT(i,j,k) + xywT(i,j,k) + xybT(i,j,k) 

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyzseb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ ON THE SOUTH-EAST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl

      dc=rh(j) + pht*rh(nx-i+1) + pht*rh(k)
      xxsT(i,j,k) = (2.*dth*xl*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              xxsT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rh(nx-i+1) + pht*rh(k)
      xxeT(nx-i+1,j,k) = (2.*dth*xl2m*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 xxeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + pht*rh(nx-i+1) + rh(k)
      xxbT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              xxbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)
 
      xx(i,j,k) = xxsT(i,j,k) + xxeT(nx-i+1,j,k) + xxbT(i,j,k) 


      dc=rh(j) + pht*rh(nx-i+1) + pht*rh(k)
      yysT(i,j,k) = (2.*dth*xl2m*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              yysT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rh(nx-i+1) + pht*rh(k)
      yyeT(nx-i+1,j,k) = (2.*dth*xl*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 yyeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + pht*rh(nx-i+1) + rh(k)
      yybT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              yybT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yy(i,j,k) = yysT(i,j,k) + yyeT(nx-i+1,j,k) + yybT(i,j,k) 


      dc=rh(j) + pht*rh(nx-i+1) + pht*rh(k)
      zzsT(i,j,k) = (2.*dth*xl*
     +              (v1(i,j,k) - v1(i,j-1,k)) +
     +              zzsT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + rh(nx-i+1) + pht*rh(k)
      zzeT(nx-i+1,j,k) = (2.*dth*xl*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 zzeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + pht*rh(nx-i+1) + rh(k)
      zzbT(i,j,k) = (2.*dth*xl2m*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              zzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)
 
      zz(i,j,k) = zzsT(i,j,k) + zzeT(nx-i+1,j,k) + zzbT(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xzseb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ ON THE SOUTH-EAST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(j) + pht*rf(nx-i+1) + pht*rf(k)
      xzsT(i,j,k) = xzsT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      dc=pht*rh(j) + rf(nx-i+1) + pht*rf(k)
      xzeT(nx-i+1,j,k) = (2.*dth*xm*
     +                 (w1(i,j,k) - w1(i-1,j,k)) +
     +                 xzeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(j) + pht*rf(nx-i+1) + rf(k)
      xzbT(i,j,k) = (2.*dth*xm*
     +              (u1(i,j,k+1) - u1(i,j,k)) +
     +              xzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xz(i,j,k) = xzsT(i,j,k) + xzeT(nx-i+1,j,k) + xzbT(i,j,k) 

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yzseb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ ON THE SOUTH-EAST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(j) + pht*rh(nx-i+1) + pht*rf(k)
      yzsT(i,j,k) = (2.*dth*xm*
     +              (w1(i,j+1,k) - w1(i,j,k)) +
     +              yzsT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + rh(nx-i+1) + pht*rf(k)
      yzeT(nx-i+1,j,k) = yzeT(nx-i+1,j,k)*
     +                 (2. - dt*dc)/(2. + dt*dc)

      dc=pht*rf(j) + pht*rh(nx-i+1) + rf(k)
      yzbT(i,j,k) = (2.*dth*xm*
     +              (v1(i,j,k+1) - v1(i,j,k)) +
     +              yzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yz(i,j,k) = yzsT(i,j,k) + yzeT(nx-i+1,j,k) + yzbT(i,j,k) 

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyseb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY ON THE SOUTH-EAST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(j) + pht*rf(nx-i+1) + pht*rh(k)
      xysT(i,j,k) = (2*dth*xm*
     +              (u1(i,j+1,k) - u1(i,j,k)) +
     +              xysT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + rf(nx-i+1) + pht*rh(k)
      xyeT(nx-i+1,j,k) = (2*dth*xm*
     +                 (v1(i,j,k) - v1(i-1,j,k)) +
     +                 xyeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(j) + pht*rf(nx-i+1) + rh(k)
      xybT(i,j,k) = xybT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      xy(i,j,k) = xysT(i,j,k) + xyeT(nx-i+1,j,k) + xybT(i,j,k) 

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyznwb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ ON THE NORTH-WEST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl

      dc=rf(ny-j+1) + pht*rf(i) + pht*rh(k)
      xxnT(i,ny-j+1,k) = (2.*dth*xl*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 xxnT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rf(i) + pht*rh(k)
      xxwT(i,j,k) = (2.*dth*xl2m*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              xxwT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + pht*rf(i) + rh(k)
      xxbT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              xxbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xx(i,j,k) = xxnT(i,ny-j+1,k) + xxwT(i,j,k) + xxbT(i,j,k) 

 
      dc=rf(ny-j+1) + pht*rf(i) + pht*rh(k)
      yynT(i,ny-j+1,k) = (2.*dth*xl2m*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 yynT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rf(i) + pht*rh(k)
      yywT(i,j,k) = (2.*dth*xl*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              yywT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + pht*rf(i) + rh(k)
      yybT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              yybT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yy(i,j,k) = yynT(i,ny-j+1,k) + yywT(i,j,k) + yybT(i,j,k) 

    
      dc=rf(ny-j+1) + pht*rf(i) + pht*rh(k)
      zznT(i,ny-j+1,k) = (2.*dth*xl*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 zznT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rf(i) + pht*rh(k)
      zzwT(i,j,k) = (2.*dth*xl*
     +              (u1(i+1,j,k) - u1(i,j,k)) +
     +              zzwT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + pht*rf(i) + rh(k)
      zzbT(i,j,k) = (2.*dth*xl2m*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              zzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      zz(i,j,k) = zznT(i,ny-j+1,k) + zzwT(i,j,k) + zzbT(i,j,k) 

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xznwb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ ON THE NORTH-WEST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(ny-j+1) + pht*rh(i) + pht*rf(k)
      xznT(i,ny-j+1,k) = xznT(i,ny-j+1,k)*
     +                 (2. - dt*dc)/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rh(i) + pht*rf(k)
      xzwT(i,j,k) = (2.*dth*xm*
     +              (w1(i,j,k) - w1(i-1,j,k)) +
     +              xzwT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + pht*rh(i) + rf(k)
      xzbT(i,j,k) = (2.*dth*xm*
     +              (u1(i,j,k+1) - u1(i,j,k)) +
     +              xzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xz(i,j,k) = xznT(i,ny-j+1,k) + xzwT(i,j,k) + xzbT(i,j,k) 

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yznwb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ ON THE NORTH-WEST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(ny-j+1) + pht*rf(i) + pht*rf(k)
      yznT(i,ny-j+1,k) = (2.*dth*xm*
     +                 (w1(i,j+1,k) - w1(i,j,k)) +
     +                 yznT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + rf(i) + pht*rf(k)
      yzwT(i,j,k) = yzwT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + pht*rf(i) + rf(k)
      yzbT(i,j,k) = (2.*dth*xm*
     +              (v1(i,j,k+1) - v1(i,j,k)) +
     +              yzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yz(i,j,k) = yznT(i,ny-j+1,k) + yzwT(i,j,k) + yzbT(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xynwb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY ON THE NORTH-WEST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(ny-j+1) + pht*rh(i) + pht*rh(k)
      xynT(i,ny-j+1,k) = (2*dth*xm*
     +                 (u1(i,j+1,k) - u1(i,j,k)) +
     +                 xynT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + rh(i) + pht*rh(k)
      xywT(i,j,k) = (2*dth*xm*
     +              (v1(i,j,k) - v1(i-1,j,k)) +
     +              xywT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + pht*rh(i) + rh(k)
      xybT(i,j,k) = xybT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      xy(i,j,k) = xynT(i,ny-j+1,k) + xywT(i,j,k) + xybT(i,j,k)

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyzneb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XYZ ON THE NORTH-EAST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xl = 1./lam(i,j,k)
      xl2m = 2./mu(i,j,k) + xl

      dc=rf(ny-j+1) + pht*rh(nx-i+1) + pht*rh(k)
      xxnT(i,ny-j+1,k) = (2.*dth*xl*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 xxnT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)
 
      dc=pht*rf(ny-j+1) + rh(nx-i+1) + pht*rh(k)
      xxeT(nx-i+1,j,k) = (2.*dth*xl2m*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 xxeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + pht*rh(nx-i+1) + rh(k)
      xxbT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              xxbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      xx(i,j,k) = xxnT(i,ny-j+1,k) + xxeT(nx-i+1,j,k) + xxbT(i,j,k)


      dc=rf(ny-j+1) + pht*rh(nx-i+1) + pht*rh(k)
      yynT(i,ny-j+1,k) = (2.*dth*xl2m*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 yynT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rh(nx-i+1) + pht*rh(k)
      yyeT(nx-i+1,j,k) = (2.*dth*xl*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 yyeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + pht*rh(nx-i+1) + rh(k)
      yybT(i,j,k) = (2.*dth*xl*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              yybT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yy(i,j,k) = yynT(i,ny-j+1,k) + yyeT(nx-i+1,j,k) + yybT(i,j,k)

 
      dc=rf(ny-j+1) + pht*rh(nx-i+1) + pht*rh(k)
      zznT(i,ny-j+1,k) = (2.*dth*xl*
     +                 (v1(i,j,k) - v1(i,j-1,k)) +
     +                 zznT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rh(nx-i+1) + pht*rh(k)
      zzeT(nx-i+1,j,k) = (2.*dth*xl*
     +                 (u1(i+1,j,k) - u1(i,j,k)) +
     +                 zzeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + pht*rh(nx-i+1) + rh(k)
      zzbT(i,j,k) = (2.*dth*xl2m*
     +              (w1(i,j,k) - w1(i,j,k-1)) +
     +              zzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      zz(i,j,k) = zznT(i,ny-j+1,k) + zzeT(nx-i+1,j,k) + zzbT(i,j,k) 

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xzneb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XZ ON THE NORTH-EAST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rf(ny-j+1) + pht*rf(nx-i+1) + pht*rf(k)
      xznT(i,ny-j+1,k) = xznT(i,ny-j+1,k)*
     +                 (2. - dt*dc)/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + rf(nx-i+1) + pht*rf(k)
      xzeT(nx-i+1,j,k) = (2.*dth*xm*
     +                 (w1(i,j,k) - w1(i-1,j,k)) +
     +                 xzeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rf(ny-j+1) + pht*rf(nx-i+1) + rf(k)
      xzbT(i,j,k) = (2.*dth*xm*
     +              (u1(i,j,k+1) - u1(i,j,k)) +
     +              xzbT(i,j,k)*
     +              (2. - dt*rf(k)))/(2. + dt*rf(k))

      xz(i,j,k) = xznT(i,ny-j+1,k) + xzeT(nx-i+1,j,k) + xzbT(i,j,k) 

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_yzneb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE YZ ON THE NORTH-EAST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(ny-j+1) + pht*rh(nx-i+1) + pht*rf(k)
      yznT(i,ny-j+1,k) = (2.*dth*xm*
     +                 (w1(i,j+1,k) - w1(i,j,k)) +
     +                 yznT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + rh(nx-i+1) + pht*rf(k)
      yzeT(nx-i+1,j,k) = yzeT(nx-i+1,j,k)*
     +                 (2. - dt*dc)/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + pht*rh(nx-i+1) + rf(k)
      yzbT(i,j,k) = (2.*dth*xm*
     +              (v1(i,j,k+1) - v1(i,j,k)) +
     +              yzbT(i,j,k)*
     +              (2. - dt*dc))/(2. + dt*dc)

      yz(i,j,k) = yznT(i,ny-j+1,k) + yzeT(nx-i+1,j,k) + yzbT(i,j,k) 

   50 continue
      enddo
      enddo
      return
      end



      subroutine m_xyneb(nxi,nxf,nyi,nyf,nzi,nzf,nx,ny,nz,dh,dt)

CCC   ROUTINE TO COMPUTE XY ON THE NORTH-EAST-BOTTOM CORNER

      use parstat

      dth = dt/dh

      do kk= nzi,nzf,kblock
      do jj= nyi,nyf,jblock
      do 50 k=kk,MIN(kk+kblock-1,nzf)
      do 50 j=jj,MIN(jj+jblock-1,nyf)
      do 50 i= nxi,nxf

      xm = 1./mu(i,j,k)

      dc=rh(ny-j+1) + pht*rf(nx-i+1) + pht*rh(k)
      xynT(i,ny-j+1,k) = (2*dth*xm*
     +                 (u1(i,j+1,k) - u1(i,j,k)) +
     +                 xynT(i,ny-j+1,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + rf(nx-i+1) + pht*rh(k)
      xyeT(nx-i+1,j,k) = (2*dth*xm*
     +                 (v1(i,j,k) - v1(i-1,j,k)) +
     +                 xyeT(nx-i+1,j,k)*
     +                 (2. - dt*dc))/(2. + dt*dc)

      dc=pht*rh(ny-j+1) + pht*rf(nx-i+1) + rh(k)
      xybT(i,j,k) = xybT(i,j,k)*
     +              (2. - dt*dc)/(2. + dt*dc)

      xy(i,j,k) = xynT(i,ny-j+1,k) + xyeT(nx-i+1,j,k) + xybT(i,j,k) 

   50 continue
      enddo
      enddo
      return
      end



C STRESS CORNERS =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=





