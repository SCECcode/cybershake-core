CCC   'cerjan.f' CERJAN ABSORBING BOUNDARY CONDITIONS


      subroutine inicrj(coords,maxdim,nxt,nyt,nzt,ndt,nx,ny,arbc)

CCC   ROUTINE TO INITIALIZE DAMPING SERIES

      use parstat

      integer :: nxt,nyt,nzt,ndt,nx,ny
      integer :: nxp,nyp,nzp

      real :: arbc

      integer, dimension(maxdim) :: coords 

C     SCALE DAMPING FOR ABC WIDTH

      alpha=sqrt(-log(arbc))/ndt

      flag_cerjanx=0
      flag_cerjany=0
      flag_cerjanz=0
      xs_crj = 1
      xe_crj = -1
      ys_crj = 1
      ye_crj = -1
      zs_crj = 1
      ze_crj = -1

C     FIND CERJAN NODES AND INITIALIZE

      nxp = coords(1)*nxt + 1
      if(nxp <= ndt) then

        flag_cerjanx=1
        xs_crj(1) = 1
        xe_crj(1) = ndt
        ys_crj(1) = 1
        ye_crj(1) = nyt
        zs_crj(1) = 1
        ze_crj(1) = nzt

        allocate(dcrjx(xs_crj(1):xe_crj(1),ys_crj(1):ye_crj(1),zs_crj(1):ze_crj(1)))       
        dcrjx=1

        do i=xs_crj(1),xe_crj(1)
        do j=ys_crj(1),ye_crj(1)
        do k=zs_crj(1),ze_crj(1)

          nxp = coords(1)*nxt + i
          dcrjx(i,j,k) = dcrjx(i,j,k)*(exp(-((alpha*(ndt-nxp+1))**2)))

        end do
        end do
        end do

      else if( (nxp+nxt-1) >= (nx-ndt+1)) then

        flag_cerjanx=2
        xs_crj(1) = nxt-ndt+1
        xe_crj(1) = nxt
        ys_crj(1) = 1
        ye_crj(1) = nyt
        zs_crj(1) = 1
        ze_crj(1) = nzt

        allocate(dcrjx(xs_crj(1):xe_crj(1),ys_crj(1):ye_crj(1),zs_crj(1):ze_crj(1)))       
        dcrjx=1

        do i=xs_crj(1),xe_crj(1)
        do j=ys_crj(1),ye_crj(1)
        do k=zs_crj(1),ze_crj(1)

          nxp = coords(1)*nxt + i
          dcrjx(i,j,k) = dcrjx(i,j,k)*(exp(-((alpha*(ndt-(nx-nxp)))**2)))

        end do
        end do
        end do

      end if

      nyp = coords(2)*nyt + 1
      if(nyp <= ndt) then

         flag_cerjany=1
         if(flag_cerjanx==1) then
           xs_crj(2) = xe_crj(1)+1
           xe_crj(2) = nxt
         else if(flag_cerjanx==2) then    
           xs_crj(2) = 1 
           xe_crj(2) = xs_crj(1)-1
         else
           xs_crj(2) = 1
           xe_crj(2) = nxt
         end if
         ys_crj(2) = 1
         ye_crj(2) = ndt
         zs_crj(2) = 1
         ze_crj(2) = nzt

         do i=xs_crj(1),xe_crj(1)
         do j=ys_crj(2),ye_crj(2)
         do k=zs_crj(2),ze_crj(2)

           nyp = coords(2)*nyt + j
           dcrjx(i,j,k) = dcrjx(i,j,k)*(exp(-((alpha*(ndt-nyp+1))**2)))

         end do
         end do
         end do

         allocate(dcrjy(xs_crj(2):xe_crj(2),ys_crj(2):ye_crj(2),zs_crj(2):ze_crj(2)))
         dcrjy=1

         do i=xs_crj(2),xe_crj(2)
         do j=ys_crj(2),ye_crj(2)
         do k=zs_crj(2),ze_crj(2)

           nyp = coords(2)*nyt + j
           dcrjy(i,j,k) = dcrjy(i,j,k)*(exp(-((alpha*(ndt-nyp+1))**2)))

         end do
         end do
         end do

      else if( (nyp+nyt-1) >= (ny-ndt+1)) then

         flag_cerjany=2
         if(flag_cerjanx==1) then
           xs_crj(2) = xe_crj(1)+1
           xe_crj(2) = nxt
         else if(flag_cerjanx==2) then
           xs_crj(2) = 1
           xe_crj(2) = xs_crj(1)-1
         else
           xs_crj(2) = 1
           xe_crj(2) = nxt
         end if
         ys_crj(2) = nyt-ndt+1
         ye_crj(2) = nyt
         zs_crj(2) = 1
         ze_crj(2) = nzt

         do i=xs_crj(1),xe_crj(1)
         do j=ys_crj(2),ye_crj(2)
         do k=zs_crj(2),ze_crj(2)

           nyp = coords(2)*nyt + j
           dcrjx(i,j,k) = dcrjx(i,j,k)*(exp(-((alpha*(ndt-(ny-nyp)))**2)))

         end do
         end do
         end do

         allocate(dcrjy(xs_crj(2):xe_crj(2),ys_crj(2):ye_crj(2),zs_crj(2):ze_crj(2)))
         dcrjy=1


         do i=xs_crj(2),xe_crj(2)
         do j=ys_crj(2),ye_crj(2)
         do k=zs_crj(2),ze_crj(2)

           nyp = coords(2)*nyt + j
           dcrjy(i,j,k) = dcrjy(i,j,k)*(exp(-((alpha*(ndt-(ny-nyp)))**2)))

         end do
         end do
         end do

      end if

      nzp = coords(3)*nzt + 1
      if(nzp <= ndt) then

         flag_cerjanz=1
         
         if(flag_cerjanx==1) then
           xs_crj(3) = xe_crj(1)+1
           xe_crj(3) = nxt
         else if(flag_cerjanx==2) then
           xs_crj(3) = 1
           xe_crj(3) = xs_crj(1)-1
         else
           xs_crj(3) = 1
           xe_crj(3) = nxt
         end if

         ys_crj(3) = 1
         ye_crj(3) = nyt
         zs_crj(3) = 1
         ze_crj(3) = ndt

         do i=xs_crj(1),xe_crj(1)
         do j=ys_crj(3),ye_crj(3)
         do k=zs_crj(3),ze_crj(3)

           nzp = coords(3)*nzt + k
           dcrjx(i,j,k) = dcrjx(i,j,k)*(exp(-((alpha*(ndt-nzp+1))**2)))

         end do
         end do
         end do

         if(flag_cerjany==1) then
           ys_crj(3) = ye_crj(2)+1
           ye_crj(3) = nyt
         else if(flag_cerjany==2) then
           ys_crj(3) = 1
           ye_crj(3) = ys_crj(2)-1
         else
           ys_crj(3) = 1
           ye_crj(3) = nyt
         end if


         do i=xs_crj(3),xe_crj(3)
         do j=ys_crj(2),ye_crj(2)
         do k=zs_crj(3),ze_crj(3)

           nzp = coords(3)*nzt + k
           dcrjy(i,j,k) = dcrjy(i,j,k)*(exp(-((alpha*(ndt-nzp+1))**2)))

         end do
         end do
         end do

         allocate(dcrjz(xs_crj(3):xe_crj(3),ys_crj(3):ye_crj(3),zs_crj(3):ze_crj(3))) 
         dcrjz=1


         do i=xs_crj(3),xe_crj(3)
         do j=ys_crj(3),ye_crj(3)
         do k=zs_crj(3),ze_crj(3)
    
           nzp = coords(3)*nzt + k
           dcrjz(i,j,k) = dcrjz(i,j,k)*(exp(-((alpha*(ndt-nzp+1))**2)))

         end do
         end do
         end do

      end if

      return
      end

      subroutine addcrjstr(nxt, nyt, nzt )
      use parstat
      integer :: nxt, nyt, nzt

      if(flag_cerjanx>0) then 
         call addcrjstr_loop(dcrjx, xs_crj(1), xe_crj(1), ys_crj(1), ye_crj(1), zs_crj(1), ze_crj(1))
      end if
      if(flag_cerjany>0) then 
         call addcrjstr_loop(dcrjy, xs_crj(2), xe_crj(2), ys_crj(2), ye_crj(2), zs_crj(2), ze_crj(2))
      end if
      if(flag_cerjanz==1) then
         call addcrjstr_loop(dcrjz, xs_crj(3), xe_crj(3), ys_crj(3), ye_crj(3), zs_crj(3), ze_crj(3))
      end if

      return
      end 


      subroutine addcrjstr_loop(dcrj, i_b, i_e, j_b, j_e, k_b, k_e)
      use parstat
      integer :: i, j, k, kk, jj
      integer :: i_b, i_e, j_b, j_e, k_b, k_e 
      real, dimension(i_b:i_e,j_b:j_e,k_b:k_e) :: dcrj                                                      

      do kk= k_b,k_e,kblock
      do jj= j_b,j_e,jblock
      do 20 k= kk,min(kk+kblock-1,k_e)
      do 20 j= jj,min(jj+jblock-1,j_e)
                                                          
      xx(i_b:i_e,j,k) = xx(i_b:i_e,j,k)*dcrj(i_b:i_e,j,k)
      yy(i_b:i_e,j,k) = yy(i_b:i_e,j,k)*dcrj(i_b:i_e,j,k)
      zz(i_b:i_e,j,k) = zz(i_b:i_e,j,k)*dcrj(i_b:i_e,j,k)
      xz(i_b:i_e,j,k) = xz(i_b:i_e,j,k)*dcrj(i_b:i_e,j,k)
      yz(i_b:i_e,j,k) = yz(i_b:i_e,j,k)*dcrj(i_b:i_e,j,k)
      xy(i_b:i_e,j,k) = xy(i_b:i_e,j,k)*dcrj(i_b:i_e,j,k)
                                                              
   20 continue
      end do
      end do
                                                         
      return
      end

      subroutine addcrjvel(nxt, nyt, nzt)
      use parstat
      integer :: nxt, nyt, nzt

      if(flag_cerjanx>0) then
         call addcrjvel_loop(dcrjx, xs_crj(1), xe_crj(1), ys_crj(1), ye_crj(1), zs_crj(1), ze_crj(1))
      end if
      if(flag_cerjany>0) then                      
         call addcrjvel_loop(dcrjy, xs_crj(2), xe_crj(2), ys_crj(2), ye_crj(2), zs_crj(2), ze_crj(2))
      end if
      if(flag_cerjanz==1) then
         call addcrjvel_loop(dcrjz, xs_crj(3), xe_crj(3), ys_crj(3), ye_crj(3), zs_crj(3), ze_crj(3))
      end if

      return
      end
    

      subroutine addcrjvel_loop(dcrj, i_b, i_e, j_b, j_e, k_b, k_e)
      use parstat
      integer :: i,j,k,kk,jj
      integer :: i_b, i_e, j_b, j_e, k_b, k_e
      real, dimension(i_b:i_e,j_b:j_e,k_b:k_e) :: dcrj                       
                                                             
      do kk= k_b,k_e,kblock
      do jj= j_b,j_e,jblock
      do 20 k= kk,min(kk+kblock-1,k_e)
      do 20 j= jj,min(jj+jblock-1,j_e)

      u1(i_b:i_e,j,k) = u1(i_b:i_e,j,k)*dcrj(i_b:i_e,j,k)
      v1(i_b:i_e,j,k) = v1(i_b:i_e,j,k)*dcrj(i_b:i_e,j,k)
      w1(i_b:i_e,j,k) = w1(i_b:i_e,j,k)*dcrj(i_b:i_e,j,k)
                                                             
   20 continue
      end do
      end do
                                                            
      return
      end
