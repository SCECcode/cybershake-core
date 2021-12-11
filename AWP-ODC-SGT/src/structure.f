CCC   There are four options:
CCC   0. homogeneous mode
CCC   1. serial IO mode(small case) - the mesh file is specifined in IN3D
CCC   2. MPI IO mode(large case) - the mesh file is specifined in IN3D (!!! NOTE: npx*npy*npz >= nz !!!)
CCC   3. partitioned mode - partioned files mediaxxxxx.bin (where xxxxx is rank) 
CCC      are located at /input_rst/mediapart/
CCC   4. post-processed partitioned mode - partitioned files mediaIxxxxx.bin and
CCC      maxmin file mediaImaxmin.bin are located at /input_rst/mediapart/
CCC
CCC   nvar is for selecting number of variables in a mesh point. 3,5 and 8 is currently possible.
CCC

      subroutine inimesh(rank,ranksize,comm,master,invel,dims,coords,maxdim,
     +                   npx,npy,npz,nxt,nyt,nzt,nx,ny,nz,vpe,vse,dde,taumin,taumax,
     +                   fl,fh,fp,bsize,y_dt,idyna,nve,readoption,nvar,iost,partdeg,SoCalQ)

CCC   ROUTINE TO INITIALIZE SPACE WITH PREPPED 3D STRUCTURE
      use parstat
      implicit none

      integer, parameter :: filen=180
      character (len=filen) :: mediafile,invel,meshfile

      integer :: maxdim
      real :: taumax,taumin
      integer(i8) ::  p,r,e,f,g,l,m,n
      integer :: npx,npy,npz, nxt,nyt,nzt, nx,ny,nz, nxn,nyn,nzn, nxp,nyp,nzp
      integer(i8) :: nxt_l,nyt_l,nzt_l,ny_l,nx_l,nvar_l
      integer(i8), dimension (3):: dims_l
      integer(i8) :: i,j,k,a,b,c
      integer :: rank,ranksize,comm,master,err
      real :: tmax,tmin,w0,ww1,w2,qpinv,qsinv,pi
      integer, dimension(maxdim) :: dims,coords

      real, dimension(2) :: vpe,vse,dde,tmpvpe,tmpvse,tmpdde
      real :: fl,fh,fp
      integer :: bsize

      real :: tmp1, tmp2
      real :: vpvs
      integer :: y_dt,idyna,nve,readoption,SoCalQ

!     for SGT
      real :: tmpo1, tmpo2 

!     for option 1
      integer :: buf, tag

!     for option 2
      integer(KIND=MPI_ADDRESS_KIND),parameter :: EIGHT = 8
      integer(KIND=MPI_ADDRESS_KIND),parameter :: FIVE = 5
      integer(KIND=MPI_ADDRESS_KIND),parameter :: WSIZE = 4

      integer(KIND=MPI_ADDRESS_KIND) :: partdeg_l
      integer :: partdeg, est_partdeg, yindex
      real :: partdeg_f

      integer :: x,y,z,nvar,act_nvar,var_offset,sourcerank,mpirank,zrank,flipped_zrank
      real*8 :: t1,t2,t3,t4,t5,t6,total
      integer :: COMM_IO, COMM_RECV

      integer :: zpos
      
      integer :: varnum(1), bsize2
      real, dimension (:), allocatable :: tmpx
      real, dimension (:,:,:,:,:), allocatable :: tmpx2
      real, dimension (:,:,:,:,:), allocatable :: xyplane
      real, dimension (:,:,:,:), allocatable :: cube

      integer :: merr
      integer(MPI_OFFSET_KIND) :: disp
      integer :: PLANEUNIT, CUBEPLANEUNIT
      integer :: mfh
      integer,dimension(MPI_STATUS_SIZE) :: mystatus

      INTEGER, dimension(:), allocatable :: requests
      INTEGER, dimension(:,:), allocatable :: mpistatus
      
!     for option 3 
      integer :: iost,iwrite

      if (rank.eq.0) print *, "inimesh - ", fl,fh,fp,bsize     
      if (rank.eq.0) print *, "inimesh - ", nve,readoption,nvar,SoCalQ     

      call MPI_BCAST(invel,filen,MPI_CHARACTER,0,comm,merr)        

      ! CASE 0: Homogeneous Mesh Information
      if (readoption.eq.0) then
        if (rank.eq.0) print *, "initmesh - CASE 0, homogeneous mesh"
        pi=4.*atan(1.)
        taumax=1./(2*pi*fl)
        taumin=1./(2*pi*fh)
        call inihomo(nxt,nyt,nzt,vpe,vse,dde,y_dt,idyna)
      else
        allocate(tmpvp(1:nxt,1:nyt,1:nzt))
        allocate(tmpvs(1:nxt,1:nyt,1:nzt))
        allocate(tmpdd(1:nxt,1:nyt,1:nzt))
        allocate(tmppq(1:nxt,1:nyt,1:nzt))
        allocate(tmpsq(1:nxt,1:nyt,1:nzt))

        tmpvp=0
        tmpvs=0
        tmpdd=0
        tmppq=0
        tmpsq=0

        if (nvar.eq.8) then
          var_offset=3
          act_nvar=5
        else if (nvar.eq.5) then
          var_offset=0
          act_nvar=5
        else if (nvar.eq.3) then
          var_offset=0
          act_nvar=3
        end if

        if (rank.eq.0) then ! master
          print *,'nxt,nyt,nzt:',nxt,nyt,nzt
          print *,'nx,ny,nz:',nx,ny,nz
          print *,'nvar:',nvar
          print *,'IOST:',iost
          print *,'partdeg:',partdeg
          print *,'invel:',invel
        end if

        select case (readoption)
          ! CASE 1: Reading Single Mesh File with serial IO (small scale)
          ! Currently nve=1 is supported only
          case (1)
            if (rank.eq.0) print *, "CASE 1, serial IO - small files"

            if (rank.ne.0) then
              call MPI_RECV(buf,1,MPI_INTEGER,rank-1,MPI_ANY_TAG,comm,mystatus,merr)
            else
              close(27)
            end if

            allocate(tmpx(1:nxt*nvar))

c            if (nve==1) then
            open(100+rank,file=invel,access='direct',
     +             form='unformatted', recl=nvar*4*nxt)
c            end if

            x=rank/(npy*npz)
            y=(rank-x*npz*npy)/npz
            z=mod(rank-x*npy*npz,npz)
            z=(npz-z-1)

            do i=1,nzt
              do j=1,nyt
                m=z*npx*ny*nzt +
     +            y*npx*nyt +
     +            x +
     +            (j-1)*npx +
     +            (i-1)*npx*ny +
     +            1

                !print *, i, j, "-", m
                read(100+rank,rec=m) tmpx

                do k=1,nxt
                  tmpvp(k,j,i) = tmpx((k-1)*nvar+var_offset+1)
                  tmpvs(k,j,i) = tmpx((k-1)*nvar+var_offset+2)
                  tmpdd(k,j,i) = tmpx((k-1)*nvar+var_offset+3)
                  if (nvar.gt.3) then
                    tmppq(k,j,i) = tmpx((k-1)*nvar+var_offset+4)
                    tmpsq(k,j,i) = tmpx((k-1)*nvar+var_offset+5)
                  end if
                end do
              end do
            end do

            close(100+rank)

            deallocate(tmpx)

            if (rank.lt.(ranksize-1)) then
              call MPI_SEND(buf,1,MPI_INTEGER,rank+1,1234,comm,merr)
            end if

            print *, "inimesh - rank", rank

          ! CASE 2: Reading Single Mesh File with MPIIO (large scale)
          case (2)
            if (rank.eq.0) print *, "CASE 2, MPI - for big files"

            nx_l=nx
            ny_l=ny
            nvar_l=nvar

            allocate(requests(1:npx*npy+nzt))
            allocate(mpistatus(MPI_STATUS_SIZE,1:npx*npy+nzt))

            zrank=rank   ! 0..(nz-1)
 
            ! maximum 1GB for temporal mesh buffers
            partdeg_f=real(nx_l*ny_l*(nvar_l+act_nvar)*4)/real(1024*1024*1024)
            if (rank.eq.0) print *,partdeg_f,"GB"
            if (partdeg_f.gt.1) then
              est_partdeg=floor(log(partdeg_f)/log(2.))+1
              est_partdeg=2**est_partdeg
            else
              est_partdeg=1
            end if

	    ! verify the ny and npy are divisible by partdeg and npx*npy*npz >= nz*partdeg
            partdeg_l=partdeg
            if (rank.eq.0) print *,"Suggested amount of partitioning", est_partdeg
            if (rank.eq.0) print *,"Actual amount of partitioning", partdeg_l

            if (zrank<nz*partdeg) then
              call MPI_COMM_SPLIT(comm,1,0,COMM_IO,merr)        
              call MPI_FILE_OPEN(COMM_IO,invel,MPI_MODE_RDONLY,MPI_INFO_NULL,mfh,merr)

              ! new data type containing a contiguous data chunks
              call MPI_Type_contiguous(nx_l*(ny_l/partdeg_l)*nvar_l,MPI_REAL,PLANEUNIT,merr)
              call MPI_Type_commit(PLANEUNIT,merr)

              allocate(tmpx2(nvar,1:nxt,1:npx,1:nyt,1:npy/partdeg))        
              allocate(xyplane(act_nvar,1:nxt,1:nyt,1:npx,1:npy/partdeg))

              t1=mpi_wtime()
              disp=int(zrank,MPI_ADDRESS_KIND)*nx_l*(ny_l/partdeg_l)*nvar_l*WSIZE ! inverted z-axis, an absolute offset in bytes
              call MPI_File_set_view(mfh,disp,MPI_REAL,PLANEUNIT,
     +                               "native",MPI_INFO_NULL,merr)
              disp=0
              if (rank.eq.0) print *, "read_at_all", rank
              call MPI_File_read_at_all(mfh,disp,tmpx2,nx*(ny/partdeg)*nvar,MPI_REAL,mystatus,merr)

              if (rank.eq.0) print *, "read_at_all passed", rank
              t2=mpi_wtime()

              !call MPI_BARRIER(COMM_IO, merr)

              if (rank.eq.0) then ! master
                !write (*,'("Reading time =",f10.3," sec")') t2-t1
                print *,"Reading time =",t2-t1," sec"
              end if

              call MPI_FILE_CLOSE(mfh,merr)
              allocate(cube(act_nvar,1:nxt,1:nyt,1:nzt))

              yindex=mod(zrank,partdeg)
              do i=1,npy/partdeg
                do j=1,npx
                  do k=1,nyt
                    do l=1,nxt
                      xyplane(1,l,k,j,i)=tmpx2(var_offset+1,l,j,k,i)
                      xyplane(2,l,k,j,i)=tmpx2(var_offset+2,l,j,k,i)
                      xyplane(3,l,k,j,i)=tmpx2(var_offset+3,l,j,k,i)
                      if (nvar.gt.3) then
                        xyplane(4,l,k,j,i)=tmpx2(var_offset+4,l,j,k,i)
                        xyplane(5,l,k,j,i)=tmpx2(var_offset+5,l,j,k,i)
                      end if
                    end do
                  end do
                end do
              end do

              t3=mpi_wtime()
              if (rank.eq.0) then ! master          
                !write (*,'("Restructuring time =",f10.3," sec")') t3-t2
                print *,"Reconstructing time =",t3-t2," sec"
              end if
  
              !call MPI_BARRIER(COMM_IO, merr)

              ! DISTRIBUTING.....
              zpos=(zrank/partdeg)/nzt
              do i=yindex*(npy/partdeg)+1,(yindex+1)*npy/partdeg
                do j=1,npx
                  k=(i-yindex*(npy/partdeg)-1)*npx+j
                  
                  ! changing coordinate system: regular X-Y-Z -> MPI Z-Y-X 
                  z=npz-zpos-1   ! flip z-axis of the cubes
                  mpirank=(j-1)*npy*npz+(i-1)*npz + z
                  
                  !if (mod(zrank,nzt).eq.0) print *,"S",rank,mpirank,i,j,z

                  call MPI_ISEND(xyplane(1,1,1,j,i-yindex*(npy/partdeg)),nxt*nyt*act_nvar,MPI_REAL,
     $                           mpirank,mod(zrank/partdeg,nzt)+1,comm,requests(k),merr)
                end do
              end do

              ! RECEIVING..... ! do we need flipping z-axis?
              x=rank/(npy*npz)
              y=(rank-x*npy*npz)/npz
              z=mod(rank-x*npy*npz,npz)
              yindex=y/(npy/partdeg)

              z=npz-z-1   ! flip z-axis of the cubes

              !print *,"SR",rank,sourcerank,x,y,z

              do i=1,nzt
                call MPI_IRECV(cube(1,1,1,i),nxt*nyt*act_nvar,MPI_REAL,
     $                        (z*nzt+i-1)*partdeg+yindex,i,comm,requests(i+npx*npy/partdeg),merr)
              end do
              CALL MPI_WAITALL(npx*npy/partdeg+nzt,requests,mpistatus,merr)

              deallocate(tmpx2)
            else
              call MPI_COMM_SPLIT(comm,2,0,COMM_RECV,merr)        
              allocate(cube(act_nvar,1:nxt,1:nyt,1:nzt))

              ! RECEIVING..... ! do wee need flipping z-axis?
              x=rank/(npy*npz)
              y=(rank-x*npy*npz)/npz
              z=mod(rank-x*npy*npz,npz)
              yindex=y/(npy/partdeg)

              z=npz-z-1   ! flip z-axis of the cubes

              !print *,"R", rank,sourcerank,x,y,z

              do i=1,nzt
                call MPI_IRECV(cube(1,1,1,i),nxt*nyt*act_nvar,MPI_REAL,
     $                        (z*nzt+i-1)*partdeg+yindex,i,comm,requests(i),merr)
              end do
              CALL MPI_WAITALL(nzt,requests,mpistatus,merr)
            end if

            call MPI_BARRIER(comm, merr)
            t4=mpi_wtime()
            if (rank.eq.0) then ! master
              !write (*,'("Reception time =",f10.3," sec")') t4-t3
              print *,"Reception time =",t4-t3," sec"
            end if

            if (zrank<nz*partdeg) then
              deallocate(xyplane)
            end if

            ! For test purpose only
!            do iwrite=0,IOST
!              if (mod(rank,IOST+1).eq.iwrite) then
!                write(mediafile, '(a,i7.7,a)') 'input_rst/meshpart/media', rank,'.bin'
!                open(9, file=mediafile, form='unformatted', access='direct',
!     $               status='replace',recl=bsize*(nxt)*(nyt)*(nzt)*act_nvar)
!                write(9, rec=1) cube
!                close(9)
!              end if
!              call MPI_BARRIER(comm,merr)
!            end do
      
            t5=mpi_wtime()
            if (rank.eq.0) then ! master
              !write (*,'("Writng file time =",f10.3," sec")') t5-t4
              print *,"Writing file time =",t5-t4," sec"
            end if

            do i=1,nzt
              do j=1,nyt
                do k=1,nxt
                  tmpvp(k,j,i)=cube(1,k,j,i)
                  tmpvs(k,j,i)=cube(2,k,j,i)
                  tmpdd(k,j,i)=cube(3,k,j,i)
                  if (nvar.gt.3) then
                    tmppq(k,j,i)=cube(4,k,j,i)
                    tmpsq(k,j,i)=cube(5,k,j,i)
                  end if
                end do
              end do
            end do

            deallocate(cube)
            deallocate(requests)
            deallocate(mpistatus)

          ! CASE 3: Reading partitioned Mesh File (mediaxxxxx.bin)
          case (3)
            if (rank.eq.0) print *, "CASE 3, partitioned files"
            allocate(cube(act_nvar,-1:nxt+2,-1:nyt+2,-1:nzt+2))
        
            ! Maximally allow IOST amount of readers
            do iwrite=0,IOST
              if (rank.eq.0) print *,iwrite,"iterations out of ", IOST
              if (mod(rank,IOST+1).eq.iwrite) then
                write(mediafile,'(a,i3.3,i3.3,i3.3)') trim(invel),coords(1),coords(2),coords(3)
                print*, 'rank: ', rank, trim(mediafile), nxt, nyt, nzt
                open(9,file=mediafile,form='unformatted',access='direct',recl=bsize*nvar*(nxt+4)*(nyt+4)*(nzt+4))
                read(9, rec=1) cube
                close(9)
              end if
              call MPI_BARRIER(comm,merr)  
            end do

            do i=1,nzt
              do j=1,nyt
                do k=1,nxt
                  tmpvp(k,j,i)=cube(1,k,j,i)
                  tmpvs(k,j,i)=cube(2,k,j,i)
                  tmpdd(k,j,i)=cube(3,k,j,i)
                  if (nvar.gt.3) then
                    tmppq(k,j,i)=cube(4,k,j,i)
                    tmpsq(k,j,i)=cube(5,k,j,i)
                  end if
                end do
              end do
            end do

            deallocate(cube)

          ! CASE 4: Reading partitioned Preprocessed Mesh File (mediaIxxxxx.bin and mediaImaxmin.bin)
          case (4)
            if (rank.eq.0) print *, "initmesh - CASE 4, pre-processed partitioned files"
            write(mediafile, '(a,i7.7,a)') 'input_rst/mediapart/mediaI', rank,'.bin'
            if(rank.eq.0) print *, 'Read media file'
            open(9, file=mediafile, form='unformatted', access='direct',
     $          recl=bsize*(nxt+4)*(nyt+4)*(nzt+4)*5, status='old',
     $          iostat=err)
            read(9, rec=1) d1,mu,lam,qp,qs
            close(9)
            pi=4.*atan(1.)
            taumax=1./(2*pi*fl)
            taumin=1./(2*pi*fh)
C           11/24/2008 Kwangyoon added code to use global maxmin file named mediamaxmin.bin
C           Local maxmin files are still supported for compatibility purpose
            write(mediafile, '(a,i7.7,a)') 'input_rst/mediapart/mediaImaxmin', rank, '.bin'
            open(9, file=mediafile, form='unformatted', access='direct',
     $           status='old',recl=bsize*6,iostat=err)
            if (err.eq.0) then ! no error
              if(rank.eq.0) write(0, *) 'Reading Local maxmin media file'
              read(9, rec=1) vse(1),vse(2),vpe(1),vpe(2),dde(1),dde(2)
              close(9)
            else
              write(mediafile, '(a)') 'input_rst/mediapart/mediaImaxmin.bin'
              if (rank.eq.0) then
                print *, 'Reading global maxmin media file'
                open(9, file=mediafile, form='unformatted', 
     $               access='direct', status='old',recl=bsize*6,iostat=err)
                read(9, rec=1) vse(1),vse(2),vpe(1),vpe(2),dde(1),dde(2)
                close(9)
              end if
              call MPI_BCAST(vse,2,MPI_REAL,master,comm,err)
              call MPI_BCAST(vpe,2,MPI_REAL,master,comm,err)
              call MPI_BCAST(dde,2,MPI_REAL,master,comm,err)
            end if

          case default
            if (rank.eq.0) print *, "unsupported option is chosen!!!!"

        end select
      end if

!     for CASES 1,2,3 and generates dl,mu,lam,qs,qp and vpe,vse,dde
      if ((readoption.gt.0).and.(readoption.lt.4)) then
        
        ! Verify this!!!!
        if (nvar.eq.3) then
          tmpsq = 0.05 * tmpvs
          tmppq = 2.0 * tmpsq
c  SGT
c          tmpsq = 10000.0
c          tmppq = 10000.0
c  SGT
        end if

        pi=4.*atan(1.)
        w0=2*pi*fp
        ww1=2*pi*fl
        w2=2*pi*fh
        taumax=1./ww1
        taumin=1./w2
        tmp1=2./pi*(log(taumax/taumin))
        tmp2=2./pi*log(w0*taumin)

        vpe(1)=99999.
        vse(1)=99999.
        dde(1)=99999.
c calculate d1,mu,lam,qp,qs
        do i=1,nxt
          do j=1,nyt
            do k=1,nzt
              ! only if SoCalQ flag is on
              !if (tmpvs(i,j,k).lt.400.) then
              !   tmpvs(i,j,k)=400.
              !   tmpvp(i,j,k)=1200.
              !endif
              if (SoCalQ.eq.1) then
                tmpvs(i,j,k) = tmpvs(i,j,k)*(1+ ( log(w2/w0) )/(pi*tmpsq(i,j,k)) )
                tmpvp(i,j,k) = tmpvp(i,j,k)*(1+ ( log(w2/w0) )/(pi*tmppq(i,j,k)) )
                vpvs=tmpvp(i,j,k)/tmpvs(i,j,k)
                if (vpvs .lt. 1.45)  tmpvs(i,j,k)=tmpvp(i,j,k)/1.45
              end if ! SoCalQ
              mu(i,j,nzt+1-k) = 1./(tmpdd(i,j,k)*tmpvs(i,j,k)*tmpvs(i,j,k))
              lam(i,j,nzt+1-k) = 1./(tmpdd(i,j,k)*(tmpvp(i,j,k)*tmpvp(i,j,k)-2.*tmpvs(i,j,k)*tmpvs(i,j,k)))
              d1(i,j,nzt+1-k) = tmpdd(i,j,k)
c SGT
c               mu(i,j,k)  = 1./(tmpdd(i,j,k)*tmpvs(i,j,k)*tmpvs(i,j,k))
c               lam(i,j,k) = 1./(tmpdd(i,j,k)*(tmpvp(i,j,k)*tmpvp(i,j,k)-2.*tmpvs(i,j,k)*tmpvs(i,j,k)))
c               d1(i,j,k)  = tmpdd(i,j,k)
c               if (lam(i,j,k).le.0 .and. rank .eq. 0) then
c                     write(*,*) 'lam(i,j,k) is less than 0.'
c                     write(*,'(3E22.14)') tmpdd(i,j,k), tmpvp(i,j,k), tmpvs(i,j,k)
c                     tmpo1=tmpvp(i,j,k)
c                     tmpo2=tmpvs(i,j,k)
c                     write(*,'(2E22.14)') dble(tmpo1)*dble(tmpo1),2.d0*dble(tmpo2)*dble(tmpo2)
c                     write(*,'(1E22.14)') dble(tmpo1)*dble(tmpo1)-2.d0*dble(tmpo2)*dble(tmpo2)
c                     write(*,'(1E22.14)') tmpvp(i,j,k)*tmpvp(i,j,k)-2.*tmpvs(i,j,k)*tmpvs(i,j,k)
c               endif
               
c SGT
              if(nve==1) then
                if (tmppq(i,j,k) .le. 0.0) then
                  qpinv=0.0
                  qsinv=0.0
                else
                  qpinv=1./tmppq(i,j,k)
                  qsinv=1./tmpsq(i,j,k)
                end if

                tmppq(i,j,k)=tmp1*qpinv/(1.0-tmp2*qpinv)
                tmpsq(i,j,k)=tmp1*qsinv/(1.0-tmp2*qsinv)

                qp(i,j,nzt+1-k) = tmppq(i,j,k)
                qs(i,j,nzt+1-k) = tmpsq(i,j,k)
c SGT            
c                qp(i,j,k) = tmppq(i,j,k)
c                qs(i,j,k) = tmpsq(i,j,k)
c SGT
              endif
              if(tmpvs(i,j,k).lt.vse(1)) vse(1) = tmpvs(i,j,k)
              if(tmpvs(i,j,k).gt.vse(2)) vse(2) = tmpvs(i,j,k)
              if(tmpvp(i,j,k).lt.vpe(1)) vpe(1) = tmpvp(i,j,k)
              if(tmpvp(i,j,k).gt.vpe(2)) vpe(2) = tmpvp(i,j,k)
              if(tmpdd(i,j,k).lt.dde(1)) dde(1) = tmpdd(i,j,k)
              if(tmpdd(i,j,k).gt.dde(2)) dde(2) = tmpdd(i,j,k)
            end do
          end do
        end do

        call padmedia(nxt,nyt,nzt)
        call afsmedia(coords,dims,maxdim,nxt,nyt,nzt,nve)

        deallocate(tmpvp)
        deallocate(tmpvs)
        deallocate(tmpdd)
        deallocate(tmppq)
        deallocate(tmpsq)

        call MPI_BARRIER(comm,err)
        call MPI_ALLREDUCE(vse,tmpvse,2,MPI_REAL,MPI_MAX,comm,err)
        call MPI_ALLREDUCE(vpe,tmpvpe,2,MPI_REAL,MPI_MAX,comm,err)
        call MPI_ALLREDUCE(dde,tmpdde,2,MPI_REAL,MPI_MAX,comm,err)
        vse(2)=tmpvse(2)
        vpe(2)=tmpvpe(2)
        dde(2)=tmpdde(2)

        call MPI_ALLREDUCE(vse,tmpvse,2,MPI_REAL,MPI_MIN,comm,err)
        call MPI_ALLREDUCE(vpe,tmpvpe,2,MPI_REAL,MPI_MIN,comm,err)
        call MPI_ALLREDUCE(dde,tmpdde,2,MPI_REAL,MPI_MIN,comm,err)
        vse(1)=tmpvse(1)
        vpe(1)=tmpvpe(1)
        dde(1)=tmpdde(1)
      end if

      return
      end

      subroutine padmedia(nxt,nyt,nzt)
      use parstat
      implicit none

      integer :: nxt,nyt,nzt
      integer :: i,j,k

C     avoids wrong initialization for the outermost cells

      ! 5 Planes (except upper XY-plane, which is duplicated in afsmedia)
      do k=1,nzt
        do j=1,nyt
          lam(0,j,k) = lam(1,j,k)
          lam(nxt+1,j,k) = lam(nxt,j,k)
          mu(0,j,k) = mu(1,j,k)
          mu(nxt+1,j,k) = mu(nxt,j,k)
          d1(0,j,k) = d1(1,j,k)
          d1(nxt+1,j,k) = d1(nxt,j,k)
        end do
      end do

      do k=1,nzt
        do i=1,nxt
          lam(i,0,k) = lam(i,1,k)
          lam(i,nyt+1,k) = lam(i,nyt,k)
          mu(i,0,k) = mu(i,1,k)
          mu(i,nyt+1,k) = mu(i,nyt,k)
          d1(i,0,k) = d1(i,1,k)
          d1(i,nyt+1,k) = d1(i,nyt,k)
        end do
      end do

      do j=1,nyt
        do i=1,nxt
          lam(i,j,0)= lam(i,j,1)
          mu(i,j,0) = mu(i,j,1)
          d1(i,j,0) = d1(i,j,1)
        end do
      end do

      ! 12 border lines
      do i=1,nxt
        lam(i,0,0)= lam(i,1,1)
        mu(i,0,0) = mu(i,1,1)
        d1(i,0,0) = d1(i,1,1)
        lam(i,nyt+1,0)= lam(i,nyt,1)
        mu(i,nyt+1,0) = mu(i,nyt,1)
        d1(i,nyt+1,0) = d1(i,nyt,1)
        lam(i,0,nzt+1)= lam(i,1,nzt)
        mu(i,0,nzt+1) = mu(i,1,nzt)
        d1(i,0,nzt+1) = d1(i,1,nzt)
        lam(i,nyt+1,nzt+1)= lam(i,nyt,nzt)
        mu(i,nyt+1,nzt+1) = mu(i,nyt,nzt)
        d1(i,nyt+1,nzt+1) = d1(i,nyt,nzt)
      end do

      do i=1,nyt
        lam(0,i,0)= lam(1,i,1)
        mu(0,i,0) = mu(1,i,1)
        d1(0,i,0) = d1(1,i,1)
        lam(nxt+1,i,0)= lam(nxt,i,1)
        mu(nxt+1,i,0) = mu(nxt,i,1)
        d1(nxt+1,i,0) = d1(nxt,i,1)
        lam(0,i,nzt+1)= lam(1,i,nzt)
        mu(0,i,nzt+1) = mu(1,i,nzt)
        d1(0,i,nzt+1) = d1(1,i,nzt)
        lam(nxt+1,i,nzt+1)= lam(nxt,i,nzt)
        mu(nxt+1,i,nzt+1) = mu(nxt,i,nzt)
        d1(nxt+1,i,nzt+1) = d1(nxt,i,nzt)
      end do

      do i=1,nzt
        lam(0,0,i)= lam(1,1,i)
        mu(0,0,i) = mu(1,1,i)
        d1(0,0,i) = d1(1,1,i)
        lam(nxt+1,0,i)= lam(nxt,1,i)
        mu(nxt+1,0,i) = mu(nxt,1,i)
        d1(nxt+1,0,i) = d1(nxt,1,i)
        lam(0,nyt+1,i)= lam(1,nyt,i)
        mu(0,nyt+1,i) = mu(1,nyt,i)
        d1(0,nyt+1,i) = d1(1,nyt,i)
        lam(nxt+1,nyt+1,i)= lam(nxt,nyt,i)
        mu(nxt+1,nyt+1,i) = mu(nxt,nyt,i)
        d1(nxt+1,nyt+1,i) = d1(nxt,nyt,i)
      end do

      ! 8 Corners
      lam(0,0,0)= lam(1,1,1)
      mu(0,0,0)= mu(1,1,1)
      d1(0,0,0)= d1(1,1,1)
      lam(nxt+1,0,0)= lam(nxt,1,1)
      mu(nxt+1,0,0)= lam(nxt,1,1)
      d1(nxt+1,0,0)= lam(nxt,1,1)
      lam(0,nyt+1,0)= lam(1,nyt,1)
      mu(0,nyt+1,0)= mu(1,nyt,1)
      d1(0,nyt+1,0)= d1(1,nyt,1)
      lam(0,0,nzt+1)= lam(1,1,nzt)
      mu(0,0,nzt+1)= mu(1,1,nzt)
      d1(0,0,nzt+1)= d1(1,1,nzt)
      lam(nxt+1,0,nzt+1)= lam(nxt,1,nzt)
      mu(nxt+1,0,nzt+1)= mu(nxt,1,nzt)
      d1(nxt+1,0,nzt+1)= d1(nxt,1,nzt)
      lam(nxt+1,nyt+1,0)= lam(nxt,nyt,1)
      mu(nxt+1,nyt+1,0)= mu(nxt,nyt,1)
      d1(nxt+1,nyt+1,0)= d1(nxt,nyt,1)
      lam(0,nyt+1,nzt+1)= lam(1,nyt,nzt)
      mu(0,nyt+1,nzt+1)= mu(1,nyt,nzt)
      d1(0,nyt+1,nzt+1)= d1(1,nyt,nzt)
      lam(nxt+1,nyt+1,nzt+1)= lam(nxt,nyt,nzt)
      mu(nxt+1,nyt+1,nzt+1)= mu(nxt,nyt,nzt)
      d1(nxt+1,nyt+1,nzt+1)= d1(nxt,nyt,nzt)
 
      return
      end


      subroutine afsmedia(coords,dims,maxdim,nxt,nyt,nzt,nve)
CCC   EXTENDS FREE SURFACE MEDIA ABOVE TWO LAYERS
      use parstat
      implicit none
      
      integer :: nxt,nyt,nzt
      integer :: nve,maxdim
      integer :: i,j,k
      
      integer, dimension(maxdim) :: coords,dims
      
      if((coords(3)+1)==dims(3)) then
        do i=1,nxt
          do j=1,nyt
            do k=nzt+1,nzt+1
              d1(i,j,k) = d1(i,j,nzt)
              mu(i,j,k) = mu(i,j,nzt)
              lam(i,j,k) = lam(i,j,nzt)
              if (nve.eq.1) then
                qp(i,j,k) = qp(i,j,nzt)
                qs(i,j,k) = qs(i,j,nzt)
              end if
            enddo
          enddo
        enddo
      endif
      
      return
      end

      subroutine inihomo(nxt,nyt,nzt,vpe,vse,dde,y_dt,idyna)
CCC   ROUTINE TO INITIALIZE HOMOGENEOUS STRUCTURE

C     VPE(1)  LOWEST P-VELOCITY               (REAL) (RETURN)
C     VPE(2)  HIGHEST P-VELOCITY              (REAL) (RETURN)
C     VSE(1)  LOWEST S-VELOCITY               (REAL) (RETURN)
C     VSE(2)  HIGHEST S-VELOCITY              (REAL) (RETURN)
C     PV      P VELOCITY                      (REAL)
C     PS      S VELOCITY                      (REAL)
C     DD      DENSITY                         (REAL)

      use parstat

      integer :: nxt,nyt,nzt,y_dt,idyna
      real :: vp,vs,dd,pq,sq
      real, dimension(2) :: vpe,vse,dde
      if (idyna==1) then
      vp=6000.
      vs=3464.
      dd=2670.
      else
      vp=4800.
      vs=2800.
      dd=2500.
      endif
      pq=100.
      sq=50.

c     vpe(1) = 0.8*vp
      vpe(1) = vp
      vpe(2) = vp
c     vse(1) = 0.8*vs
      vse(1) = vs
      vse(2) = vs
      dde(1) = dd
      dde(2) = dd

      do 60 k=0,nzt+1
      do 60 j=0,nyt+1
        if (idyna==1) then
           vp=vpe(2)
           vs=vse(2)
           dd=dde(2)
          if(j.lt.y_dt) then
            vp=vpe(1)
            vs=vse(1)
            dd=dde(1)
          end if
        endif
      do 60 i=0,nxt+1
        lam(i,j,k)=1./(dd*(vp*vp - 2.*vs*vs))
        mu(i,j,k)=1./(dd*vs*vs)
        d1(i,j,k)=dd
   60 continue

      return
      end
