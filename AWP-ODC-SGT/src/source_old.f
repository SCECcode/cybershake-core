CCC   'source.f' SOURCE INITIALIZATION, LOCATION, ADDITION

      subroutine inisource(rank,ranksize,comm,master,insrc,dims,coords,maxdim,
     +                     npx,npy,npz,nxt,nyt,nzt,nx,ny,nz,
     +                     nsrc,read_step,nst,cur_step,srcproc,npsrc,
     +                     ix_d,fx_d,y_d,iz_d,fz_d,y_dt,
     +                     ifault,idyna)


      use parstat
      implicit none

      integer, parameter :: filen=180
      character (len=filen) :: insrc,sourcefile

      integer :: maxdim
      integer :: npx,npy,npz, nxt,nyt,nzt, nx,ny,nz, nxn,nyn,nzn, nxp,nyp,nzp
      integer :: i,j,k
      integer :: rank,ranksize,comm,master,err
      integer, dimension(maxdim) :: dims,coords
      integer :: nsrc,read_step,nst,cur_step,srcproc,npsrc
      integer :: ix_d,fx_d,y_d,iz_d,fz_d
      integer :: y_dt   !=-1
      integer :: ifault,idyna

CCC   ALLOCATE SOURCE LOCATION BUFFER (tpsrc)
      allocate(tpsrc(nsrc,maxdim))

CCC   ALLOCATE TEMPORARY SOURCE AND RECEIVER MEMORY
      select case (ifault)
        case (0:1)
        ! CASE 0: text based source file (small scale)
        ! CASE 1: binary source file (small scale)
          allocate(taxx(nsrc, read_step),tayy(nsrc, read_step),tazz(nsrc, read_step))
          allocate(taxz(nsrc, read_step),tayz(nsrc, read_step),taxy(nsrc, read_step))

        case (2)
        ! CASE 2: partitioned source file

        case (3:4)
        ! CASE 3,4: dynamic mode
        ! IDYNA=1 : SGSN mode
          if (idyna == 0) then
            allocate(tstrinxx(nsrc))
            allocate(tstrisxx(nsrc))
            allocate(tdmaxx(nsrc))
          else if(idyna == 1) then
            allocate(tstriny(nsrc))
            allocate(tstrisx(nsrc))
            allocate(tstrisz(nsrc))
            allocate(tdmax(nsrc))
            allocate(tmu_s(nsrc))
            allocate(tmu_d(nsrc))
          else
            print *, 'error!!, no dynamic model selected, choose idyna=1 for SGSN, or 0 for SG '
            stop
          end if
 
        case default

      end select

C     SOURCE PROCS, LOCATION AND INITIALIZATION
      select case (ifault)
        case (0:2) ! Wave propagation mode
        ! CASE 0: text based source file (small scale)
        ! CASE 1: binary source file (small scale)
        ! CASE 2: partitioned source file
          if ((ifault.lt.2).or.(nst==read_step)) then ! for CASE 0,1 or no time spliting (CASE 0,1,2)
            do cur_step= 1, nst, read_step
              if(rank==master) call rdsrcpos(nsrc,nst,nx,ny,nz,ifault,filen,idyna,read_step,cur_step)
              call MPI_BCAST(tpsrc,nsrc*maxdim,MPI_INTEGER,master,comm,err)
              if ((ifault==0).or.(ifault==1)) then ! for CASE 1,2
                call MPI_BCAST(taxx,nsrc*read_step,MPI_REAL,master,comm,err)
                call MPI_BCAST(tayy,nsrc*read_step,MPI_REAL,master,comm,err)
                call MPI_BCAST(tazz,nsrc*read_step,MPI_REAL,master,comm,err)
                call MPI_BCAST(taxz,nsrc*read_step,MPI_REAL,master,comm,err)
                call MPI_BCAST(tayz,nsrc*read_step,MPI_REAL,master,comm,err)
                call MPI_BCAST(taxy,nsrc*read_step,MPI_REAL,master,comm,err)
              endif

              call getsrc(coords,maxdim,rank,nxt,nyt,nzt,srcproc,nsrc,
     +             nst,npsrc,ifault,ix_d,fx_d,y_d,iz_d,fz_d,nz,y_dt,idyna,read_step,cur_step)
            end do 

            call MPI_BARRIER(comm,err)
          end if

          if (ifault==2 .and. (nst.gt.read_step)) then ! For CASE 2 and (nst>read_step)
            if (rank==master) print *, "call read_src_ifault_2  with idx 1"
            call read_src_ifault_2(coords,maxdim,comm,rank,master,nxt,nyt,nzt,srcproc,
     +                      nsrc,nst,npsrc,nz,1,read_step)
          end if   ! match 'if (ifault != 2)'  

        case (3:4)
        ! CASE 3,4: dynamic mode
        ! IDYNA=1 : SGSN mode
          if (rank==master) call rdsrcpos(nsrc,nst,nx,ny,nz,ifault,filen,idyna,nst,1)
          call MPI_BCAST(tpsrc,nsrc*maxdim,MPI_INTEGER,master,comm,err)

          if(idyna == 0) then
            call MPI_BCAST(tstrinxx,nsrc,MPI_REAL,master,comm,err)
            call MPI_BCAST(tstrisxx,nsrc,MPI_REAL,master,comm,err)
            call MPI_BCAST(tdmaxx,nsrc,MPI_REAL,master,comm,err)
          else
            call MPI_BCAST(tstrisx,nsrc,MPI_REAL,master,comm,err)
            call MPI_BCAST(tdmax,nsrc,MPI_REAL,master,comm,err)
            call MPI_BCAST(tstriny,nsrc,MPI_REAL,master,comm,err)
            call MPI_BCAST(tstrisz,nsrc,MPI_REAL,master,comm,err)
            call MPI_BCAST(tmu_s,nsrc,MPI_REAL,master,comm,err)
            call MPI_BCAST(tmu_d,nsrc,MPI_REAL,master,comm,err)
          end if

          call getsrc(coords,maxdim,rank,nxt,nyt,nzt,srcproc,nsrc,
     +         nst,npsrc,ifault,ix_d,fx_d,y_d,iz_d,fz_d,nz,y_dt,idyna,nst,1)

          call MPI_BARRIER(comm,err)

          if (ifault==2 .and. (nst.gt.read_step)) then
            if (rank==master) print *, "call read_src_ifault_2  with idx 1"
            call read_src_ifault_2(coords,maxdim,comm,rank,master,nxt,nyt,nzt,srcproc,
     +                             nsrc,nst,npsrc,nz,1,read_step)
          end if   ! match 'if (ifault != 2)'   

        case default

      end select
   
      return
      end
  
      subroutine rdsrcpos(nsrc,nst,nx,ny,nz,ifault,filen,idyna,read_step,cur_step)
CCC   ACQUIRE FAULT NODES AND S(T)
      use parstat
      implicit none
     
C     NUMBER OF FAULT NODES AND RUPTURE FUNCTION STEPS
      integer :: i,j,ifault,nsrc,nst,nx,ny,nz,filen,idyna,read_step,cur_step
      real :: taxx1,tayy1,tazz1,taxz1,tayz1,taxy1
      character (len=filen) :: tpsrcfile
      integer :: bsize, bsize1
      integer :: varnum(6)
      integer :: varnum1(1)
        
      real, dimension (:), allocatable :: tmpxx, tmpyy,tmpzz,tmpxz,tmpyz,tmpxy
        
      bsize=4*size(varnum)
      bsize1=4*size(varnum1)
      inquire(iolength=bsize)varnum
      inquire(iolength=bsize1)varnum1

C     in the case of ifault=2, basically skip taxx etc, not needed
      rewind(25)
      
      print *, 'ifault=', ifault
      if (ifault==3 .or. ifault==4) then ! DYNAMIC MODE AND TEXT SOURCE FILE
        if(idyna == 0)  then
          do i=1,nsrc
            read(33,*) tpsrc(i,1),tpsrc(i,2),tpsrc(i,3),tstrinxx(i),tstrisxx(i),tdmaxx(i)
            tpsrc(i,3) = nz + 1 - tpsrc(i,3)
          end do
        else
          do i=1,nsrc
            read(33,*) tpsrc(i,1),tpsrc(i,2),tpsrc(i,3),
     +                 tstriny(i),tstrisx(i),tdmax(i),tstrisz(i),tmu_s(i),tmu_d(i)
            tpsrc(i,3) = nz + 1 - tpsrc(i,3)
          end do
        end if
      else    ! WAV mode
        if (ifault==0) then ! TEXT SOURCE FILE and need enhancement!!!!
           do i=1,nsrc
              read(25,*) tpsrc(i,1),tpsrc(i,2),tpsrc(i,3)
c	Check bounds on source coordinates to make sure they're in a valid range
c	This is actually just a check to make sure we've included the coordinates
c	of the source on the first line and we didn't forget them.
              if((tpsrc(i,1)<1).OR.(tpsrc(i,1)>nx)) then
                  print *, "X source coordinate ", tpsrc(i,1), " is outside range 1-", nx
                  stop
              end if
              if((tpsrc(i,2)<1).OR.(tpsrc(i,2)>ny)) then
                  print *, "Y source coordinate ", tpsrc(i,2), " is outside range 1-", ny
                  stop
              end if
              if((tpsrc(i,3)<1).OR.(tpsrc(i,3)>nz)) then
                  print *, "Z source coordinate ", tpsrc(i,3), " is outside range 1-", nz
                  stop
              end if
              tpsrc(i,3) = nz + 1 - tpsrc(i,3)
              do j=1,cur_step-1
                read(25,*) taxx1,tayy1,tazz1,taxz1,tayz1,taxy1
              end do
              do j=1,read_step
                 read(25,*) taxx(i,j),tayy(i,j),tazz(i,j),taxz(i,j),tayz(i,j),taxy(i,j)
              end do
              do j=cur_step+read_step,nst
                read(25,*) taxx1,tayy1,tazz1,taxz1,tayz1,taxy1
              end do
           end do
        else if (ifault==1) then ! SMALL BINARY FILE and WRITE a tpsrc_tmp file
          allocate(tmpxx(nst))
          allocate(tmpyy(nst))
          allocate(tmpzz(nst))
          allocate(tmpxz(nst))
          allocate(tmpyz(nst))
          allocate(tmpxy(nst))

          write(tpsrcfile, '(a)') 'input_rst/srcpart/tpsrc/tpsrc_tmp'
          print *,'create subdomain tpsrc files '
          open(17, file=tpsrcfile,form='unformatted', access='direct',
     $      recl=bsize1*3)

          do i=1,nsrc
            read(25,rec=i) tpsrc(i,1),tpsrc(i,2),tpsrc(i,3),(tmpxx(j),tmpyy(j),
     +                     tmpzz(j),tmpxz(j),tmpyz(j),tmpxy(j),j=1,nst)
            write(17,rec=i) tpsrc(i,1),tpsrc(i,2),tpsrc(i,3)
            tpsrc(i,3) = nz + 1 - tpsrc(i,3)
            do j=1,read_step
              taxx(i,j)=tmpxx(cur_step+j-1)
              tayy(i,j)=tmpyy(cur_step+j-1)
              tazz(i,j)=tmpzz(cur_step+j-1)
              taxz(i,j)=tmpxz(cur_step+j-1)
              tayz(i,j)=tmpyz(cur_step+j-1)
              taxy(i,j)=tmpxy(cur_step+j-1)
            end do
          end do

          close(17)
          deallocate(tmpxx)
          deallocate(tmpyy)
          deallocate(tmpzz)
          deallocate(tmpxz)
          deallocate(tmpyz)
          deallocate(tmpxy)
        else if (ifault==2) then ! PARTITIONED BINARY FILE and read from tpsrc_tmp
          write(tpsrcfile, '(a)') 'input_rst/srcpart/tpsrc/tpsrc_tmp'
          open(17, file=tpsrcfile,form='unformatted', access='direct',
     $       recl=bsize1*3)
          do i=1,nsrc
            read(17, rec=i) tpsrc(i,1), tpsrc(i,2), tpsrc(i,3)
            tpsrc(i,3) = nz + 1 - tpsrc(i,3)
          end do
        endif
      end if

      return
      end


      subroutine getsrc(coords,maxdim,rank,nxt,nyt,nzt,srcproc,
     +       nsrc,nst,npsrc,ifault,ix_d,fx_d,y_d,iz_d,fz_d,nz,y_dt,idyna,read_step,cur_step)
            
CCC   FINDS SOURCE PROC(S) AND COORDINATES FOR FAULT NODE(S)
      use parstat
      implicit none
        
      integer, parameter :: filen=180
      integer :: i,j,k,nxt,nyt,nzt,maxdim,rank,srcproc,nz,cur_step,read_step
      integer :: nsrc,nst,npsrc,ifault,idyna
      integer, dimension(maxdim) :: coords
      integer :: varnum(1),bsize
      character (len=filen) :: faultfile
      
      integer :: ix_d,fx_d,y_d,iz_d,fz_d
      integer :: ix_dt,fx_dt,y_dt,iz_dt,fz_dt
      integer, dimension(3) :: floc    

      integer :: nxp,nyp,nzp
      integer :: itemp,igx,jgx,igz,jgz,ipx,jpx,ipy,ipz,jpz,jj

      bsize=4*size(varnum)
      inquire(iolength=bsize)varnum
      j=1   

C     SEARCH OVER ALL SOURCE(FAULT) NODES
      if (ifault==3 .or. ifault==4) then
        if(idyna == 1) then
          floc=maxval(tpsrc, dim=1)
          fx_dt=floc(1)+2
          y_dt=floc(2)
          fz_dt=floc(3)
          if(fz_dt.le.nz-2) fz_dt=fz_dt+2
          if(fz_dt.eq.nz-1) fz_dt=fz_dt+1
          floc=minval(tpsrc, dim=1)
          ix_dt=floc(1)-2
          iz_dt=floc(3)-2
          if(floc(2).ne.y_dt) then
            print *, 'ERROR, fault is no planar!!, y location of fault should be constant everywhere'
            stop
          endif

C         Y borders
          if(mod(y_dt,nyt)/=0) then
            itemp = mod(y_dt,nyt)
          else
            itemp = nyt
          end if
          if(itemp.le.2.or.itemp.ge.nyt-1) then
             print *, 'ERROR!!, fault is less than two grids from procesor border'
             stop
          endif

C         X borders
          if(mod(ix_dt,nxt)/=0) then
            itemp = mod(ix_dt,nxt)
          else
            itemp = nxt
          end if
          if(itemp.le.3) ix_dt=ix_dt-itemp-3
          if(itemp.ge.nxt-2) ix_dt=ix_dt+(nxt-itemp)-3

          if(mod(fx_dt,nxt)/=0) then
            itemp = mod(fx_dt,nxt)
          else
            itemp = nxt
          end if
          if(itemp.ge.nxt-2) fx_dt=fx_dt+(nxt-itemp)+4
          if(itemp.le.3) fx_dt=fx_dt-itemp+4

C         Z borders
          if(mod(iz_dt,nzt)/=0) then
            itemp = mod(iz_dt,nzt)
          else
            itemp = nzt
          end if
          if(itemp.le.3) iz_dt=iz_dt-itemp-3
          if(itemp.ge.nzt-2) iz_dt=iz_dt+(nzt-itemp)-3

          if(mod(fz_dt,nzt)/=0) then
            itemp = mod(fz_dt,nzt)
          else
            itemp = nzt
          end if
          if(itemp.ge.nzt-2.and.fz_dt.lt.nz-2) fz_dt=fz_dt+(nzt-itemp)+4
          if(itemp.le.3) fz_dt=fz_dt-itemp+4

C         determine local borders of fault on each processor
          fx_d=0
          y_d=0
          fz_d=0
          ix_d=0
          iz_d=0
          igx=(coords(1))*nxt + 1
          jgx=(coords(1)+1)*nxt

          igz=(coords(3))*nzt + 1
          jgz=(coords(3)+1)*nzt

          ipx=int((ix_dt-1)/nxt)
          jpx=int((fx_dt-1)/nxt)

          ipy=int((y_dt-1)/nyt)

          ipz=int((iz_dt-1)/nzt)
          jpz=int((fz_dt-1)/nzt)

          jj=ipy
          do k=ipz,jpz
            do i=ipx,jpx
              if(coords(1)==i .and. coords(2)==jj .and. coords(3)==k) then
                srcproc = rank
                fx_d=nxt
                y_d=mod(y_dt,nyt)
                fz_d=nzt
                ix_d=1
                iz_d=1

                if(igx.lt.ix_dt) ix_d=mod(ix_dt,nxt)
                if(jgx.gt.fx_dt) fx_d=mod(fx_dt,nxt)

                if(igz.lt.iz_dt) iz_d=mod(iz_dt,nzt)
                if(jgz.gt.fz_dt) fz_d=mod(fz_dt,nzt)
              end if
            end do
          end do
        end if    ! match if (idyna==1)
      end if     ! match if (ifault==3 .or. ifault==4)

      j=1
      do i=1,nsrc
        nxp = int((tpsrc(i,1)-1)/nxt)
        nyp = int((tpsrc(i,2)-1)/nyt)
        nzp = int((tpsrc(i,3)-1)/nzt)

C       DETERMINE IF NODE BELONGS TO PROC
        if(coords(1)==nxp .and. coords(2)==nyp .and. coords(3)==nzp) then
          if(mod(tpsrc(i,1),nxt)/=0) then
            tpsrc(j,1) = mod(tpsrc(i,1),nxt)
          else
            tpsrc(j,1) = nxt
          end if

          if(mod(tpsrc(i,2),nyt)/=0) then
            tpsrc(j,2) = mod(tpsrc(i,2),nyt)
          else
            tpsrc(j,2) = nyt
          end if

          if(mod(tpsrc(i,3),nzt)/=0) then
            tpsrc(j,3) = mod(tpsrc(i,3),nzt)
          else
            tpsrc(j,3) = nzt
          end if

          if (ifault==3 .or. ifault==4) then
            if(idyna == 0) then
              tstrisxx(j)=tstrisxx(i)
              tstrinxx(j)=tstrinxx(i)
              tdmaxx(j)=tdmaxx(i)
            else
              tstrisx(j)=tstrisx(i)
              tstrisz(j)=tstrisz(i)
              tstriny(j)=tstriny(i)
              tdmax(j)=tdmax(i)
              tmu_s(j)=tmu_s(i)
              tmu_d(j)=tmu_d(i)
            end if
          else
            if ((ifault==0).or.(ifault==1))  then
              do k=1,read_step
                taxx(j,k)=taxx(i,k)
                tayy(j,k)=tayy(i,k)
                tazz(j,k)=tazz(i,k)
                taxz(j,k)=taxz(i,k)
                tayz(j,k)=tayz(i,k)
                taxy(j,k)=taxy(i,k)
              end do
            endif
          endif      ! match if (ifault==3 .or. ifault==4)

C         TAG PROC AS FAULT BEARING
          srcproc = rank
          j=j+1
        end if
      end do

C     SET NUMBER OF FAULT NODES IN PROC
      npsrc=j-1

CCC   ALLOCATE SRC ARRAYS AND DEALLOCATE TMP ARRAYS
      if(cur_step == 1) then
        allocate(psrc(npsrc,maxdim))
      end if

      if (ifault==3 .or. ifault==4) then
        if(idyna == 0) then
          allocate(strinxx(npsrc))
          allocate(strisxx(npsrc))
          allocate(dmaxx(npsrc))
          strinxx=0
          strisxx=0
          dmaxx=0
        end if
      else    ! match if (ifault==3|4)
        if(cur_step == 1) then
          allocate(axx(npsrc,nst))
          allocate(ayy(npsrc,nst))
          allocate(azz(npsrc,nst))
          allocate(axz(npsrc,nst))
          allocate(ayz(npsrc,nst))
          allocate(axy(npsrc,nst))
        end if
      endif

C     TRANSFER SRC INFO TO PERMANENT ARRAYS

      do i=1,npsrc
        psrc(i,1) = tpsrc(i,1)
        psrc(i,2) = tpsrc(i,2)
        psrc(i,3) = tpsrc(i,3)
      end do

      fx_dt=0
      fz_dt=0
      if(srcproc == rank) then
        fx_dt=nxt
        fz_dt=nzt
      end if

      if (ifault==3 .or. ifault ==4) then
C       allocate a y_dt value to distinguish if node is in (+)/(-) side of the fault
C       it is used in subroutine inihomo to define bimaterial properties (fast initialization)
        if(idyna == 0)  then
           do i=1,npsrc
             strisxx(i) = tstrisxx(i)
             strinxx(i) = tstrinxx(i)
             dmaxx(i) = tdmaxx(i)
           end do
        else
          if(coords(2)==ipy) y_dt=mod(y_dt,nyt)
          if(coords(2).gt.ipy) y_dt=-1

C         allocate dynamic source array
          if (cur_step == 1) then
            allocate(strisx(-1:fx_dt+2,-1:fz_dt+2))
            allocate(striny(-1:fx_dt+2,-1:fz_dt+2))
            allocate(strisz(-1:fx_dt+2,-1:fz_dt+2))
            allocate(u_f(-1:fx_dt+2,4,-1:fz_dt+2))
            allocate(v_f(-1:fx_dt+2,4,-1:fz_dt+2))
            allocate(w_f(-1:fx_dt+2,4,-1:fz_dt+2))
            allocate(xx_f(-1:fx_dt+2,4,-1:fz_dt+2))
            allocate(yy_f(-1:fx_dt+2,4,-1:fz_dt+2))
            allocate(zz_f(-1:fx_dt+2,4,-1:fz_dt+2))
            allocate(xy_f(-1:fx_dt+2,4,-1:fz_dt+2))
            allocate(xz_f(-1:fx_dt+2,4,-1:fz_dt+2))
            allocate(yz_f(-1:fx_dt+2,4,-1:fz_dt+2))
            allocate(trupt(-1:fx_dt+2,-1:fz_dt+2))
            allocate(brokenu(-1:fx_dt+2,-1:fz_dt+2))
            allocate(brokenw(-1:fx_dt+2,-1:fz_dt+2))
            allocate(d_u(-1:fx_dt+2,-1:fz_dt+2))
            allocate(d_w(-1:fx_dt+2,-1:fz_dt+2))
            allocate(rateu_u(-1:fx_dt+2,-1:fz_dt+2))
            allocate(rateu_w(-1:fx_dt+2,-1:fz_dt+2))
            allocate(rateu(-1:fx_dt+2,-1:fz_dt+2))
            allocate(ratew_u(-1:fx_dt+2,-1:fz_dt+2))
            allocate(ratew_w(-1:fx_dt+2,-1:fz_dt+2))
            allocate(ratew(-1:fx_dt+2,-1:fz_dt+2))
            allocate(Rplusu(-1:fx_dt+2,-1:fz_dt+2))
            allocate(Rplusw(-1:fx_dt+2,-1:fz_dt+2))
            allocate(Rplusue(-1:fx_dt+2,-1:fz_dt+2))
            allocate(Rpluswe(-1:fx_dt+2,-1:fz_dt+2))
            allocate(Rminusu(-1:fx_dt+2,-1:fz_dt+2))
            allocate(Rminusw(-1:fx_dt+2,-1:fz_dt+2))
            allocate(Rminusue(-1:fx_dt+2,-1:fz_dt+2))
            allocate(Rminuswe(-1:fx_dt+2,-1:fz_dt+2))

            allocate(slipu_u(-1:fx_dt+2,-1:fz_dt+2))
            allocate(slipu_w(-1:fx_dt+2,-1:fz_dt+2))
            allocate(trv1(-1:fx_dt+2,-1:fz_dt+2))
            allocate(vplus(-1:fx_dt+2,-1:fz_dt+2))
            allocate(wplusi(-1:fx_dt+2,-1:fz_dt+2))
            allocate(umin(-1:fx_dt+2,-1:fz_dt+2))
            allocate(vmin(-1:fx_dt+2,-1:fz_dt+2))
            allocate(wmin(-1:fx_dt+2,-1:fz_dt+2))
            allocate(duplus(-1:fx_dt+2,-1:fz_dt+2))
            allocate(dvplus(-1:fx_dt+2,-1:fz_dt+2))
            allocate(dwplus(-1:fx_dt+2,-1:fz_dt+2))
            allocate(dumin(-1:fx_dt+2,-1:fz_dt+2))
            allocate(dvmin(-1:fx_dt+2,-1:fz_dt+2))
            allocate(dwmin(-1:fx_dt+2,-1:fz_dt+2))
          end if    ! match 'if (cur_step == 1)'

          slipu_u = 0.0
          slipu_w = 0.0
          trv1 = 0.0
          vplus = 0.0
          wplusi = 0.0
          umin = 0.0
          vmin = 0.0
          wmin = 0.0
          duplus = 0.0
          dvplus = 0.0
          dwplus = 0.0
          dumin = 0.0
          dvmin = 0.0
          dwmin = 0.0

          strisx=0.
          striny=120.0e6
          strisz=0.
          u_f=0.
          v_f=0.
          w_f=0.
          xx_f=0.
          yy_f=0.
          zz_f=0.
          xy_f=0.
          xz_f=0.
          yz_f=0.
          trupt=0.
          brokenu=0
          brokenw=0
          d_u=0.
          d_w=0.
          rateu_u=0.
          rateu_w=0.
          rateu=0.
          ratew_u=0.
          ratew_w=0.
          ratew=0.
          Rplusu=0.
          Rplusw=0.
          Rplusue=0.
          Rpluswe=0.
          Rminusu=0.
          Rminusw=0.
          Rminusue=0.
          Rminuswe=0.

C         and set to overlap two nodes between procs
          if(cur_step == 1) then
            allocate(d0(-1:fx_dt+2,-1:fz_dt+2))
            allocate(mu_s(-1:fx_dt+2,-1:fz_dt+2))
            allocate(mu_d(-1:fx_dt+2,-1:fz_dt+2))
            allocate(trw1(-1:fx_dt+2,-1:fz_dt+2))
            allocate(tru1(-1:fx_dt+2,-1:fz_dt+2))
            allocate(tru2(-1:fx_dt+2,-1:fz_dt+2))
            allocate(trw2(-1:fx_dt+2,-1:fz_dt+2))
            allocate(uplus(-1:fx_dt+2,-1:fz_dt+2))
            allocate(wplus(-1:fx_dt+2,-1:fz_dt+2))
            allocate(xxplus(-1:fx_dt+2,-1:fz_dt+2))
            allocate(zzplus(-1:fx_dt+2,-1:fz_dt+2))
            allocate(xzplus(-1:fx_dt+2,-1:fz_dt+2))
          end if

          d0=1.0
          mu_s=1000.
          mu_d=1000.
          trw1=0.
          tru1=0.
          trw2=0.
          tru2=0.
          uplus=0.
          wplus=0.
          xxplus=0.
          zzplus=0.
          xzplus=0.

          do i=1,npsrc
            strisx(psrc(i,1),psrc(i,3)) = tstrisx(i)
            tru1(psrc(i,1),psrc(i,3)) = tstrisx(i)
            tru2(psrc(i,1),psrc(i,3)) = tstrisx(i)
            strisz(psrc(i,1),psrc(i,3)) = tstrisz(i)
            trw1(psrc(i,1),psrc(i,3)) = tstrisz(i)
            trw2(psrc(i,1),psrc(i,3)) = tstrisz(i)
            striny(psrc(i,1),psrc(i,3)) = tstriny(i)
            d0(psrc(i,1),psrc(i,3)) = tdmax(i)
            mu_s(psrc(i,1),psrc(i,3)) = tmu_s(i)
            mu_d(psrc(i,1),psrc(i,3)) = tmu_d(i)
          end do
        end if     ! match 'if(idyna == 0)'
      else if (ifault.ne.2)  then
        do i=1,npsrc
          do j=1,read_step
            axx(i,j+cur_step-1) = taxx(i,j)
            ayy(i,j+cur_step-1) = tayy(i,j)
            azz(i,j+cur_step-1) = tazz(i,j)
            axz(i,j+cur_step-1) = taxz(i,j)
            ayz(i,j+cur_step-1) = tayz(i,j)
            axy(i,j+cur_step-1) = taxy(i,j)
          end do
        end do
      end if

C     DEALLOCATE TEMPORAY SOURCE ARRAYS
      if(cur_step + read_step > nst) then
        deallocate(tpsrc)
      end if

      if (ifault==3 .or. ifault==4) then
        if(idyna == 0)  then
          deallocate(tstrisxx)
          deallocate(tstrinxx)
          deallocate(tdmaxx)
        else
C ==ADDED FOR SGSN DYNAMIC MODEL
          deallocate(tstrisx)
          deallocate(tstriny)
          deallocate(tstrisz)
          deallocate(tdmax)
          deallocate(tmu_s)
          deallocate(tmu_d)
        end if
C ==
      else if (ifault.ne.2) then
        if(cur_step + read_step > nst) then
          deallocate(taxx)
          deallocate(tayy)
          deallocate(tazz)
          deallocate(taxz)
          deallocate(tayz)
          deallocate(taxy)
        end if
      endif

      ! for ifault 1, create partitioned fault file but this is only for small files
      ! For bigger files, SourcePartitioner should be used and ifault=2 should be set
      if ((ifault==1).and.(npsrc.ne.0))  then
        write( faultfile, '(a,i7.7)') 
     $         'input_rst/srcpart/part_faults/fault', rank
        if (rank==0) print *,'create subdomain fault files '
        open( 19, file=faultfile,form='unformatted', access='direct',
     $       recl=bsize*npsrc*nst*3)
        write(19,rec=1) axx,ayy,azz
        write(19,rec=2) axz,ayz,axy
        close(19)
      end if

      if ((ifault==2).and.(npsrc.ne.0))  then
        write( faultfile, '(a,i7.7)') 
     $         'input_rst/srcpart/part_faults/fault', rank
        if (rank==0)  print *,'use subdomain fault files '
        open( 19, file=faultfile, form='unformatted', access='direct',
     $        recl=bsize*npsrc*nst*3, status='old')
        read( 19, rec=1 ) axx,ayy,azz
        read( 19, rec=2 ) axz,ayz,axy
        close( 19 )
      end if

      return
      end

      subroutine read_src_ifault_2(coords,maxdim,comm,rank,master,nxt,nyt,nzt,srcproc,
     +       nsrc,nst,npsrc,nz,idx,read_step)

      use parstat
      implicit none

      integer, parameter :: filen=180
      integer :: maxdim,comm,rank,master,nxt,nyt,nzt,srcproc
      integer :: nsrc,nst,npsrc,nz,idx, read_step
      integer, dimension(maxdim) :: coords

      integer :: nxp,nyp,nzp,i,j,err
      integer :: varnum(1),bsize
      character (len=filen) :: faultfile, tpsrcfile

      real :: dummy

      bsize=4*size(varnum)
      inquire(iolength=bsize)varnum

      ! load individual tpsrc files to each rank
      if (idx==1) then
C.1   Each process READ tpsrc
        write(tpsrcfile, '(a,i7.7)') 'input_rst/srcpart/tpsrc/tpsrc',rank
        open(17, file=tpsrcfile,form='unformatted', access='direct',recl=bsize*3,
     +      status='old', iostat=err)
        if (err.eq.0) then ! no error
          read(17, rec=1) npsrc, dummy, dummy
          print *, "Reading ", tpsrcfile, '-', npsrc
          do i=1,npsrc
            read(17, rec=i+1) tpsrc(i,1), tpsrc(i,2), tpsrc(i,3)
            tpsrc(i,3) = nz + 1 - tpsrc(i,3)
          end do
          close(17)
        else
          ! no tpsrc_tmp files
          npsrc=0
        end if

C.2   DETERMINE IF NODE BELONGS TO PROC
        j=1
        do i=1,npsrc
          nxp = int((tpsrc(i,1)-1)/nxt)
          nyp = int((tpsrc(i,2)-1)/nyt)
          nzp = int((tpsrc(i,3)-1)/nzt)

          !if(coords(1)==nxp .and. coords(2)==nyp .and. coords(3)==nzp) then
            if(mod(tpsrc(i,1),nxt)/=0) then
              tpsrc(j,1) = mod(tpsrc(i,1),nxt)
            else
              tpsrc(j,1) = nxt
            end if

            if(mod(tpsrc(i,2),nyt)/=0) then
              tpsrc(j,2) = mod(tpsrc(i,2),nyt)
            else
              tpsrc(j,2) = nyt
            end if

            if(mod(tpsrc(i,3),nzt)/=0) then
              tpsrc(j,3) = mod(tpsrc(i,3),nzt)
            else
              tpsrc(j,3) = nzt
            end if

            srcproc = rank
            j = j+1
          !end if
        end do

        !npsrc = j-1    ! SET NUMBER OF FAULT NODES IN PROC
        !print *, 'rank, npsrc: ', rank, npsrc

        if (npsrc.ne.0) then
          allocate(psrc(npsrc,maxdim))

          allocate(axx(npsrc,read_step))
          allocate(ayy(npsrc,read_step))
          allocate(azz(npsrc,read_step))
          allocate(axz(npsrc,read_step))
          allocate(ayz(npsrc,read_step))
          allocate(axy(npsrc,read_step))

C.3     TRANSFER SRC INFO TO PERMANENT ARRAYS     ! ?? not needed actually for ifault==2 ?
          do i=1,npsrc
            psrc(i,1) = tpsrc(i,1)
            psrc(i,2) = tpsrc(i,2)
            psrc(i,3) = tpsrc(i,3)
          end do

          deallocate(tpsrc)    ! deallocate tpsrc() now
        end if
      end if    ! match 'if (idx==1)'

!!    Read axx,...axy from corresponding splitted sub-grid fault file

      if (npsrc.ne.0)  then
        write(faultfile, '(a,i7.7,a,i3.3)') 
     +       'input_rst/srcpart/split_faults/fault', rank, '_', idx
        if (rank==master)  print *,'Read source from ', faultfile,bsize,npsrc,read_step
        open(19, file=faultfile, form='unformatted', access='direct',
     +         recl=bsize*npsrc*read_step*3, status='old')
        read(19, rec=1) axx,ayy,azz
        read(19, rec=2) axz,ayz,axy
        close(19)
      end if

      return
      end

      subroutine addsrc(i,dx,dt,nst,npsrc,read_step,ifault,rank,igreen)

CCC   ROUTINE TO ADD STAGGERED-GRID SOURCE TO STRESS TENSOR
      use parstat
      implicit none

      integer :: i,nst,npsrc,read_step,ifault,rank,igreen
      real :: dt,dx,vtst
      integer :: j,k

      if(i .le. nst) then
C       MARCH OVER NUMBER OF FAULT NODES  IN PROC
        if (ifault==2) then
          if (mod(i,read_step) .ne. 0) then
            k = mod(i,read_step)
          else
            k = read_step
          end if
        else
          k=i
        end if

        vtst = dt/(dx*dx*dx)     

        do j=1,npsrc

          if (igreen.eq.-1) then

          xx(psrc(j,1),psrc(j,2),psrc(j,3)) =
     +    xx(psrc(j,1),psrc(j,2),psrc(j,3)) - vtst*axx(j,k)

          yy(psrc(j,1),psrc(j,2),psrc(j,3)) =
     +    yy(psrc(j,1),psrc(j,2),psrc(j,3)) - vtst*ayy(j,k)

          zz(psrc(j,1),psrc(j,2),psrc(j,3)) =
     +    zz(psrc(j,1),psrc(j,2),psrc(j,3)) - vtst*azz(j,k)

          xz(psrc(j,1),psrc(j,2),psrc(j,3)) =
     +    xz(psrc(j,1),psrc(j,2),psrc(j,3)) - vtst*axz(j,k)

          yz(psrc(j,1),psrc(j,2),psrc(j,3)) =
     +    yz(psrc(j,1),psrc(j,2),psrc(j,3)) - vtst*ayz(j,k)

          xy(psrc(j,1),psrc(j,2),psrc(j,3)) =
     +    xy(psrc(j,1),psrc(j,2),psrc(j,3)) - vtst*axy(j,k)

          elseif (igreen.eq.1) then

          u1(psrc(j,1),psrc(j,2),psrc(j,3)) =
     +    u1(psrc(j,1),psrc(j,2),psrc(j,3)) + axx(j,k)*vtst/d1(psrc(j,1),psrc(j,2),psrc(j,3))

          elseif (igreen.eq.2) then

          v1(psrc(j,1),psrc(j,2),psrc(j,3)) =
     +    v1(psrc(j,1),psrc(j,2),psrc(j,3)) + ayy(j,k)*vtst/d1(psrc(j,1),psrc(j,2),psrc(j,3))

          elseif (igreen.eq.3) then

          w1(psrc(j,1),psrc(j,2),psrc(j,3)) =
     +    w1(psrc(j,1),psrc(j,2),psrc(j,3)) + azz(j,k)*vtst/d1(psrc(j,1),psrc(j,2),psrc(j,3))

          elseif (igreen.eq.-2) then

          u1(psrc(j,1),psrc(j,2),psrc(j,3)) =
     +    u1(psrc(j,1),psrc(j,2),psrc(j,3)) + axx(j,k)*vtst/d1(psrc(j,1),psrc(j,2),psrc(j,3))

          v1(psrc(j,1),psrc(j,2),psrc(j,3)) =
     +    v1(psrc(j,1),psrc(j,2),psrc(j,3)) + ayy(j,k)*vtst/d1(psrc(j,1),psrc(j,2),psrc(j,3))

          w1(psrc(j,1),psrc(j,2),psrc(j,3)) =
     +    w1(psrc(j,1),psrc(j,2),psrc(j,3)) + azz(j,k)*vtst/d1(psrc(j,1),psrc(j,2),psrc(j,3))

          end if
        end do
      end if

      return
      end

      subroutine dyna(i,mu_dd,mu_ss,dh,dt,npsrc,nxt,nyt,nzt,
     +                nsrc,nst,rank)

c     subroutine to compute dynamic slip and sliprate on
c     a vertical planar fault surface located at Y=nbgx=psrc(1:npsrc,2).
c     Should be called during each time step in FD code

c     i     current time step in FD code (int)
c     npsrc number of subfaults (int)

      use parstat

      real :: mu_dd,mu_ss,tmpu2

      integer :: rank,itmpbroken

      do j=1,npsrc
        deltad=0.0

        tmpu2=u2(psrc(j,1),psrc(j,2),psrc(j,3))
        itmpbroken=broken(psrc(j,1),psrc(j,2),psrc(j,3))

        tot1=xy(psrc(j,1),psrc(j,2),psrc(j,3))+strisxx(j)
        tot=strinxx(j)*mu_ss

        if (tmpu2.le.dmaxx(j).and.dmaxx(j).gt.0.0) then
          tot=strinxx(j)*((1.-tmpu2/dmaxx(j))*mu_ss+mu_dd*tmpu2/dmaxx(j))
        else
          tot=strinxx(j)*mu_dd
        endif

        if (tot1.ge.tot) then
          broken(psrc(j,1),psrc(j,2),psrc(j,3))=1
          rupt(psrc(j,1),psrc(j,2),psrc(j,3))=i*dt*(1-itmpbroken)
     +                                        + rupt(psrc(j,1),psrc(j,2),psrc(j,3))*itmpbroken
          xm1 = mu (psrc(j,1),psrc(j,2),psrc(j,3))
          xm2 = mu (psrc(j,1),psrc(j,2)+1,psrc(j,3))
          xm3 = mu (psrc(j,1)-1,psrc(j,2)+1,psrc(j,3))
          xm4 = mu (psrc(j,1)-1,psrc(j,2),psrc(j,3))
          xm5 = mu (psrc(j,1),psrc(j,2),psrc(j,3)-1)
          xm6 = mu (psrc(j,1),psrc(j,2)+1,psrc(j,3)-1)
          xm7 = mu (psrc(j,1)-1,psrc(j,2)+1,psrc(j,3)-1)
          xm8 = mu (psrc(j,1)-1,psrc(j,2),psrc(j,3)-1)
          xm = 0.125*(xm1+xm2+xm3+xm4+xm5+xm6+xm7+xm8)
          strainp=(tot1-tot)/xm
          deltad=strainp*dh
          u2(psrc(j,1),psrc(j,2),psrc(j,3))=u2(psrc(j,1),psrc(j,2),psrc(j,3))+deltad
          xy(psrc(j,1),psrc(j,2),psrc(j,3))=tot-strisxx(j)
        endif

        u3(psrc(j,1),psrc(j,2),psrc(j,3))=deltad/dt
      enddo

      return
      end

