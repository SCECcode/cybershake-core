CCC   'station.f' RECORDS SEISMOGRAMS AND WAVEFIELDS


      subroutine getrec(igreen,coords,maxdim,rank,nxt,nyt,nzt,nz,recproc,nprec,nrec,
     +           sgt_recproc,sgt_nprec,
     +           nbgx,nedx,nskpx,nbgy,nedy,nskpy,nbgz,nedz,nskpz,casenumber)

CCC   DETERMINES STATION DISTRIBUTION AMONG PROCS

      use parstat
      integer :: igreen, i,j,k
      integer*8 :: r
      integer :: l
      integer :: nprec,sgt_nprec
      integer(i8) :: nrec
      integer :: nxt,nyt,nzt, nz, rank,recproc, sgt_recproc, casenumber
      integer :: nbgx,nedx,nskpx,nbgy,nedy,nskpy,nbgz,nedz,nskpz
      integer, dimension(maxdim) :: coords
      integer, dimension(:,:), allocatable :: tpsgt_local_sta_coord
      integer*8 :: Li, Lj, Lk
      integer*8 :: nbgx_l,nedx_l,nskpx_l,nbgy_l,nedy_l,nskpy_l,nbgz_l,nedz_l,nskpz_l

C     SEARCH OVER ALL RECEIVER NODES

C     DETERMINE IF NODE BELONGS TO PROC AND COMPUTE FILE POSITION

cSGT -------start---------------------------------------
c     SGT's extended source output mode
c     In this mode, only selected stations' SGT tensors are output.
c     Determine if station belongs to proc and compute its position
      if (sgt_io_out .eq. 1) then
c        if (rank==0) write (*,*) 'getrec start'
        allocate(tpsgt_local_sta_coord(nxt*nyt*nzt,4))
        l=1;
        do i=1,sgt_numsta
            nzp = int((nz-sgt_sta_coord(i,3))/nzt)
            if (coords(3)==nzp) then
                nyp = int((sgt_sta_coord(i,2)-1)/nyt)
                if (coords(2)==nyp) then
                    nxp = int((sgt_sta_coord(i,1)-1)/nxt)
                    if (coords(1)==nxp) then
                        a = mod(sgt_sta_coord(i,1),nxt)
                        b = ANINT(a/(a+1))
                        tprec(l,1) = a + (1-b) * nxt

                        a = mod(sgt_sta_coord(i,2),nyt)
                        b = ANINT(a/(a+1))
                        tprec(l,2) = a + (1-b) * nyt

                        a = mod(nz-sgt_sta_coord(i,3)+1,nzt)
                        b = ANINT(a/(a+1))
                        tprec(l,3) = a + (1-b) * nzt

                        tpsgt_local_sta_coord(l,1)=sgt_sta_coord(i,1)
                        tpsgt_local_sta_coord(l,2)=sgt_sta_coord(i,2)
                        tpsgt_local_sta_coord(l,3)=sgt_sta_coord(i,3)
                        tpsgt_local_sta_coord(l,4)=i-1  !global index
                        l=l+1
                        sgt_recproc=rank
                    endif
                endif
            endif
        end do

        sgt_nprec=l-1;
        if (casenumber .eq. 1) then
            allocate(sgt_prec(sgt_nprec,maxdim))
            allocate(sgt_local_sta_coord(sgt_nprec,4))
            do i=1,sgt_nprec
                sgt_prec(i,1) = tprec(i,1)
                sgt_prec(i,2) = tprec(i,2)
                sgt_prec(i,3) = tprec(i,3)
                sgt_local_sta_coord(i,1)=tpsgt_local_sta_coord(i,1)
                sgt_local_sta_coord(i,2)=tpsgt_local_sta_coord(i,2)
                sgt_local_sta_coord(i,3)=tpsgt_local_sta_coord(i,3)
                sgt_local_sta_coord(i,4)=tpsgt_local_sta_coord(i,4) !global index
            end do
        endif
        deallocate(tpsgt_local_sta_coord)
      endif
cSGT --------end----------------------------------------


      l = 1
      r = -1

      nbgx_l=nbgx
      nedx_l=nedx
      nskpx_l=nskpx
      nbgy_l=nbgy
      nedy_l=nedy
      nskpy_l=nskpy
      nbgz_l=nbgz
      nedz_l=nedz
      nskpz_l=nskpz

      do k=nbgz,nedz,nskpz
        if(k==0) then
           nzp = int((nz-1)/nzt);
        else
           nzp = int((nz-k)/nzt);
        end if

        if (coords(3)==nzp) then
          do j=nbgy,nedy,nskpy
            nyp = int((j-1)/nyt)

            if (coords(2)==nyp) then
              do i=nbgx,nedx,nskpx
                nxp = int((i-1)/nxt)

                if (coords(1)==nxp) then
                  a = mod(i,nxt)
                  b = ANINT(a/(a+1))
                  tprec(l,1) = a + (1-b) * nxt

                  a = mod(j,nyt)
                  b = ANINT(a/(a+1))
                  tprec(l,2) = a + (1-b) * nyt

                  if(k == 0) then
                     tprec(l,3) = nzt + 1;
                  else
                     a = mod(nz-k+1,nzt)
                     b = ANINT(a/(a+1))
                     tprec(l,3) = a + (1-b) * nzt
                  end if

C    TAG PROC AS STATION BEARING
c    Calculate r first,
                  Li = int((nedx-nbgx)/nskpx) + 1
                  Lj = int((nedy-nbgy)/nskpy) + 1
                  Lk = int((nedz-nbgz)/nskpz) + 1
                  r = Li*Lj*int((k-nbgz)/nskpz) + Li*int((j-nbgy)/nskpy) + int((i-nbgx)/nskpx)

                  tpmap(l) = r
                  recproc = rank
                  l=l+1
                endif
              end do
            endif
          end do
        endif
      end do


C     SET NUMBER OF STATION NODES IN PROC

      nprec=l-1
      nrec=(((nedx_l-nbgx_l)/nskpx_l)+1) * (((nedy_l-nbgy_l)/nskpy_l)+1) * (((nedz_l-nbgz_l)/nskpz_l)+1)

C     ALLOCATE RECEIVER ARRAYS AND DEALLOCATE TMP ARRAYS
c      allocate(prec(nprec,maxdim))
c      allocate(pmap(nprec))

C     TRANSFER RECEIVER INFO TO PERMANENT ARRAYS
c      mm=1
CC     use casenumber to assign to different prec* and pmap* variables
CC     1: ; 2: original getrec2; 3: for rupt output; 4: for all other 19 variables in SGSN

      SELECT case (casenumber)
        case (1)
          allocate(prec(nprec,maxdim))
          allocate(pmap(nprec))
    
          do i=1,nprec
            prec(i,1) = tprec(i,1)
            prec(i,2) = tprec(i,2)
            prec(i,3) = tprec(i,3)

            pmap(i) = tpmap(i)
          end do

        case (2)
          allocate(prec2(nprec,maxdim))
          allocate(pmap2(nprec))

         do i=1,nprec
            prec2(i,1) = tprec(i,1)
            prec2(i,2) = tprec(i,2)
            prec2(i,3) = tprec(i,3)

            pmap2(i) = tpmap(i)
          end do

        case (3)
          allocate(prec3(nprec,maxdim))
          allocate(pmap3(nprec))

         do i=1,nprec
            prec3(i,1) = tprec(i,1)
            prec3(i,2) = tprec(i,2)
            prec3(i,3) = tprec(i,3)

            pmap3(i) = tpmap(i)
         end do

      END SELECT

C     DEALLOCATE TEMPORARY RECEIVER ARRAYS

        deallocate(tprec)
        deallocate(tpmap)
 
        return
      end

c     wrtrec is not in use
      subroutine wrtrec(ii,ntiskp,fhx,fhy,fhz,nprec,nrec)
             
CCC   WRITE OUT SEISMOGRAMS IN SPECIFIED VOLUME
      use parstat

      !include 'mpif.h'
c      integer,parameter :: rtype=selected_int_kind( 18 )
      integer :: nprec
      integer :: nrec
c      integer(i8) :: nrec
c      integer(rtype) :: nprec,nrec
      integer :: ii,j, ntiskp,it
      integer :: err,filetype
      integer :: xn,yn,zn
      integer :: fhx,fhy,fhz

      integer (kind=MPI_OFFSET_KIND) :: disp
      integer :: status(MPI_STATUS_SIZE)
      real, dimension(nprec) :: bufx,bufy,bufz

      disp=0

      it = ntiskp*(ii/ntiskp)
      if(ii==it) then
        call MPI_TYPE_CREATE_INDEXED_BLOCK(nprec,1,pmap,MPI_REAL,
     +                                   filetype,err)

        call MPI_TYPE_COMMIT(filetype,err) 

        call MPI_FILE_SET_VIEW(fhx,disp,MPI_REAL,filetype,'native',
     +                       MPI_INFO_NULL,err)

        call MPI_FILE_SET_VIEW(fhy,disp,MPI_REAL,filetype,'native',
     +                       MPI_INFO_NULL,err)

        call MPI_FILE_SET_VIEW(fhz,disp,MPI_REAL,filetype,'native',
     +                       MPI_INFO_NULL,err)

        do j=1,nprec

          xn = prec(j,1)
          yn = prec(j,2)
          zn = prec(j,3)

          bufx(j) = u1(xn,yn,zn)
          bufy(j) = v1(xn,yn,zn)
          bufz(j) = w1(xn,yn,zn)


        end do

        call MPI_FILE_WRITE_ALL(fhx,bufx,nprec,MPI_REAL,status,err)
        call MPI_FILE_WRITE_ALL(fhy,bufy,nprec,MPI_REAL,status,err)
        call MPI_FILE_WRITE_ALL(fhz,bufz,nprec,MPI_REAL,status,err)

        pmap = pmap + nrec
    
      end if
      return
      end

