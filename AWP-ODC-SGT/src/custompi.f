CCC     'custompi.f' CONTAINS CUSTOM MPI-LIKE ROUTINES (CMPI)


      subroutine CMPI_CART_JUMP(comm,coords,dims,maxdim,rank,shew,shns,shud,newrank)

CCC   FINDS RANK OF A NEARBY (ADJACENT, EDGE, CORNER) LOCATED PROC

      include "mpif.h"

      integer :: comm,rank,newrank
      integer :: shew,shns,shud

C     define following variables  

      integer :: maxdim
      integer :: ierr
       
      integer, dimension(maxdim) :: coords,dims,newcoords

      newcoords(1) = coords(1) + shew
      newcoords(2) = coords(2) + shns
      newcoords(3) = coords(3) + shud

      newrank = MPI_PROC_NULL

C     DETERMINE IF NEW COORDINATES ARE IN MODEL SPACE

      if( (newcoords(1).lt.dims(1) .and. newcoords(1).gt.-1) .and.
     +    (newcoords(2).lt.dims(2) .and. newcoords(2).gt.-1) .and.
     +    (newcoords(3).lt.dims(3) .and. newcoords(3).gt.-1) ) 
     +      call MPI_CART_RANK(comm,newcoords,newrank,ierr)

      return
      end

 
