CCC   'bound.f' SETS BOUNDS FOR EACH COMPUTATIONAL ZONE FOR EACH PROC

C     PROPAGATION SPACE COMPRISED OF 18 ZONES (17 PML AND 1 INTERNAL)
C     WITH EACH ZONE DESCRIBED BY 3 (6 POSSIBLE) COORDINATE RANGES

C     RANGE 1 = EDGE TO PML INTERFACE
C     RANGE 2 = PML INTERFACE TO PROC BOUNDARY
C     RANGE 3 = FULL PROC RANGE BETWEEN PML ZONES 
C     RANGE 4 = PROC BOUNDARY TO PML INTERFACE
C     RANGE 5 = PML INTERFACE TO EDGE
C     RANGE 6 = PML INTERFACE TO PML INTERFACE


      subroutine x1(i,nx,nd)

CCC   SETS X BOUNDS FOR X BOUND RANGE 1

      use parstat
  
      integer :: i,nx,nd

      xi(i) = 1
      xf(i) = nd

      return
      end


      subroutine x2(i,nx,nd)

CCC   SETS X BOUNDS FOR X BOUND RANGE 2

      use parstat
  
      integer :: i,nx,nd

      xi(i) = nd+1
      xf(i) = nx

      return
      end


      subroutine x3(i,nx,nd)

CCC   SETS X BOUNDS FOR X BOUND RANGE 3

      use parstat
  
      integer :: i,nx,nd

      xi(i) = 1
      xf(i) = nx

      return
      end


      subroutine x4(i,nx,nd)

CCC   SETS X BOUNDS FOR X BOUND RANGE 4

      use parstat
  
      integer :: i,nx,nd

      xi(i) = 1
      xf(i) = nx-nd

      return
      end


      subroutine x5(i,nx,nd)

CCC   SETS X BOUNDS FOR X BOUND RANGE 5

      use parstat
  
      integer :: i,nx,nd

      xi(i) = nx-nd+1
      xf(i) = nx

      return
      end


      subroutine x6(i,nx,nd)

CCC   SETS X BOUNDS FOR X BOUND RANGE 6

      use parstat
  
      integer :: i,nx,nd

      xi(i) = nd+1 
      xf(i) = nx-nd

      return
      end


      subroutine y1(i,ny,nd)

CCC   SETS Y BOUNDS FOR Y BOUND RANGE 1

      use parstat
  
      integer :: i,ny,nd

      yi(i) = 1
      yf(i) = nd

      return
      end


      subroutine y2(i,ny,nd)

CCC   SETS Y BOUNDS FOR Y BOUND RANGE 2

      use parstat
  
      integer :: i,ny,nd

      yi(i) = nd+1
      yf(i) = ny

      return
      end


      subroutine y3(i,ny,nd)

CCC   SETS Y BOUNDS FOR Y BOUND RANGE 3

      use parstat
  
      integer :: i,ny,nd

      yi(i) = 1
      yf(i) = ny
 
      return
      end


      subroutine y4(i,ny,nd)

CCC   SETS Y BOUNDS FOR Y BOUND RANGE 4

      use parstat
  
      integer :: i,ny,nd

      yi(i) = 1
      yf(i) = ny-nd

      return
      end


      subroutine y5(i,ny,nd)

CCC   SETS Y BOUNDS FOR Y BOUND RANGE 5

      use parstat

      integer :: i,ny,nd

      yi(i) = ny-nd+1
      yf(i) = ny

      return
      end


      subroutine y6(i,ny,nd)

CCC   SETS Y BOUNDS FOR Y BOUND RANGE 6

      use parstat
  
      integer :: i,ny,nd

      yi(i) = nd+1
      yf(i) = ny-nd

      return
      end


      subroutine z1(i,nz,nd)

CCC   SETS Z BOUNDS FOR Z BOUND RANGE 1

      use parstat
  
      integer :: i,nz,nd

      zi(i) = 1
      zf(i) = nd

      return
      end


      subroutine z2(i,nz,nd)

CCC   SETS Z BOUNDS FOR Z BOUND RANGE 2

      use parstat
  
      integer :: i,nz,nd

      zi(i) = nd+1
      zf(i) = nz

      return
      end


      subroutine z3(i,nz,nd)

CCC   SETS Z BOUNDS FOR Z BOUND RANGE 3

      use parstat
  
      integer :: i,nz,nd

      zi(i) = 1
      zf(i) = nz

      return
      end



 
      subroutine bounds(nx,ny,nz,nd,fs,fw,fn,fe,fb)

CCC   SET BOUNDS FOR EACH REGION AND PROC

      integer :: fs,fw,fn,fe,fb, nx,ny,nz,nd

C     PLANE REGION PROCS

      if     (fs==1 .and. fw==0 .and. fn==0 .and. fe==0 .and. fb==0) then
        call x3(1,nx,nd)
        call y1(1,ny,nd)
        call z3(1,nz,nd) 
        call x3(18,nx,nd)
        call y2(18,ny,nd)
        call z3(18,nz,nd)
      
      else if(fs==0 .and. fw==1 .and. fn==0 .and. fe==0 .and. fb==0) then
        call x1(2,nx,nd) 
        call y3(2,ny,nd)
        call z3(2,nz,nd)
        call x2(18,nx,nd)
        call y3(18,ny,nd)
        call z3(18,nz,nd)

      else if(fs==0 .and. fw==0 .and. fn==1 .and. fe==0 .and. fb==0) then
        call x3(3,nx,nd) 
        call y5(3,ny,nd)
        call z3(3,nz,nd)
        call x3(18,nx,nd)
        call y4(18,ny,nd)
        call z3(18,nz,nd)

      else if(fs==0 .and. fw==0 .and. fn==0 .and. fe==1 .and. fb==0) then
        call x5(4,nx,nd) 
        call y3(4,ny,nd)
        call z3(4,nz,nd)
        call x4(18,nx,nd)
        call y3(18,ny,nd)
        call z3(18,nz,nd)

      else if(fs==0 .and. fw==0 .and. fn==0 .and. fe==0 .and. fb==1) then
        call x3(5,nx,nd) 
        call y3(5,ny,nd)
        call z1(5,nz,nd)
        call x3(18,nx,nd)
        call y3(18,ny,nd)
        call z2(18,nz,nd)

C     VERTICAL EDGE PROCS

      else if(fs==1 .and. fw==1 .and. fn==0 .and. fe==0 .and. fb==0) then
        call x1(6,nx,nd)
        call y1(6,ny,nd)
        call z3(6,nz,nd) 
        call x2(1,nx,nd) 
        call y1(1,ny,nd)
        call z3(1,nz,nd)
        call x1(2,nx,nd) 
        call y2(2,ny,nd)
        call z3(2,nz,nd)
        call x2(18,nx,nd)
        call y2(18,ny,nd)
        call z3(18,nz,nd)

      else if(fs==0 .and. fw==1 .and. fn==1 .and. fe==0 .and. fb==0) then
        call x1(7,nx,nd) 
        call y5(7,ny,nd)
        call z3(7,nz,nd) 
        call x1(2,nx,nd) 
        call y4(2,ny,nd)
        call z3(2,nz,nd)
        call x2(3,nx,nd) 
        call y5(3,ny,nd)
        call z3(3,nz,nd)
        call x2(18,nx,nd)
        call y4(18,ny,nd)
        call z3(18,nz,nd)

      else if(fs==1 .and. fw==0 .and. fn==0 .and. fe==1 .and. fb==0) then
        call x5(8,nx,nd) 
        call y1(8,ny,nd)
        call z3(8,nz,nd) 
        call x4(1,nx,nd) 
        call y1(1,ny,nd)
        call z3(1,nz,nd)
        call x5(4,nx,nd) 
        call y2(4,ny,nd)
        call z3(4,nz,nd)
        call x4(18,nx,nd)
        call y2(18,ny,nd)
        call z3(18,nz,nd)

      else if(fs==0 .and. fw==0 .and. fn==1 .and. fe==1 .and. fb==0) then
        call x5(9,nx,nd) 
        call y5(9,ny,nd)
        call z3(9,nz,nd) 
        call x4(3,nx,nd) 
        call y5(3,ny,nd)
        call z3(3,nz,nd)
        call x5(4,nx,nd) 
        call y4(4,ny,nd)
        call z3(4,nz,nd)
        call x4(18,nx,nd)
        call y4(18,ny,nd)
        call z3(18,nz,nd)

C     HORIZONTAL EDGE PROCS

      else if(fs==1 .and. fw==0 .and. fn==0 .and. fe==0 .and. fb==1) then
        call x3(10,nx,nd)
        call y1(10,ny,nd)
        call z1(10,nz,nd)
        call x3(1,nx,nd)
        call y1(1,ny,nd)
        call z2(1,nz,nd)
        call x3(5,nx,nd)
        call y2(5,ny,nd)
        call z1(5,nz,nd)
        call x3(18,nx,nd)
        call y2(18,ny,nd)
        call z2(18,nz,nd)

      else if(fs==0 .and. fw==1 .and. fn==0 .and. fe==0 .and. fb==1) then
        call x1(11,nx,nd)
        call y3(11,ny,nd)
        call z1(11,nz,nd)
        call x1(2,nx,nd)
        call y3(2,ny,nd)
        call z2(2,nz,nd)
        call x2(5,nx,nd)
        call y3(5,ny,nd)
        call z1(5,nz,nd)
        call x2(18,nx,nd)
        call y3(18,ny,nd)
        call z2(18,nz,nd)

      else if(fs==0 .and. fw==0 .and. fn==1 .and. fe==0 .and. fb==1) then
        call x3(12,nx,nd)
        call y5(12,ny,nd)
        call z1(12,nz,nd)
        call x3(3,nx,nd)
        call y5(3,ny,nd)
        call z2(3,nz,nd)
        call x3(5,nx,nd)
        call y4(5,ny,nd)
        call z1(5,nz,nd)
        call x3(18,nx,nd)
        call y4(18,ny,nd)
        call z2(18,nz,nd)

      else if(fs==0 .and. fw==0 .and. fn==0 .and. fe==1 .and. fb==1) then
        call x5(13,nx,nd)
        call y3(13,ny,nd)
        call z1(13,nz,nd)
        call x5(4,nx,nd)
        call y3(4,ny,nd)
        call z2(4,nz,nd)
        call x4(5,nx,nd)
        call y3(5,ny,nd)
        call z1(5,nz,nd)
        call x4(18,nx,nd)
        call y3(18,ny,nd)
        call z2(18,nz,nd)

C     CORNER PROCS

      else if(fs==1 .and. fw==1 .and. fn==0 .and. fe==0 .and. fb==1) then
        call x1(14,nx,nd)
        call y1(14,ny,nd)
        call z1(14,nz,nd)
        call x1(6,nx,nd)
        call y1(6,ny,nd)
        call z2(6,nz,nd)
        call x2(10,nx,nd)
        call y1(10,ny,nd)
        call z1(10,nz,nd)
        call x1(11,nx,nd)
        call y2(11,ny,nd)
        call z1(11,nz,nd)
        call x2(1,nx,nd)
        call y1(1,ny,nd)
        call z2(1,nz,nd)
        call x1(2,nx,nd)
        call y2(2,ny,nd)
        call z2(2,nz,nd)
        call x2(5,nx,nd)
        call y2(5,ny,nd)
        call z1(5,nz,nd)
        call x2(18,nx,nd)
        call y2(18,ny,nd)
        call z2(18,nz,nd)

      else if(fs==0 .and. fw==1 .and. fn==1 .and. fe==0 .and. fb==1) then
        call x1(15,nx,nd)
        call y5(15,ny,nd)
        call z1(15,nz,nd)
        call x1(7,nx,nd)
        call y5(7,ny,nd)
        call z2(7,nz,nd)
        call x1(11,nx,nd)
        call y4(11,ny,nd)
        call z1(11,nz,nd)
        call x2(12,nx,nd)
        call y5(12,ny,nd)
        call z1(12,nz,nd)
        call x1(2,nx,nd)
        call y4(2,ny,nd)
        call z2(2,nz,nd)
        call x2(3,nx,nd)
        call y5(3,ny,nd)
        call z2(3,nz,nd)
        call x2(5,nx,nd)
        call y4(5,ny,nd)
        call z1(5,nz,nd)
        call x2(18,nx,nd)
        call y4(18,ny,nd)
        call z2(18,nz,nd)

      else if(fs==1 .and. fw==0 .and. fn==0 .and. fe==1 .and. fb==1) then
        call x5(16,nx,nd)
        call y1(16,ny,nd)
        call z1(16,nz,nd)
        call x5(8,nx,nd)
        call y1(8,ny,nd)
        call z2(8,nz,nd)
        call x4(10,nx,nd)
        call y1(10,ny,nd)
        call z1(10,nz,nd)
        call x5(13,nx,nd)
        call y2(13,ny,nd)
        call z1(13,nz,nd)
        call x4(1,nx,nd)
        call y1(1,ny,nd)
        call z2(1,nz,nd)
        call x5(4,nx,nd)
        call y2(4,ny,nd)
        call z2(4,nz,nd)
        call x4(5,nx,nd)
        call y2(5,ny,nd)
        call z1(5,nz,nd)
        call x4(18,nx,nd)
        call y2(18,ny,nd)
        call z2(18,nz,nd)

      else if(fs==0 .and. fw==0 .and. fn==1 .and. fe==1 .and. fb==1) then
        call x5(17,nx,nd)
        call y5(17,ny,nd)
        call z1(17,nz,nd)
        call x5(9,nx,nd)
        call y5(9,ny,nd)
        call z2(9,nz,nd)
        call x4(12,nx,nd)
        call y5(12,ny,nd)
        call z1(12,nz,nd)
        call x5(13,nx,nd)
        call y4(13,ny,nd)
        call z1(13,nz,nd)
        call x4(3,nx,nd)
        call y5(3,ny,nd)
        call z2(3,nz,nd)
        call x5(4,nx,nd)
        call y4(4,ny,nd)
        call z2(4,nz,nd)
        call x4(5,nx,nd)
        call y4(5,ny,nd)
        call z1(5,nz,nd)
        call x4(18,nx,nd)
        call y4(18,ny,nd)
        call z2(18,nz,nd)

C     X NORMAL WHOLE SPACE SLICES

      else if(fs==1 .and. fw==1 .and. fn==1 .and. fe==0 .and. fb==1) then
        call x2(1,nx,nd)
        call y1(1,ny,nd)
        call z2(1,nz,nd)
        call x1(2,nx,nd)
        call y6(2,ny,nd)
        call z2(2,nz,nd)
        call x2(3,nx,nd)
        call y5(3,ny,nd)
        call z2(3,nz,nd)
        call x2(5,nx,nd)
        call y6(5,ny,nd)
        call z1(5,nz,nd)
        call x1(6,nx,nd)
        call y1(6,ny,nd)
        call z2(6,nz,nd)
        call x1(7,nx,nd)
        call y5(7,ny,nd)
        call z2(7,nz,nd)
        call x2(10,nx,nd)
        call y1(10,ny,nd)
        call z1(10,nz,nd)
        call x1(11,nx,nd)
        call y6(11,ny,nd)
        call z1(11,nz,nd)
        call x2(12,nx,nd)
        call y5(12,ny,nd)
        call z1(12,nz,nd)
        call x1(14,nx,nd)
        call y1(14,ny,nd)
        call z1(14,nz,nd)
        call x1(15,nx,nd)
        call y5(15,ny,nd)
        call z1(15,nz,nd)
        call x2(18,nx,nd)
        call y6(18,ny,nd)
        call z2(18,nz,nd)

      else if(fs==1 .and. fw==0 .and. fn==1 .and. fe==0 .and. fb==1) then
        call x3(1,nx,nd)
        call y1(1,ny,nd)
        call z2(1,nz,nd)
        call x3(3,nx,nd)
        call y5(3,ny,nd)
        call z2(3,nz,nd)
        call x3(5,nx,nd)
        call y6(5,ny,nd)
        call z1(5,nz,nd)
        call x3(10,nx,nd)
        call y1(10,ny,nd)
        call z1(10,nz,nd)
        call x3(12,nx,nd)
        call y5(12,ny,nd)
        call z1(12,nz,nd)
        call x3(18,nx,nd)
        call y6(18,ny,nd)
        call z2(18,nz,nd)

      else if(fs==1 .and. fw==0 .and. fn==1 .and. fe==1 .and. fb==1) then
        call x4(1,nx,nd)
        call y1(1,ny,nd)
        call z2(1,nz,nd)
        call x4(3,nx,nd)
        call y5(3,ny,nd)
        call z2(3,nz,nd)
        call x5(4,nx,nd)
        call y6(4,ny,nd)
        call z2(4,nz,nd)
        call x4(5,nx,nd)
        call y6(5,ny,nd)
        call z1(5,nz,nd)
        call x5(8,nx,nd)
        call y1(8,ny,nd)
        call z2(8,nz,nd)
        call x5(9,nx,nd)
        call y5(9,ny,nd)
        call z2(9,nz,nd)
        call x4(10,nx,nd)
        call y1(10,ny,nd)
        call z1(10,nz,nd)
        call x4(12,nx,nd)
        call y5(12,ny,nd)
        call z1(12,nz,nd)
        call x5(13,nx,nd)
        call y6(13,ny,nd)
        call z1(13,nz,nd)
        call x5(16,nx,nd)
        call y1(16,ny,nd)
        call z1(16,nz,nd)
        call x5(17,nx,nd)
        call y5(17,ny,nd)
        call z1(17,nz,nd)
        call x4(18,nx,nd)
        call y6(18,ny,nd)
        call z2(18,nz,nd)

C     corner cases when single process for each direction

      else if(fs==1 .and. fw==1 .and. fn==1 .and. fe==0 .and. fb==0) then
        call x2(1,nx,nd)
        call y1(1,ny,nd)
        call z3(1,nz,nd)
        call x1(2,nx,nd)
        call y6(2,ny,nd)
        call z3(2,nz,nd)
        call x2(3,nx,nd)
        call y5(3,ny,nd)
        call z3(3,nz,nd)
        call x1(6,nx,nd)
        call y1(6,ny,nd)
        call z3(6,nz,nd)
        call x1(7,nx,nd)
        call y5(7,ny,nd)
        call z3(7,nz,nd)
        call x2(18,nx,nd)
        call y6(18,ny,nd)
        call z3(18,nz,nd)

      else if(fs==1 .and. fw==0 .and. fn==1 .and. fe==0 .and. fb==0) then
        call x3(1,nx,nd)
        call y1(1,ny,nd)
        call z3(1,nz,nd)
        call x3(3,nx,nd)
        call y5(3,ny,nd)
        call z3(3,nz,nd)
        call x3(18,nx,nd)
        call y6(18,ny,nd)
        call z3(18,nz,nd)

      else if(fs==1 .and. fw==0 .and. fn==1 .and. fe==1 .and. fb==0) then
        call x4(1,nx,nd)
        call y1(1,ny,nd)
        call z3(1,nz,nd)
        call x4(3,nx,nd)
        call y5(3,ny,nd)
        call z3(3,nz,nd)
        call x5(4,nx,nd)
        call y6(4,ny,nd)
        call z3(4,nz,nd)
        call x5(8,nx,nd)
        call y1(8,ny,nd)
        call z3(8,nz,nd)
        call x5(9,nx,nd)
        call y5(9,ny,nd)
        call z3(9,nz,nd)
        call x4(18,nx,nd)
        call y6(18,ny,nd)
        call z3(18,nz,nd)

C     Y NORMAL WHOLE SPACE SLICES
      else if(fs==1 .and. fw==1 .and. fn==0 .and. fe==1 .and. fb==1) then
        call x6(1,nx,nd)
        call y1(1,ny,nd)
        call z2(1,nz,nd)
        call x1(2,nx,nd)
        call y2(2,ny,nd)
        call z2(2,nz,nd)
        call x6(5,nx,nd)
        call y2(5,ny,nd)
        call z1(5,nz,nd)
        call x5(4,nx,nd)
        call y2(4,ny,nd)
        call z2(4,nz,nd)
        call x1(6,nx,nd)
        call y1(6,ny,nd)
        call z2(6,nz,nd)
        call x5(8,nx,nd)
        call y1(8,ny,nd)
        call z2(8,nz,nd)
        call x6(10,nx,nd)
        call y1(10,ny,nd)
        call z1(10,nz,nd) 
        call x1(11,nx,nd)
        call y2(11,ny,nd)
        call z1(11,nz,nd) 
        call x5(13,nx,nd)
        call y2(13,ny,nd)
        call z1(13,nz,nd)
        call x1(14,nx,nd)
        call y1(14,ny,nd)
        call z1(14,nz,nd)
        call x5(16,nx,nd)
        call y1(16,ny,nd)
        call z1(16,nz,nd)
        call x6(18,nx,nd)
        call y2(18,ny,nd)
        call z2(18,nz,nd)

      else if(fs==0 .and. fw==1 .and. fn==1 .and. fe==1 .and. fb==1) then
        call x1(2,nx,nd)
        call y4(2,ny,nd)
        call z2(2,nz,nd)
        call x6(3,nx,nd)
        call y5(3,ny,nd)
        call z2(3,nz,nd)
        call x5(4,nx,nd)
        call y4(4,ny,nd)
        call z2(4,nz,nd)
        call x6(5,nx,nd)
        call y4(5,ny,nd)
        call z1(5,nz,nd)
        call x1(7,nx,nd)
        call y5(7,ny,nd)
        call z2(7,nz,nd)
        call x5(9,nx,nd)
        call y5(9,ny,nd)
        call z2(9,nz,nd)
        call x1(11,nx,nd)
        call y4(11,ny,nd)
        call z1(11,nz,nd)
        call x6(12,nx,nd)
        call y5(12,ny,nd)
        call z1(12,nz,nd)
        call x5(13,nx,nd)
        call y4(13,ny,nd)
        call z1(13,nz,nd)
        call x1(15,nx,nd)
        call y5(15,ny,nd)
        call z1(15,nz,nd)
        call x5(17,nx,nd)
        call y5(17,ny,nd)
        call z1(17,nz,nd)
        call x6(18,nx,nd)
        call y4(18,ny,nd)
        call z2(18,nz,nd)

      else if(fs==0 .and. fw==1 .and. fn==0 .and. fe==1 .and. fb==1) then
        call x1(2,nx,nd)
        call y3(2,ny,nd)
        call z2(2,nz,nd)
        call x5(4,nx,nd)
        call y3(4,ny,nd)
        call z2(4,nz,nd)
        call x6(5,nx,nd)
        call y3(5,ny,nd)
        call z1(5,nz,nd)
        call x1(11,nx,nd)
        call y3(11,ny,nd)
        call z1(11,nz,nd)
        call x5(13,nx,nd)
        call y3(13,ny,nd)
        call z1(13,nz,nd)
        call x6(18,nx,nd) 
        call y3(18,ny,nd)
        call z2(18,nz,nd)

      else if(fs==1 .and. fw==1 .and. fn==0 .and. fe==1 .and. fb==0) then
        call x6(1,nx,nd)
        call y1(1,ny,nd)
        call z3(1,nz,nd)
        call x1(2,nx,nd)
        call y2(2,ny,nd)
        call z3(2,nz,nd)
        call x5(4,nx,nd)
        call y2(4,ny,nd)
        call z3(4,nz,nd)
        call x1(6,nx,nd)
        call y1(6,ny,nd)
        call z3(6,nz,nd)
        call x5(8,nx,nd)
        call y1(8,ny,nd)
        call z3(8,nz,nd)
        call x6(18,nx,nd)
        call y2(18,ny,nd)
        call z3(18,nz,nd)

      else if(fs==0 .and. fw==1 .and. fn==1 .and. fe==1 .and. fb==0) then
        call x1(2,nx,nd)
        call y4(2,ny,nd)
        call z3(2,nz,nd)
        call x6(3,nx,nd)
        call y5(3,ny,nd)
        call z3(3,nz,nd)
        call x5(4,nx,nd)
        call y4(4,ny,nd)
        call z3(4,nz,nd)
        call x1(7,nx,nd)
        call y5(7,ny,nd)
        call z3(7,nz,nd)
        call x5(9,nx,nd)
        call y5(9,ny,nd)
        call z3(9,nz,nd)
        call x6(18,nx,nd)
        call y4(18,ny,nd)
        call z3(18,nz,nd)

      else if(fs==0 .and. fw==1 .and. fn==0 .and. fe==1 .and. fb==0) then
        call x1(2,nx,nd)
        call y3(2,ny,nd)
        call z3(2,nz,nd)
        call x5(4,nx,nd)
        call y3(4,ny,nd)
        call z3(4,nz,nd)
        call x6(18,nx,nd) 
        call y3(18,ny,nd)
        call z3(18,nz,nd)

      else if(fs==1 .and. fw==1 .and. fn==1 .and. fe==1 .and. fb==0) then
        call x6(1,nx,nd)
        call y1(1,ny,nd)
        call z3(1,nz,nd)
        call x1(2,nx,nd)
        call y6(2,ny,nd)
        call z3(2,nz,nd)
        call x6(3,nx,nd)
        call y5(3,ny,nd)
        call z3(3,nz,nd)
        call x5(4,nx,nd)
        call y6(4,ny,nd)
        call z3(4,nz,nd)
        call x6(18,nx,nd)
        call y6(18,ny,nd)
        call z3(18,nz,nd)


C     INTERNAL PROCS
 
      else if(fs==0 .and. fw==0 .and. fn==0 .and. fe==0 .and. fb==0) then
        call x3(18,nx,nd)
        call y3(18,ny,nd)
        call z3(18,nz,nd)

C     WHOLE MODEL SPACE

      else if(fs==1 .and. fw==1 .and. fn==1 .and. fe==1 .and. fb==1) then
        call x6(1,nx,nd)
        call y1(1,ny,nd)
        call z2(1,nz,nd)
        call x1(2,nx,nd)
        call y6(2,ny,nd)
        call z2(2,nz,nd)
        call x6(3,nx,nd)
        call y5(3,ny,nd)
        call z2(3,nz,nd)
        call x5(4,nx,nd)
        call y6(4,ny,nd)
        call z2(4,nz,nd)
        call x6(5,nx,nd)
        call y6(5,ny,nd)
        call z1(5,nz,nd)
        call x1(6,nx,nd)
        call y1(6,ny,nd)
        call z2(6,nz,nd)
        call x1(7,nx,nd)
        call y5(7,ny,nd)
        call z2(7,nz,nd)
        call x5(8,nx,nd)
        call y1(8,ny,nd)
        call z2(8,nz,nd)
        call x5(9,nx,nd)
        call y5(9,ny,nd)
        call z2(9,nz,nd)
        call x6(10,nx,nd)
        call y1(10,ny,nd)
        call z1(10,nz,nd)
        call x1(11,nx,nd)
        call y6(11,ny,nd)
        call z1(11,nz,nd)
        call x6(12,nx,nd)
        call y5(12,ny,nd)
        call z1(12,nz,nd)
        call x5(13,nx,nd)
        call y6(13,ny,nd)
        call z1(13,nz,nd)
        call x1(14,nx,nd)
        call y1(14,ny,nd)
        call z1(14,nz,nd)
        call x1(15,nx,nd)
        call y5(15,ny,nd)
        call z1(15,nz,nd)
        call x5(16,nx,nd)
        call y1(16,ny,nd)
        call z1(16,nz,nd)
        call x5(17,nx,nd)
        call y5(17,ny,nd)
        call z1(17,nz,nd)
        call x6(18,nx,nd)
        call y6(18,ny,nd)
        call z2(18,nz,nd)
      end if

      return
      end



