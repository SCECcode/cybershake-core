CCC   'relax.f' CALCULATES RELAXATION ARRAY AND MODULI


      subroutine tausub(taumin,taumax)

CCC   CALCULATES THE RELAXATION TIME ARRAY

      use parstat

      integer, parameter :: ntau=2

      real, dimension(ntau,ntau,ntau) :: tautem
      real taumin,taumax

         tautem(1,1,1)=1.
         tautem(2,1,1)=6.
         tautem(1,2,1)=7.
         tautem(2,2,1)=4.
         tautem(1,1,2)=8.
         tautem(2,1,2)=3.
         tautem(1,2,2)=2.
         tautem(2,2,2)=5.

      do 55 idx = 1,ntau
        do 55 idy = 1,ntau
          do 55 idz = 1,ntau

            tmp=tautem(idx,idy,idz)

             tmp = (tmp-0.5)/8
             tmp = 2.0*tmp - 1.0

             tau(idx,idy,idz) = exp(0.5*(log(taumax*taumin) +
     +                          log(taumax/taumin)*tmp))

55     continue

       return
       end



      subroutine relaxmod(fp,fl,fh,vp,vs,qqp,qqs,tmin,tmax)

CCC   CALCULATE THE EFFECTIVE RELAXATION MODULI (BRADLEY 7/98)
CCC   REPLACES QS, QP ARRAYS WITH THE CALCULATED MODULI

        real tmax,tmin,w0,w1,w2,qpinv,qsinv

        pi=4.*atan(1.)

C radial frequency for the center, upper and lower bounds of the Q range

        w0=2*pi*fp
        w1=2*pi*fl
        w2=2*pi*fh

        tmax=1./w1
        tmin=1./w2

        if (qqp .le. 0.0) then
          qpinv=0.0
          qsinv=0.0
        else
          qpinv=1./qqp
          qsinv=1./qqs
        end if

        qqp = 2./pi*qpinv*(log(tmax/tmin))/
     +       (1.0-2./pi*qpinv*log(w0*tmin))

        qqs = 2./pi*qsinv*(log(tmax/tmin))/
     +       (1.0-2./pi*qsinv*log(w0*tmin))

        return
        end


