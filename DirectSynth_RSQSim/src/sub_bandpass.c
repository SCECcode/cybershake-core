/* subroutine to perform Butterworth filter.
Original code is from Stanford SEP group (David Hale).

For this subroutine, make sure the incoming and outgoing
arrays are allocated with two extra samples:
	y[nt+2] and x[nt+2].

For most other applications the exra two samples will never be
used.  Unfortunately this code needs them.

stanford_bandpass_(y_in, x_out, nt,dt, flo, fhi, npolelo, npolehi, phase)
*/

/* bandpass Butterworth filter 

bandpass nt= [ dt=1. nx=all ny=1 in=stdin out=stdout flo=0. fhi=0.5/dt 
	     nplo=6 nphi=6 phase=0 ]

nt, dt, nx, ny	trace length, sampling interval, number of traces
in, out		input, output filenames
flo, fhi	low and high cutoff frequencies; 
		flo=0 for lowpass filter; fhi=0.5/dt for highpass;
		gain at flo and fhi is:	-6db for zero-phase filter
					-3db for minimum-phase filter
nplo, nphi	number of poles for low and high cutoffs;
		filter rolloff rate is proportional to number of poles
phase		=0 for zero-phase; =1 for minimum-phase

Dave Hale, 2/25/82
modified to run in core, by Jan Morton, 1/17/84
Technical reference: Oppenheim, A.V., and Schafer, R.W., 1975, 
Digital signal processing, Prentice Hall, Inc.
*/
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h> 
void stanford_bandpass_(y,x, nt_in,dt_in, flow,fhigh, npolelo,npolehi,phase_in)
float	y[],x[],*dt_in,*flow,*fhigh;
int	*nt_in,*npolelo,*npolehi,*phase_in;
  {
	int ok,nt2,j,n,nb,lopass,hipass,nd, nt,nplo,nphi,phase;
	int m,i,nodd;
	float b[50][5],d[50][5];
	float dt,flo,fhi,c,a,aa,aap4,dtheta,theta0,pi=3.141593;
	float e,ee,b1,b2,b3,den,fno2;
	double sin(),cos();

	nt = *nt_in;		 dt = *dt_in;
	nplo = *npolelo;	 nphi = *npolehi;
	phase = *phase_in;


	if (nplo < 1) nplo = 1; if (nplo > 98) nplo = 98;
	if (nphi < 1) nphi = 1; if (nphi > 98) nphi = 98;

	flo = *flow * dt;
	fhi = *fhigh * dt;
	hipass = (flo > 0.0001);
	lopass = (fhi < 0.4999);

	nt2 = nt+2;

	/* compute lowpass filter coefficients if required */
	if (lopass)
	  {
		nodd = (nphi%2);
		if (phase==0)
		  {
			if (nodd) nphi = (nphi+1)/2;
			else nphi = nphi/2;
			nodd = (nphi%2);
		  }
		nb = (nphi+1)/2;

		a = 2.*sin(pi*fhi)/cos(pi*fhi); /* radius of poles in s-plane */
		aa = a*a;
		aap4 = aa+4;
		dtheta = pi/nphi;		/* angular separation of poles */
		theta0 = (nodd)?0:dtheta/2;	/* pole closest to real s axis */
		if (nodd)
		  {
			b[0][0] = a/(a+2); 
			b[0][1] = b[0][0]; 
			b[0][2] = 0.;
			b[0][3] = (a-2)/(a+2);	
			b[0][4] = 0.;
		  }
		else
		  {
			c = 4.*a*cos(theta0);
			b[0][0] = aa/(aap4+c);
			b[0][1] = 2.*b[0][0];
			b[0][2] = b[0][0];
			b[0][3] = (2.*aa-8.)/(aap4+c);
			b[0][4] = (aap4-c)/(aap4+c);
		  }
		for (j = 1; j < nb; j++)
		  {
			c = 4.*a*cos(theta0+j*dtheta);
			b[j][0] = aa/(aap4+c);
			b[j][1] = 2.*b[j][0];
			b[j][2] = b[j][0];
			b[j][3] = (2.*aa-8.)/(aap4+c);
			b[j][4] = (aap4-c)/(aap4+c);
		  }
	  }

	/* compute highpass filter coefficients if required by
	   transforming a lowpass filter with cutoff at half Nyquist */
	if (hipass)
	  {
		nodd = (nplo%2);
		if (phase==0)
		  {
			if (nodd) nplo = (nplo+1)/2;
			else nplo = nplo/2;
			nodd = (nplo%2);
		  }
		nd = (nplo+1)/2;
		fno2 = 0.25;
		a = 2.*sin(pi*fno2)/cos(pi*fno2); aa = a*a; aap4 = aa+4;
		e = -cos(pi*(flo+fno2))/cos(pi*(flo-fno2)); ee = e*e;
		dtheta = pi/nplo; theta0 = (nodd)?0:dtheta/2;
		if (nodd)
		  {
			b1 = a/(a+2);
			b2 = (a-2)/(a+2);
			den = 1.-b2*e;
			d[0][0] = b1*(1.-e)/den;
			d[0][1] = -d[0][0];
			d[0][2] = 0.;
			d[0][3] = (e-b2)/den;
			d[0][4] = 0.;
		  }
		else
		  {
			c = 4.*a*cos(theta0);
			b1 = aa/(aap4+c);
			b2 = (2.*aa-8.)/(aap4+c);
			b3 = (aap4-c)/(aap4+c);
			den = 1.-b2*e+b3*ee;
			d[0][0] = b1*(1.-e)*(1.-e)/den;
			d[0][1] = -2.*d[0][0];
			d[0][2] = d[0][0];
			d[0][3] = (2.*e*(1.+b3)-b2*(1.+ee))/den;
			d[0][4] = (ee-b2*e+b3)/den;
		  }
		for (j = 1; j < nd; j++)
		  {
			c = 4.*a*cos(theta0+j*dtheta);
			b1 = aa/(aap4+c);
			b2 = (2.*aa-8.)/(aap4+c);
			b3 = (aap4-c)/(aap4+c);
			den = 1.-b2*e+b3*ee;
			d[j][0] = b1*(1.-e)*(1.-e)/den;
			d[j][1] = -2.*d[j][0];
			d[j][2] = d[j][0];
			d[j][3] = (2.*e*(1.+b3)-b2*(1.+ee))/den;
			d[j][4] = (ee-b2*e+b3)/den;
		  }
	  }


/* Now do seismogram: */

	/* apply filters */
		for (m=2;m<nt2;m++)
			x[m]=y[m-2];
		x[0]=0.;
		y[0]=0.;
		x[1]=0.;
		y[1]=0.;
		if (lopass)
		{
			for (i=nb; i>0; i--)
			{
				j=nb-i;
				for (m=2; m<nt2; m++)
				y[m]=b[j][0]*x[m]+b[j][1]*x[m-1]+b[j][2]*x[m-2]-b[j][3]*y[m-1]-b[j][4]*y[m-2];
				for (m=2;m<nt2;m++)
					x[m]=y[m];
			}
		if (phase==0)
			{
			for (m=2;m<nt2;m++)
				y[m]=x[nt2+1-m];
			for (i=nb; i>0; i--)
			{
				j=nb-i;
				for (m=2; m<nt2; m++)
				x[m]=b[j][0]*y[m]+b[j][1]*y[m-1]+b[j][2]*y[m-2]-b[j][3]*x[m-1]-b[j][4]*x[m-2];
				for (m=2; m<nt2; m++)
					y[m]=x[m];
			}
			for (m=2; m<nt2; m++)
				x[m]=y[nt2+1-m];
			}
		}
		if (hipass)
		{
			for (i=nd; i>0; i--)
			{
				j=nd-i;
				for (m=2; m<nt2; m++)
				y[m]=d[j][0]*x[m]+d[j][1]*x[m-1]+d[j][2]*x[m-2]-d[j][3]*y[m-1]-d[j][4]*y[m-2];
				for (m=2;m<nt2;m++)
					x[m]=y[m];
			}
		if (phase==0)
			{
			for (m=2;m<nt2;m++)
				y[m]=x[nt2+1-m];
			for (i=nd; i>0; i--)
			{
				j=nd-i;
				for (m=2; m<nt2; m++)
				x[m]=d[j][0]*y[m]+d[j][1]*y[m-1]+d[j][2]*y[m-2]-d[j][3]*x[m-1]-d[j][4]*x[m-2];
				for (m=2; m<nt2; m++)
					y[m]=x[m];
			}
			for (m=2; m<nt2; m++)
				x[m]=y[nt2+1-m];
			}
		}

  }


