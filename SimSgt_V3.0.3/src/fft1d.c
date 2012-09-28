
float sin_table[] =
   {
	1.0000000e+00,	/* sin(pi/2) */
	7.0710678e-01,	/* sin(pi/4) */
	3.8268343e-01,	/* sin(pi/8) */
	1.9509032e-01,	/* sin(pi/16) */
	9.8017140e-02,	/* sin(pi/32) */
	4.9067674e-02,	/* sin(pi/64) */
	2.4541228e-02,	/* sin(pi/128) */
	1.2271538e-02,	/* sin(pi/256) */
	6.1358846e-03,	/* sin(pi/512) */
	3.0679568e-03,	/* sin(pi/1024) */
	1.5339802e-03,	/* sin(pi/2048) */
	7.6699032e-04,	/* sin(pi/4096) */
	3.8349519e-04,	/* sin(pi/8192) */
	1.9174760e-04,	/* sin(pi/16384) */
	9.5873799e-05, /* sin(pi/32768) */
	4.7936899e-05, /* sin(pi/65536) */
	2.3968499e-05, /* sin(pi/131072) */
	1.1984224e-05, /* sin(pi/262144) */
	5.9921125e-06, /* sin(pi/524288) */
	2.9960562e-06, /* sin(pi/1048576) */
	1.4980281e-06  /* sin(pi/2097152) */
   };
/*
 * cfft - radix 2 FFT for complex data
 *
 * n is the number of complex points.
 * x is the single precision (not double!!) complex vector.
 * isign is the sign of the transform.
 *
 * the routine does no normalization
 */
struct complex { float re; float im; };
cfft(x,n,isign)
struct complex *x;
int n,isign;
   {
	register struct complex *px, *qx, *rx;
	struct complex *limit, *qlimit, dtemp;
	float cn, sn, cd, sd, temp, real, imag;
	int m, j, istep;
	float *psintab;
	extern float sin_table[];

	limit= x + n;
	j= 0;
	for(px=x; px<limit; px++)
	   {
		if(px < (qx= x+j))
		   {	dtemp= *qx; *qx= *px; *px= dtemp;	}
		m = n>>1;
		while( m>=1 && j>=m )
		   { j-= m; m>>= 1;    }
		j+= m;
	   }
	rx= x+1;
	for(px=x; px<limit; px+= 2, rx+= 2)
	   {
		temp= rx->re;
		rx->re= px->re -temp;
		px->re += temp;
		temp= rx->im;
		rx->im= px->im -temp;
		px->im += temp;
	   }
	j=2;
	psintab= sin_table;
	while( j < n )
	   {
		istep= j<<1;
		sd= *psintab++;
		temp= *psintab;
		cd= 2.0 * temp * temp;
		cn= 1.0;
		sn= 0.0;
		if( isign < 0 ) sd= -sd;
		qlimit= x+j;
		for(qx=x; qx< qlimit; qx++)
		   {
			for(px=qx; px<limit; px+= istep)
			   {
				rx= px + j;
				real= cn * rx->re - sn * rx->im;
				imag= sn * rx->re + cn * rx->im;
				rx->re = px->re - real;
				rx->im = px->im - imag;
				px->re += real;
				px->im += imag;
			   }
			temp= cd * cn + sd * sn;
			sn += (sd * cn - cd * sn);
			cn -= temp;
		   }
		j= istep;
	   }
	return;
   }

forfft(x,n,isign)
register struct complex *x;
int n, isign;
   {
	register struct complex *px, *rx;
	float cn, sn, cd, sd, arg;
	float are, aim, bre, bim, real, imag;
	float *psintab;
	float half = 0.5;
	extern float sin_table[];
	int k;

	cfft_r(x,n/2,isign);

	/* do DC and Nyquist */
	real= x[0].re;
	imag= x[0].im;
	x[0].re= real + imag;
	x[0].im= real - imag;

	/* set up for sine recurrsion */
	psintab= sin_table;
	for(k=4; k<n; k <<= 1) psintab++;
	sd= *psintab++;
	real= *psintab;
	cd= 2.0 * real * real;

	sn= 0.0;
	cn= 1.0;
	if(isign < 0) sd= -sd;
	px= x + 1;
	rx= x + n/2 -1;
	while( px <= rx )
	   {
		real = cd*cn + sd*sn;
		imag = sd*cn - cd*sn;
		cn -= real;
		sn += imag;

		are= half*(px->re + rx->re);
		aim= half*(px->im - rx->im);
		bre= half*(px->im + rx->im);
		bim= half*(rx->re - px->re);

		real= bre*cn - bim*sn;
		imag= bre*sn + bim*cn;

		px->re = are + real;
		px->im = aim + imag;
		rx->re = are - real;
		rx->im = imag - aim;

		px++;
		rx--;
	   }
	if(abs(isign) > 1)
	   {
		x[n/2].re= x[0].im;
		x[0].im= x[n/2].im = 0.0;
	   }
   }

invfft(x,n,isign)
register struct complex *x;
int n, isign;
   {
	register struct complex *px, *rx;
	float cn, sn, cd, sd;
	float are, aim, bre, bim, real, imag;
	float *psintab;
	extern float sin_table[];
	int k;

	if(abs(isign) > 1) x[0].im= x[n/2].re;

	/* do DC and Nyquist */
	real= x[0].re;
	imag= x[0].im;
	x[0].re= real + imag;
	x[0].im= real - imag;

	/* set up for sine recurrsion */
	psintab= sin_table;
	for(k=4; k<n; k <<= 1) psintab++;
	sd= *psintab++;
	real= *psintab;
	cd= 2.0 * real * real;

	sn= 0.0;
	cn= 1.0;
	if(isign < 0) sd= -sd;
	px= x + 1;
	rx= x + n/2 -1;
	while( px <= rx )
	   {
		real = cd*cn + sd*sn;
		imag = sd*cn - cd*sn;
		cn -= real;
		sn += imag;

		are= (px->re + rx->re);
		aim= (px->im - rx->im);
		bre= (px->re - rx->re);
		bim= (px->im + rx->im);

		real= bre*cn - bim*sn;
		imag= bre*sn + bim*cn;

		px->re = are - imag;
		px->im = aim + real;
		rx->re = are + imag;
		rx->im = real - aim;

		px++;
		rx--;
	   }
	cfft_r(x,n/2,isign);
   }

cfft_r(x,n,isign)
struct complex *x;
int n,isign;
   {
	register struct complex *px, *qx, *rx;
	struct complex *limit, *qlimit, dtemp;
	float cn, sn, cd, sd, temp, real, imag;
	int m, j, istep;
	float *psintab;
	extern float sin_table[];

	limit= x + n;
	j= 0;
	for(px=x; px<limit; px++)
	   {
		if(px < (qx= x+j))
		   {	dtemp= *qx; *qx= *px; *px= dtemp;	}
		m = n>>1;
		while( m>=1 && j>=m )
		   { j-= m; m>>= 1;    }
		j+= m;
	   }
	rx= x+1;
	for(px=x; px<limit; px+= 2, rx+= 2)
	   {
		temp= rx->re;
		rx->re= px->re -temp;
		px->re += temp;
		temp= rx->im;
		rx->im= px->im -temp;
		px->im += temp;
	   }
	j=2;
	psintab= sin_table;
	while( j < n )
	   {
		istep= j<<1;
		sd= *psintab++;
		temp= *psintab;
		cd= 2.0 * temp * temp;
		cn= 1.0;
		sn= 0.0;
		if( isign < 0 ) sd= -sd;
		qlimit= x+j;
		for(qx=x; qx< qlimit; qx++)
		   {
			for(px=qx; px<limit; px+= istep)
			   {
				rx= px + j;
				real= cn * rx->re - sn * rx->im;
				imag= sn * rx->re + cn * rx->im;
				rx->re = px->re - real;
				rx->im = px->im - imag;
				px->re += real;
				px->im += imag;
			   }
			temp= cd * cn + sd * sn;
			sn += (sd * cn - cd * sn);
			cn -= temp;
		   }
		j= istep;
	   }
	return;
   }
