#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

#define         RPERD           0.017453293
#define         DPERR           57.29577951

void srf2stoch(char* rup_geom_file, int slip_id, int hypo_id, float dx, float dy, int inbin, float avgstk, struct slipfile* sfile)
{
FILE *fopfile(), *fpr;
float len, wid, xp, yp;
float dtop, strike, dip;
float rr, ss, tt, tl, ravg, savg, *slip, *trise, *tinit, *sp, *tr, *ti;
float s1, s2, xx, yy, rake, cosA, sinA, sum;
int i, is, id, ip, nstk, ndip, nseg;
int ix, iy, j, k, m, kp, noff;
int nx, ny, nxdiv, nydiv, nxsum, nysum;
float shypo, dhypo, elon, elat, fac;
char string[1024];

struct standrupformat srf;

int revstk;

if(avgstk > -1.0e+14)
   {
   while(avgstk >= 360.0)
      avgstk = avgstk - 360.0;
   while(avgstk < 0.0)
      avgstk = avgstk + 360.0;
   }

//read_srf(&srf,infile,inbin);
rg_stats_t stats;
rupgen_genslip(rup_geom_file, slip_id, hypo_id, &stats, &srf);
char outfile[256];
sprintf(outfile, "%d_%d.srf", slip_id, hypo_id);
write_srf(&srf, outfile, 0);

//allocate
sfile->nseg = srf.srf_prect.nseg;

noff = 0;
for(i=0;i<srf.srf_prect.nseg;i++)
   {
   elon = srf.srf_prect.prectseg[i].elon;
   printf("%f\n", elon);
   elat = srf.srf_prect.prectseg[i].elat;
   nstk = srf.srf_prect.prectseg[i].nstk;
   ndip = srf.srf_prect.prectseg[i].ndip;
   len = srf.srf_prect.prectseg[i].flen;
   wid = srf.srf_prect.prectseg[i].fwid;

   strike = srf.srf_prect.prectseg[i].stk;
   dip = srf.srf_prect.prectseg[i].dip;
   dtop = srf.srf_prect.prectseg[i].dtop;
   shypo = srf.srf_prect.prectseg[i].shyp;
   dhypo = srf.srf_prect.prectseg[i].dhyp;

   while(strike >= 360.0)
      strike = strike - 360.0;
   while(strike < 0.0)
      strike = strike + 360.0;

   nx = (len/dx + 0.5);
   if(nx > nstk)
      {
      nx = nstk;
      dx = len/nx;
      }

   nxdiv = 1;
   while((nstk*nxdiv)%nx)
      nxdiv++;

   nxsum = (nstk*nxdiv)/nx;

   ny = (wid/dy + 0.5);
   if(ny > ndip)
      {
      ny = ndip;
      dy = wid/ny;
      }

   nydiv = 1;
   while((ndip*nydiv)%ny)
      nydiv++;

   nysum = (ndip*nydiv)/ny;

   fprintf(stderr,"seg= %d\n",i);
   fprintf(stderr,"nstk= %d nx= %d nxdiv= %d nxsum= %d\n",nstk,nx,nxdiv,nxsum);
   fprintf(stderr,"ndip= %d ny= %d nydiv= %d nysum= %d\n",ndip,ny,nydiv,nysum);

   slip = (float *) check_malloc (nxdiv*nydiv*nstk*ndip*sizeof(float));
   trise = (float *) check_malloc (nxdiv*nydiv*nstk*ndip*sizeof(float));
   tinit = (float *) check_malloc (nxdiv*nydiv*nstk*ndip*sizeof(float));
   sp = (float *) check_malloc (nx*ny*sizeof(float));
   tr = (float *) check_malloc (nx*ny*sizeof(float));
   ti = (float *) check_malloc (nx*ny*sizeof(float));

   savg = 0.0;
   ravg = 0.0;
   for(id=0;id<ndip;id++)
      {
      for(is=0;is<nstk;is++)
         {
	 kp = noff + is + id*nstk;

	 tt = srf.srf_apnts.apntvals[kp].tinit;
	 s1 = srf.srf_apnts.apntvals[kp].slip1;
	 s2 = srf.srf_apnts.apntvals[kp].slip2;

	 rake = srf.srf_apnts.apntvals[kp].rake;

	 cosA = cos(rake*RPERD);
	 sinA = sin(rake*RPERD);

	 xx = -s2*sinA + s1*cosA;
	 yy =  s2*cosA + s1*sinA;

	 ss = sqrt(xx*xx + yy*yy);

	 rr = 90;
	 if(yy < 0.0)
	    rr = 270;

	 if(xx != 0.0)
	    {
	    rr = DPERR*atan(yy/xx);
	    if(xx < 0.0)
	       rr = rr + 180;
	    }

	 while(rr < 0.0)
	    rr = rr + 360.0;
	 while(rr > 360.0)
	    rr = rr - 360.0;

	 tl = (srf.srf_apnts.apntvals[kp].dt)*(srf.srf_apnts.apntvals[kp].nt1 - 1);
	 if(srf.srf_apnts.apntvals[kp].nt2 > srf.srf_apnts.apntvals[kp].nt1)
	    tl = (srf.srf_apnts.apntvals[kp].dt)*(srf.srf_apnts.apntvals[kp].nt2 - 1);
	 if(tl < 0.0)
	    tl = 0.0;
	 for(k=0;k<nydiv;k++)
	    {
	    for(j=0;j<nxdiv;j++)
	       {
	       ip = is*nxdiv + j + (id*nydiv + k)*nstk*nxdiv;

	       slip[ip] = ss;
	       trise[ip] = tl;
	       tinit[ip] = tt;
	       }
	    }

         ravg = ravg + rr*ss;
         savg = savg + ss;
         }
      }

   ravg = ravg/savg;
   savg = savg/(nstk*ndip);

   fac = 1.0/(float)(nxsum*nysum);

   for(iy=0;iy<ny;iy++)
      {
      for(ix=0;ix<nx;ix++)
         {
	 ip = ix + iy*nx;

         sp[ip] = 0.0;
         tr[ip] = 0.0;
         ti[ip] = 0.0;
	 sum = 0.0;
	 for(k=0;k<nysum;k++)
	    {
	    for(j=0;j<nxsum;j++)
	       {
	       m = ix*nxsum + j + (iy*nysum + k)*nx*nxsum;

	       sp[ip] = sp[ip] + slip[m];
	       ti[ip] = ti[ip] + tinit[m];
	       tr[ip] = tr[ip] + trise[m]*slip[m];
	       sum = sum + slip[m];
	       }
	    }

         sp[ip] = sp[ip]*fac;
         ti[ip] = ti[ip]*fac;

	 if(sum > 0.0)
            tr[ip] = tr[ip]/sum;
	 else
            tr[ip] = 1.0e-05;
         }
      }

   sfile->elon[i] = elon;
   sfile->elat[i] = elat;
   sfile->nx[i] = nx;
   sfile->ny[i] = ny;
   sfile->dx[i] = dx;
   sfile->dy[i] = dy;
   sfile->strike[i] = strike;
   sfile->dip[i] = dip;
   sfile->ravg[i] = ravg;
   sfile->dtop[i] = dtop;
   sfile->shypo[i] = shypo;
   sfile->dhypo[i] = dhypo;
//   fprintf(fpw,"%10.4f %10.4f %5d %5d %8.2f %8.2f\n",elon,elat,nx,ny,dx,dy);
//   fprintf(fpw,"%4.0f %4.0f %4.0f %8.2f %8.2f %8.2f\n",strike,dip,ravg,dtop,shypo,dhypo);

   revstk = 0;
   if(avgstk > -1.0e+14 && (strike > avgstk + 90.0 || strike < avgstk - 90.0))
      revstk = 1;

   printf("revstk: %d\n", revstk);

//In for loop, need to reverse indices for fortran compatibility
//so it's y, x, segment
   for(iy=0;iy<ny;iy++)
      {
      if(revstk)
         {
         for(ix=nx-1;ix>=0;ix--)
	      sfile->sp[iy*NQ*LV + (nx-1-ix)*LV + i] = sp[ix + iy*nx];
//	      sfile->sp[iy][nx-1-ix][i] = sp[ix + iy*nx];
//	      sfile->sp[i][nx-1-ix + iy*nx] = sp[ix + iy*nx];
//            fprintf(fpw," %5.0f",sp[ix + iy*nx]);
	 }
      else
        {
         for(ix=0;ix<nx;ix++)
	      sfile->sp[iy*NQ*LV + ix*LV + i] = sp[ix + iy*nx];
//	      sfile->sp[iy][ix][i] = sp[ix + iy*nx];
//            fprintf(fpw," %5.0f",sp[ix + iy*nx]);
	 }

//      fprintf(fpw,"\n");
      }

   for(iy=0;iy<ny;iy++)
      {
      if(revstk)
         {
         for(ix=nx-1;ix>=0;ix--)
	      sfile->tr[iy*NQ*LV + (nx-1-ix)*LV + i] = tr[ix + iy*nx];
//	      sfile->tr[iy][nx-1-ix][i] = tr[ix + iy*nx];
//            sfile->tr[i][nx-1-ix + iy*nx] = tr[ix + iy*nx];
//            fprintf(fpw," %5.2f",tr[ix + iy*nx]);
	 }
      else
         {
         for(ix=0;ix<nx;ix++)
	      sfile->tr[iy*NQ*LV + ix*LV + i] = tr[ix + iy*nx];
//	      sfile->tr[iy][ix][i] = tr[ix + iy*nx];
//            fprintf(fpw," %5.2f",tr[ix + iy*nx]);
	 }

//      fprintf(fpw,"\n");
      }

   for(iy=0;iy<ny;iy++)
      {
      if(revstk)
         {
         for(ix=nx-1;ix>=0;ix--)
	      sfile->ti[iy*NQ*LV + (nx-1-ix)*LV + i] = ti[ix + iy*nx];
//	      sfile->ti[iy][nx-1-ix][i] = ti[ix + iy*nx];
//	      sfile->ti[i][nx-1-ix + iy*nx] = ti[ix + iy*nx];
 //           fprintf(fpw," %5.2f",ti[ix + iy*nx]);
	 }
      else
         {
         for(ix=0;ix<nx;ix++)
	      sfile->ti[iy*NQ*LV + ix*LV + i] = ti[ix + iy*nx];
//	      sfile->ti[iy][ix][i] = ti[ix + iy*nx];
//            fprintf(fpw," %5.2f",ti[ix + iy*nx]);
	 }

//      fprintf(fpw,"\n");
      }

   free(slip);
   free(trise);
   free(tinit);
   free(sp);
   free(tr);
   free(ti);

   noff = noff + nstk*ndip;
   }
}
