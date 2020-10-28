#include "include.h"

int main(int argc, char** argv) {
	int nt = 60000;
	float dt = 0.005;
	struct pntsrcs srcs;
	int nbounce = 1;
	char stype[256] = "cos";
	char name[256] = "source";
	char stfdir[256] = "Stf";
	char stf_file[256] = "fx_src_v3.0.8";
	int bfilt=4;
	float flo=2.0;
	float fhi=0.0;
	float tdelay = 0.0;
	float tzero = 0.1/(flo/2.0);
	struct runparamsP3 rpars_p3;
	struct nodeinfo ninfo;
	rpars_p3.geoproj = 1;
	rpars_p3.xshift = -1.0e+15;
	rpars_p3.yshift = -1.0e+15;

	rpars_p3.modelrot = 0.0;  /* rotation of y-axis from south (clockwise positive) */
	rpars_p3.modellat = 34.0;    /* latitude of model origin */
	rpars_p3.modellon = -118.0;    /* longitude of model origin */

	rpars_p3.xmom = 0.0;
	rpars_p3.ymom = 0.0;
	rpars_p3.zmom = 0.0;

	srcs.nsource = 1;
	srcs.eqsrc = 0;
	srcs.expl = 0;
	srcs.psrc = 0;
	srcs.bforce = 1;
	srcs.dblcpl = 0;
	srcs.pointmt = 0;
	srcs.ffault = 0;
	srcs.adjoint = 0;
	srcs.relative_slip = 0;
	srcs.absolute_slip = 0;
	srcs.area = 1.0;
	float it = 3.0/(dt*flo);
	tdelay = it*dt;
	rpars_p3.dt = dt;
	rpars_p3.tdelay = tdelay;
	setpar(argc, argv);
	get_bforce_parP3(&srcs,&tzero,&rpars_p3);
	//init_bforceP3(&srcs,&rpars_p3);
	int intsrc = 0;
        float* stfunc = check_malloc(nt*sizeof(float));
	getsource(stfunc,&dt,nt,&srcs.rtime[0],nbounce,intsrc,stype,name,bfilt,&flo,&fhi,0,0,&tdelay,stfdir,stf_file);
	endpar();
	return 0;
}
