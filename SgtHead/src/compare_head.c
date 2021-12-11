#include "include.h"
#include "function.h"

/* This function does a field-by-field comparison of an emod3d-generated header
* and a write_head-generated header.
*/

int compare_sgtmaster(int sgt_in, int head_in) {
	struct sgtmaster mast1, mast2;
	reed(sgt_in, &mast1, sizeof(struct sgtmaster));
	reed(head_in, &mast2, sizeof(struct sgtmaster));
	if (mast1.geoproj!=mast2.geoproj) {
		printf("Difference in geoproj: %d vs %d\n", mast1.geoproj, mast2.geoproj);
	}
	if (mast1.modellon!=mast2.modellon) {
                printf("Difference in modellon: %f vs %f\n", mast1.modellon, mast2.modellon);
	}
	if (mast1.modellat!=mast2.modellat) {
                printf("Difference in modellat: %f vs %f\n", mast1.modellat, mast2.modellat);
	}
	if (mast1.modelrot!=mast2.modelrot) {
                printf("Difference in modelrot: %f vs %f\n", mast1.modelrot, mast2.modelrot);
	}
	if (mast1.xshift!=mast2.xshift) {
               printf("Difference in xshift: %f vs %f\n", mast1.xshift, mast2.xshift);
	}
	if (mast1.yshift!=mast2.yshift) {
                printf("Difference in yshift: %f vs %f\n", mast1.yshift, mast2.yshift);
	}
	if (mast1.globnp!=mast2.globnp) {
                printf("Difference in globnp: %d vs %d\n", mast1.globnp, mast2.globnp);
	}
	if (mast1.localnp!=mast2.localnp) {
                printf("Difference in localnp: %d vs %d\n", mast1.localnp, mast2.localnp);
	}
	if (mast1.nt!=mast2.nt) {
                printf("Difference in nt: %d vs %d\n", mast1.nt, mast2.nt);
	}
	return mast1.globnp;
}

void compare_sgtindex(int sgt_in, int head_in, int num_pts) {
	struct sgtindex indx1, indx2;
	int i;
	printf("Comparing %d points.\n", num_pts);
	for (i=0; i<num_pts; i++) {
		reed(sgt_in, &indx1, sizeof(struct sgtindex));
		reed(head_in, &indx2, sizeof(struct sgtindex));
		if (indx1.indx!=indx2.indx) {
                	printf("Difference in indx: %lld vs %lld\n", indx1.indx, indx2.indx);
        	}
                if (indx1.xsgt!=indx2.xsgt) {
                        printf("Difference in xsgt: %d vs %d\n", indx1.xsgt, indx2.xsgt);
                }
                if (indx1.ysgt!=indx2.ysgt) {
                        printf("Difference in ysgt: %d vs %d\n", indx1.ysgt, indx2.ysgt);
                }
                if (indx1.zsgt!=indx2.zsgt) {
                        printf("Difference in zsgt: %d vs %d\n", indx1.zsgt, indx2.zsgt);
                }
                if (indx1.h!=indx2.h) {
                        printf("Difference in h: %f vs %f\n", indx1.h, indx2.h);
                }
	}
}

void compare_sgtheader(int sgt_in, int head_in, int num_pts, int nt, int stand_alone) {
	struct sgtheader head1, head2;
	int i;
	for (i=0; i<num_pts; i++) {
//	for (i=0; i<10; i++) {
		reed(sgt_in, &head1, sizeof(struct sgtheader));
		reed(head_in, &head2, sizeof(struct sgtheader));
		//printf("%d %d %d\n", head1.xsgt, head1.ysgt, head1.zsgt);
                //lseek(sgt_in, sizeof(float)*6*nt, SEEK_CUR);
		//continue;
		if (head1.indx!=head2.indx) {
			printf("Difference in indx: %lld vs %lld\n", head1.indx, head2.indx);
		}
                if (head1.geoproj!=head2.geoproj) {
                        printf("Difference in geoproj: %d vs %d\n", head1.geoproj, head2.geoproj);
                }
                if (head1.modellon!=head2.modellon) {
                        printf("Difference in modellon: %f vs %f\n", head1.modellon, head2.modellon);
                }
                if (head1.modellat!=head2.modellat) {
                        printf("Difference in modellat: %f vs %f\n", head1.modellat, head2.modellat);
                }
                if (head1.modelrot!=head2.modelrot) {
                        printf("Difference in modelrot: %f vs %f\n", head1.modelrot, head2.modelrot);
                }
                if (head1.xshift!=head2.xshift) {
                        printf("Difference in xshift: %f vs %f\n", head1.xshift, head2.xshift);
                }
                if (head1.yshift!=head2.yshift) {
                        printf("Difference in yshift: %f vs %f\n", head1.yshift, head2.yshift);
                }
                if (head1.nt!=head2.nt) {
                        printf("Difference in nt: %d vs %d\n", head1.nt, head2.nt);
                }
                if (head1.xazim!=head2.xazim) {
                        printf("Difference in xazim: %f vs %f\n", head1.xazim, head2.xazim);
                }
                if (head1.dt!=head2.dt) {
                        printf("Difference in dt: %f vs %f\n", head1.dt, head2.dt);
                }
                if (head1.tst!=head2.tst) {
                        printf("Difference in tst: %f vs %f\n", head1.tst, head2.tst);
                }
                if (head1.h!=head2.h) {
                        printf("Difference in h: %f vs %f\n", head1.h, head2.h);
                }
                if (head1.src_lat!=head2.src_lat) {
                        printf("Difference in src_lat: %f vs %f\n", head1.src_lat, head2.src_lat);
                }
                if (head1.src_lon!=head2.src_lon) {
                        printf("Difference in src_lon: %f vs %f\n", head1.src_lon, head2.src_lon);
                }
                if (head1.src_dep!=head2.src_dep) {
                        printf("Difference in src_dep: %f vs %f\n", head1.src_dep, head2.src_dep);
                }
                if (head1.xsrc!=head2.xsrc) {
                        printf("Difference in xsrc: %d vs %d\n", head1.xsrc, head2.xsrc);
                }
                if (head1.ysrc!=head2.ysrc) {
                        printf("Difference in ysrc: %d vs %d\n", head1.ysrc, head2.ysrc);
                }
                if (head1.zsrc!=head2.zsrc) {
                        printf("Difference in zsrc: %d vs %d\n", head1.zsrc, head2.zsrc);
                }
                if (fabs(head1.sgt_lat-head2.sgt_lat)>0.00001) {
                        printf("Difference in sgt_lat: %f vs %f\n", head1.sgt_lat, head2.sgt_lat);
                }
                if (fabs(head1.sgt_lon-head2.sgt_lon)>0.00001) {
                        printf("Difference in sgt_lon: %f vs %f\n", head1.sgt_lon, head2.sgt_lon);
                }
                if (head1.sgt_dep!=head2.sgt_dep) {
                        printf("Difference in sgt_dep: %f vs %f\n", head1.sgt_dep, head2.sgt_dep);
                }
                if (head1.xsgt!=head2.xsgt) {
                        printf("Difference in xsgt: %d vs %d\n", head1.xsgt, head2.xsgt);
                }
                if (head1.ysgt!=head2.ysgt) {
                        printf("Difference in ysgt: %d vs %d\n", head1.ysgt, head2.ysgt);
                }
                if (head1.zsgt!=head2.zsgt) {
                        printf("Difference in zsgt: %d vs %d\n", head1.zsgt, head2.zsgt);
                }
                if (head1.cdist!=head2.cdist) {
                        printf("Difference in cdist: %f vs %f\n", head1.cdist, head2.cdist);
                }
                if (head1.lam!=head2.lam) {
                        printf("Difference in lam: %f vs %f\n", head1.lam, head2.lam);
                }
		/*float mu_diff = head1.mu - head2.mu;
                if (fabs(mu_diff/head1.mu) > 0.00001) {
                        printf("Difference in mu: %f vs %f\n", head1.mu, head2.mu);
                }
		float rho_diff = head1.rho - head2.rho;
                if (fabs(rho_diff/head1.rho)>0.0002) {
                        printf("Point %d %d %d: difference in rho: %f vs %f\n", head1.xsgt, head1.ysgt, head1.zsgt, head1.rho, head2.rho);
                }*/
                if (head1.mu!=head2.mu) {
                        printf("Difference in mu: %f vs %f\n", head1.mu, head2.mu);
                }
                if (fabs(head1.rho-head2.rho)>0.00001) {
                        printf("Point %d %d %d: difference in rho: %f vs %f\n", head1.xsgt, head1.ysgt, head1.zsgt, head1.rho, head2.rho);
		}

                if (head1.xmom!=head2.xmom) {
                        printf("Difference in xmom: %f vs %f\n", head1.xmom, head2.xmom);
                }
                if (head1.ymom!=head2.ymom) {
                        printf("Difference in ymom: %f vs %f\n", head1.ymom, head2.ymom);
                }
                if (head1.zmom!=head2.zmom) {
                        printf("Difference in zmom: %f vs %f\n", head1.zmom, head2.zmom);
                }
		//skip ahead record bytes in SGT file
		if (stand_alone==1) {
			lseek(sgt_in, sizeof(float)*6*nt, SEEK_CUR);
		}
		if (stand_alone==0) {
                        lseek(sgt_in, sizeof(float)*6*nt, SEEK_CUR);
			lseek(head_in, sizeof(float)*6*nt, SEEK_CUR);
		}
	}
}

int main(int argc, char** argv) {
	if (argc<4) {
		printf("Usage: %s <SGT file> <header file> <nt_saved> <mode>\n", argv[0]);
		exit(1);
	}
	
	int sgt_in, head_in, nt;
	sgt_in = opfile_ro(argv[1]);
	head_in = opfile_ro(argv[2]);
	nt = atoi(argv[3]);
	int stand_alone=1;

	if (argc==5) {
		stand_alone=atoi(argv[4]);
	}
	printf("Stand alone: %d\n", stand_alone);

	int num_pts = compare_sgtmaster(sgt_in, head_in);
	compare_sgtindex(sgt_in, head_in, num_pts);
	compare_sgtheader(sgt_in, head_in, num_pts, nt, stand_alone);
	return 0;
}
