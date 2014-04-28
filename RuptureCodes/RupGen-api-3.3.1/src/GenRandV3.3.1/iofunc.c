#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

void _write_field(char *file,struct pointsource *ps,char *type,
		 int nx,int ny,float *dx,float *dy)
{
FILE *fpw;
float xx, yy;
int i, j;

fpw = _fopfile(file,"w");
for(j=0;j<ny;j++)
   {
   yy = (j + 0.5)*(*dy);
   for(i=0;i<nx;i++)
      {
      xx = (i+0.5)*(*dx);

      if(strcmp(type,"slip") == 0)
         fprintf(fpw,"%12.5e %12.5e %12.5e\n",xx,yy,ps[i+j*nx].slip);
      else if(strcmp(type,"rupt") == 0)
         fprintf(fpw,"%12.5e %12.5e %12.5e\n",xx,yy,ps[i+j*nx].rupt);
      }
   }
fclose(fpw);
}

/*
void _write_spec(char *file,float *as,struct complex *slip,int nx,int ny,float *dx,float *dy,float *dkx,float *dky,float *xl,float *yl,int kflag)
{
FILE *fpw;
float kx, ky, amp, amp0, lamp;
float xl2, yl2, fac;
int i, j, ip;

float hcoef = 1.8;  // H=0.8, hcoef = H + 1 

xl2 = (*xl)*(*xl);
yl2 = (*yl)*(*yl);

fft2d(slip,nx,ny,-1,dx,dy);

fpw = _fopfile(file,"w");

// this will normalize max spectrum to one 
amp0 = sqrt(slip[0].re*slip[0].re + slip[0].im*slip[0].im);

for(j=0;j<=ny/2;j++)
   {
   if(j<=ny/2)
      ky = j*(*dky);
   else
      ky = (j-ny)*(*dky);

   for(i=0;i<=nx/2;i++)
      {
      if(i<=nx/2)
         kx = i*(*dkx);
      else
         kx = (i-nx)*(*dkx);

      ip = i + j*nx;

      amp = kx*kx*xl2 + ky*ky*yl2;

      // default is somerville scaling 
      fac = 1.0/sqrt(1.0 + amp*amp);

      if(kflag == MAI_FLAG) // mai scaling 
         {
         fac = exp((hcoef)*log(1.0+amp));
         fac = 1.0/sqrt(fac);
	 }

      if(kflag == SOMERVILLE_FLAG)      // somerville scaling 
         fac = 1.0/sqrt(1.0 + amp*amp);

      amp = sqrt(slip[ip].re*slip[ip].re + slip[ip].im*slip[ip].im);
      amp = amp/amp0;

      lamp = -1.e+20;
      if(amp > 0.0)
	 lamp = log10(amp/fac);

      fprintf(fpw,"%13.5e %13.5e %12.5e %12.5e\n",kx,ky,lamp,180*atan(slip[ip].im/slip[ip].re)/3.14159);

      as[ip] = as[ip] + lamp;
      }
   }
fclose(fpw);

fft2d(slip,nx,ny,1,dkx,dky);
}
*/

void _write_avgspec(char *file,float *as,int ns,int nx,int ny,float *dkx,float *dky)
{
FILE *fpw;
float kx, ky, fac;
int i, j, ip;

fac = 1.0/(float)(ns);

fpw = _fopfile(file,"w");

for(j=0;j<=ny/2;j++)
   {
   ky = j*(*dky);

   for(i=0;i<=nx/2;i++)
      {
      kx = i*(*dkx);
      ip = i + j*nx;
      fprintf(fpw,"%13.5e %13.5e %12.5e\n",kx,ky,fac*as[ip]);
      }
   }
fclose(fpw);
}

FILE *_fopfile(char *name,char *mode)
{
FILE *fp;

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s\n", name, mode);
   exit(-1);
   }
return(fp);
}

int _opfile_ro(char *name)
{
int fd;
if ((fd = open (name, RDONLY_FLAGS, 0444)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

int _opfile(char *name)
{
int fd;
if ((fd = open (name, RDWR_FLAGS, 0664)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

int _croptrfile(char *name)
{
int fd;
if ((fd = open (name, CROPTR_FLAGS, 0664)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

int _reed(int fd, void *pntr, int length)
{
int temp;
if ((temp = read(fd, pntr, length)) < length)
   {
   fprintf (stderr, "READ ERROR\n");
   fprintf (stderr, "%d attempted  %d read\n", length, temp);
   exit(-1);
   }
return(temp);
}

int _rite(int fd, void *pntr, int length)
{
int temp;
if ((temp = write(fd, pntr, length)) < length)
   {
   fprintf (stderr, "WRITE ERROR\n");
   fprintf (stderr, "%d attempted  %d written\n", length, temp);
   exit(-1);
   }
return(temp);
}

#ifdef _USE_MEMCACHED
void _mc_add_file(char* filename, char** file_buffer, memcached_st* mst) {
	memcached_return ret;
	uint32_t flags = 0;

	FILE* fpr = _fopfile(filename, "rb");
        //Get file size so we know how big to make the file buffer
        fseek(fpr, 0, SEEK_END);
        long file_length = ftell(fpr);
	fseek(fpr, 0, SEEK_SET);
        *file_buffer = _check_malloc(sizeof(char) * (file_length+1));
	memset(*file_buffer, '\0', file_length+1);
        int part_id = 0;
	long read = fread(*file_buffer, 1, file_length, fpr);
	if (read!=file_length) {
		fprintf(stderr, "READ ERROR - %ld attempted %ld read.\n", file_length, read);
		exit(1);
	}
        char* file_key = _check_malloc(sizeof(char) * (strlen(filename) + 4));
	char buffer[1000*1000];
	//Memcached thinks this is 1 MB
	int ONE_MB = 1000*1000;
	int i=0;
	while (i*ONE_MB < file_length) {
		int to_copy = ONE_MB;
		if (file_length - i*ONE_MB < ONE_MB) {
			to_copy = file_length - i*ONE_MB;
		}
                sprintf(file_key, "%s.%d", filename, i);
		memcpy(buffer, (*file_buffer)+i*ONE_MB, to_copy);

                ret = memcached_set(mst, file_key, strlen(file_key), buffer, strlen(buffer), 0, flags);
                if (ret!=MEMCACHED_SUCCESS) {
                	printf("Caching of %s failed.\n", file_key);
                        printf("%s\n", memcached_strerror(mst, ret));
                }
                i++;
                memset(buffer, '\0', ONE_MB);
	}
        //Commit number of parts
        int num_parts = i;
	char* parts_key = _check_malloc(sizeof(char) * (strlen(filename) + 1));
	strcpy(parts_key, filename);
        ret = memcached_set(mst, parts_key, strlen(parts_key), (char*)(&num_parts), sizeof(int), 0, flags);
        if (ret!=MEMCACHED_SUCCESS) {
                printf("Caching of %s failed.\n", parts_key);
                printf("%s\n", memcached_strerror(mst, ret));
        }
        free(file_key);
	free(parts_key);
	fclose(fpr);
}


struct pointsource* _mc_read_ruppars(char *file,struct pointsource *psrc,float *mag,int *nx,int *ny,float *dx,float *dy,float *dtop,float *stk,float *dip,float *elon,float *elat, char* mc_server)
{
float area;
int i, nn;
char str[1024];

double rperd = 0.017453293;
int ONE_MB = 1000*1000; //This is what memcached thinks 1 MB is

*dtop = 1.0e+15;
*stk = 0.0;
*dip = 0.0;
*elon = 0.0;
*elat = 0.0;

char* buffer;
char* file_data;
int* num_parts;

	printf("Using %s as memcached server.\n", mc_server);

	memcached_return ret;
	memcached_st* mst;

	mst = memcached_create(NULL);
	memcached_server_st *server = memcached_server_list_append(NULL, mc_server, 11211, &ret);

	ret = memcached_server_push(mst, server);
	memcached_server_free(server);

	uint32_t flags = 0;

	//See if file parts is in cache
	//key is filename + "_" + part
	char* parts_key = _check_malloc(sizeof(char) * (strlen(file) + 1));
	strcpy(parts_key, file);
	
	size_t incoming_size;
	num_parts = (int*) memcached_get(mst, parts_key, strlen(parts_key), &incoming_size, &flags, &ret);
	if (num_parts==NULL) {
                printf("Key %s hasn't been cached, adding.\n", parts_key);
		_mc_add_file(file, &file_data, mst);
	} else {
		printf("Found key %s, num parts is %d.\n", parts_key, *num_parts);
		printf("Looking for parts.\n");
		//See if all the parts are actually around
		file_data = _check_malloc(sizeof(char) * (ONE_MB * *num_parts + 1));
		memset(file_data, '\0', ONE_MB* *num_parts + 1);
		char* file_key = _check_malloc(sizeof(char) * (strlen(file) + 4));
		size_t incoming;
		for (i=0; i<*num_parts; i++) {
			sprintf(file_key, "%s.%d", file, i);
			buffer = (char*) memcached_get(mst, file_key, strlen(file_key), &incoming, &flags, &ret);
			if (buffer==NULL) {
				printf("Piece %s was missing, adding file %s.\n", file_key, file);
				//This piece is missing
				//Easiest to add the whole file
				free(file_data);
				_mc_add_file(file, &file_data, mst);
				break;
			} else {
				memcpy(file_data+i*ONE_MB, buffer, incoming);
			}
		}
		free(file_key);
	}

	free(parts_key);

	char* tok;
	tok = strtok(file_data, "\n"); /* Probability = <float> */

	tok = strtok(NULL, "\n"); /* Magnitude = <float> */
	sscanf(tok,"%*s %*s %f",mag);

	tok = strtok(NULL, "\n");   /* GridSpacing = <float> */
	sscanf(tok,"%*s %*s %f",dx);

	if(*dx == (float)(0.0))
	   {
	   fprintf(stderr,"***** input error\n");
	   fprintf(stderr,"      GridSpacing = 0.0, exiting...\n");
	   exit(-1);
	   }

	*dy = *dx;

	tok = strtok(NULL, "\n");   /* NumRows = <int> */
	sscanf(tok,"%*s %*s %d",ny);

        tok = strtok(NULL, "\n");   /* NumCols = <int> */
	sscanf(tok,"%*s %*s %d",nx);

        tok = strtok(NULL, "\n");/* header comment */
	psrc = (struct pointsource *)_check_realloc(psrc,(*nx)*(*ny)*sizeof(struct pointsource));

	area = (*dx)*(*dy)*1.0e+10;  /* km -> cm */

	for(i=0;i<(*nx)*(*ny);i++)
	   {
           tok = strtok(NULL, "\n");   /* Lat , Lon , Depth , Rake , Dip , Strike */
	   sscanf(tok,"%f %f %f %f %f %f",&psrc[i].lat,
					  &psrc[i].lon,
					  &psrc[i].dep,
					  &psrc[i].rak,
					  &psrc[i].dip,
					  &psrc[i].stk);

	   psrc[i].area = area;

	   if(psrc[i].dep < *dtop)
	      *dtop = psrc[i].dep;

	   *stk = *stk + psrc[i].stk;
	   *dip = *dip + psrc[i].dip;
	   }
	*stk = *stk/((*nx)*(*ny));
	*dip = *dip/((*nx)*(*ny));

	nn = 0;
	for(i=0;i<(*nx)*(*ny);i++)
	   {
	   if(psrc[i].dep < (*dtop + 0.01))
	      {
	      *elon = *elon + psrc[i].lon;
	      *elat = *elat + psrc[i].lat;
	      nn++;
	      }
	   }

	/* adjust for half subfault width */
	*dtop = (*dtop) - 0.5*(*dx)*sin((*dip)*rperd);

	if(nn == 0)
	   nn++;

	*elon = *elon/nn;
	*elat = *elat/nn;

	free(file_data);

	return(psrc);
	}
#endif


	struct pointsource *_read_ruppars(char *file,struct pointsource *psrc,float *mag,int *nx,int *ny,float *dx,float *dy,float *dtop,float *stk,float *dip,float *elon,float *elat)
	{
	FILE *fpr;
	float area;
	int i, nn;
	char str[1024];

	double rperd = 0.017453293;

	*dtop = 1.0e+15;
	*stk = 0.0;
	*dip = 0.0;
	*elon = 0.0;
	*elat = 0.0;

	if(strcmp(file,"stdin") == 0)
	   fpr = stdin;
	else
	   fpr = _fopfile(file,"r");

	fgets(str,1024,fpr);   /* Probability = <float> */

	fgets(str,1024,fpr);   /* Magnitude = <float> */
	sscanf(str,"%*s %*s %f",mag);

	fgets(str,1024,fpr);   /* GridSpacing = <float> */
	sscanf(str,"%*s %*s %f",dx);

	if(*dx == (float)(0.0))
	   {
	   fprintf(stderr,"***** input error\n");
	   fprintf(stderr,"      GridSpacing = 0.0, exiting...\n");
	   exit(-1);
	   }

	*dy = *dx;

	fgets(str,1024,fpr);   /* NumRows = <int> */
	sscanf(str,"%*s %*s %d",ny);

	fgets(str,1024,fpr);   /* NumCols = <int> */
	sscanf(str,"%*s %*s %d",nx);

	fgets(str,1024,fpr);   /* header comment */

	psrc = (struct pointsource *)_check_realloc(psrc,(*nx)*(*ny)*sizeof(struct pointsource));

	area = (*dx)*(*dy)*1.0e+10;  /* km -> cm */

	for(i=0;i<(*nx)*(*ny);i++)
	   {
	   fgets(str,1024,fpr);   /* Lat , Lon , Depth , Rake , Dip , Strike */
	   sscanf(str,"%f %f %f %f %f %f",&psrc[i].lat,
					  &psrc[i].lon,
					  &psrc[i].dep,
					  &psrc[i].rak,
					  &psrc[i].dip,
					  &psrc[i].stk);

	   psrc[i].area = area;

	   if(psrc[i].dep < *dtop)
	      *dtop = psrc[i].dep;

	   *stk = *stk + psrc[i].stk;
	   *dip = *dip + psrc[i].dip;
	   }
	fclose(fpr);

	*stk = *stk/((*nx)*(*ny));
	*dip = *dip/((*nx)*(*ny));

	nn = 0;
	for(i=0;i<(*nx)*(*ny);i++)
	   {
	   if(psrc[i].dep < (*dtop + 0.01))
	      {
	      *elon = *elon + psrc[i].lon;
	      *elat = *elat + psrc[i].lat;
	      nn++;
	      }
	   }

	/* adjust for half subfault width */
	*dtop = (*dtop) - 0.5*(*dx)*sin((*dip)*rperd);

	if(nn == 0)
	   nn++;

	*elon = *elon/nn;
	*elat = *elat/nn;

	return(psrc);
	}

	struct pointsource *_set_ruppars(struct pointsource *psrc,float *mag,int *nx,int *ny,float *dx,float *dy,float *dtop,float *stk,float *dip,float *rak,float *elon,float *elat)
	{
	float area;
	float cosA, sinA, cosD, sinD, fwid, flen;
	float xx, yy, zz, dd, sn, se;
	int i, j, ip;

	double rperd = 0.017453293;

	psrc = (struct pointsource *)_check_realloc(psrc,(*nx)*(*ny)*sizeof(struct pointsource));

	area = (*dx)*(*dy)*1.0e+10;  /* km -> cm */
	flen = (*nx)*(*dx);
	fwid = (*ny)*(*dy);

	cosA = cos((*stk)*rperd);
	sinA = sin((*stk)*rperd);
	cosD = cos((*dip)*rperd);
	sinD = sin((*dip)*rperd);
	for(j=0;j<(*ny);j++)
	   {
	   dd = (j + 0.5)*(*dy);
	   yy = dd*cosD;
	   zz = (*dtop) + dd*sinD;

	   for(i=0;i<(*nx);i++)
	      {
	      ip = i + j*(*nx);
	      xx = (i+0.5)*(*dx) - 0.5*flen;

	      se = xx*sinA + yy*cosA;
	      sn = xx*cosA - yy*sinA;
	      _set_ll(elon,elat,&psrc[ip].lon,&psrc[ip].lat,&sn,&se);

	      psrc[ip].dep = zz;
	      psrc[ip].stk = (*stk);
	      psrc[ip].dip = (*dip);
	      psrc[ip].rak = (*rak);
	      psrc[ip].area = area;
	      }
	   }

	return(psrc);
	}

	struct pointsource *_read_gsfpars(char *file,struct pointsource *psrc,struct generic_slip *gslip,float *dx,float *dy,float *dtop,float *dip)
	{
	FILE *fpr;
	int i, nn;
	char str[1024];
	struct slippars *spar;

	double rperd = 0.017453293;

	*dtop = 1.0e+15;
	*dx = 0.0;
	*dy = 0.0;
	*dip = 0.0;

	if(strcmp(file,"stdin") == 0)
	   fpr = stdin;
	else
	   fpr = _fopfile(file,"r");

	fgets(str,1024,fpr);
	while(strncmp(str,"#",1) == 0)
	   fgets(str,1024,fpr);

	sscanf(str,"%d",&gslip->np);

	psrc = (struct pointsource *)_check_realloc(psrc,gslip->np*sizeof(struct pointsource));
	gslip->spar = (struct slippars *)_check_realloc(gslip->spar,gslip->np*sizeof(struct slippars));
	spar = gslip->spar;

	i = 0;
	while(fgets(str,1024,fpr) != NULL)
	   {
	   sscanf(str,"%f %f %f %f %f %f %f %f %f %f %d",&spar[i].lon,
					     &spar[i].lat,
					     &spar[i].dep,
					     &spar[i].ds,
					     &spar[i].dw,
					     &spar[i].stk,
					     &spar[i].dip,
					     &spar[i].rake,
					     &spar[i].slip,
					     &spar[i].tinit,
					     &spar[i].segno);

	   psrc[i].lon = spar[i].lon;
	   psrc[i].lat = spar[i].lat;
	   psrc[i].dep = spar[i].dep;
	   psrc[i].stk = spar[i].stk;
	   psrc[i].dip = spar[i].dip;
	   psrc[i].rak = spar[i].rake;
	   psrc[i].area = spar[i].ds*spar[i].dw*1.0e+10;

	   if(psrc[i].dep < *dtop)
	      *dtop = psrc[i].dep;

	   *dip = *dip + psrc[i].dip;
	   *dx = *dx + spar[i].ds;
	   *dy = *dy + spar[i].dw;

	   i++;
	   }
	fclose(fpr);

	*dip = *dip/(gslip->np);
	*dx = *dx/(gslip->np);
	*dy = *dy/(gslip->np);

	/* adjust for half subfault width */
	*dtop = (*dtop) - 0.5*(*dy)*sin((*dip)*rperd);
	if(*dtop < 0.0)
	   *dtop = 0.0;

	return(psrc);
	}


	int _cp(const char *to, const char *from) 
	{ 
	    FILE *fp_from, *fp_to;
	    char buf[32768]; 
	    ssize_t nread; 
	 
	    fp_from = fopen(from, "r");
	    if (fp_from == NULL) {
	      return(1);
	    }

	    fp_to = fopen(to, "w");
	    if (fp_to == NULL) {
	      fclose(fp_from);
	      return(1);
	    }

	    while (nread = fread(buf, 1, sizeof(buf), fp_from), 
		   nread > 0)  { 
	      char *out_ptr = buf; 
	      ssize_t nwritten; 
	 
	      do { 
		nwritten = fwrite(out_ptr, 1, nread, fp_to); 	
		if (nwritten == nread) { 
		  nread -= nwritten; 
		  out_ptr += nwritten; 
		} else {
		  fclose(fp_from);
		  fclose(fp_to);
		  return(1);
		} 
	      } while (nread > 0); 
	    }
	    
	    fclose(fp_from);
	    fclose(fp_to);
	    return(0);
	} 


	/* Check if file exists */
	int _file_exists(const char *file)
	{
	  struct stat st;

	  if (stat(file, &st) == 0) {
	    return(1);
	  } else {
	    return(0);
	  }
	}


	/* Returns true if path is a file */
	int _rg_is_file(const char *path)
	{
	  struct stat st;

	  if (stat(path, &st) == 0) {
	    if ((S_ISREG(st.st_mode)) && (!S_ISDIR(st.st_mode))) {
	      return(1);
	    } else {
	      return(0);
	    }
	  }

	  return(0);
	}


	/* Safe string copy */
	int _rg_strcpy(char *str1, const char *str2, int str1len)
	{
	  if (str1 == NULL) {
	    return(1);
	  }

	  if (str2 == NULL) {
	    strcpy(str1, "");
	    return(1);
	  }

	  if (snprintf(str1, str1len, "%s", str2) > str1len) {
	    fprintf(stderr, "Warning (ucvm_strcpy): String %s truncated to %s\n",
		    str2, str1);
	  }
	  return(0);
	}

	void _write2gsf(struct generic_slip *gslip,struct pointsource *ps,char *ifile,char *ofile)
	{
	FILE *fpr, *fpw;
	int ip;
	char str[1024];

	if(strcmp(ofile,"stdout") == 0)
	   fpw = stdout;
	else
	   fpw = _fopfile(ofile,"w");

	fpr = _fopfile(ifile,"r");
	fgets(str,1024,fpr);
	while(strncmp(str,"#",1) == 0)
	   {
	   fprintf(fpw,"%s",str);
	   fgets(str,1024,fpr);
	   } 
	fclose(fpr);
	   
	fprintf(fpw,"%d\n",gslip->np);
	      
	for(ip=0;ip<gslip->np;ip++)
	   {
	   fprintf(fpw,"%11.5f %11.5f %8.4f %8.4f %8.4f %6.1f %6.1f %6.1f %8.2f %8.3f %3d\n",
								       gslip->spar[ip].lon,
								       gslip->spar[ip].lat,
								       gslip->spar[ip].dep,
								       gslip->spar[ip].ds,
								       gslip->spar[ip].dw,
								       gslip->spar[ip].stk,
								       gslip->spar[ip].dip,
								       ps[ip].rak,
								       ps[ip].slip,
								       ps[ip].rupt,
								       gslip->spar[ip].segno);
	   }
	fclose(fpw);
	}


	//#define MAX_FWRITE_BUF 10000000
	//char fwrite_buf[MAX_FWRITE_BUF];
	//int fwrite_len = 0;

	//int fwrite_buffered(FILE *fd, char *pntr, int length)
	//{
	//  if (length > MAX_FWRITE_BUF) {
	//    fprintf(stderr, "String exceeds max buffer size\n");
	//    exit(-1);
	//  }

	//  if ((length + fwrite_len) > MAX_FWRITE_BUF) {
	//    /* Flush buffer */
	//    if (fwrite(fwrite_buf, sizeof(char), fwrite_len, fd) != fwrite_len) {
	//      fprintf(stderr, "Buffered write failure\n");
	//      exit(-1);
	//    }
	//    fwrite_len = 0;
	//  }
	//
	//  memcpy(fwrite_buf + fwrite_len, pntr, length);
//  fwrite_len += length;
//  return(0);
//}

//int fwrite_flush(FILE *fd) 
//{
//  if (fwrite_len > 0) {
//    /* Flush buffer */
//    if (fwrite(fwrite_buf, sizeof(char), fwrite_len, fd) != fwrite_len) {
//      fprintf(stderr, "Buffered flush failure\n");
//    exit(-1);
//    }
//    fwrite_len = 0;
//  }
//
//  return(0);
//}
