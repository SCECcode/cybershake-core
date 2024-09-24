#include        "include.h"
#include        "defs.h"
#include        "structure.h"
#include        "functions.h"

FILE *fopfile(char *name,char *mode)
{
FILE *fp;

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s\n", name, mode);
   exit(-1);
   }
//printf("Opened %s as fp %d\n", name, fp);
return(fp);
}

void fopfile_ro(char* name, FILE** fp) {
   if ((*fp = fopen(name, "rb"))==NULL) {
	fprintf(stderr,"CAN'T FOPEN FILE = %s for reading.\n", name);
	exit(-1);
   }
}

int opfile_ro(char *name)
{
int fd;
if ((fd = open (name, O_RDONLY, 0444)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
//printf("Opened %s as fp %d\n", name, fd);
return (fd);
}

int opfile(char *name)
{
int fd;
if ((fd = open (name, O_RDWR, 0664)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
//printf("Opened %s as fp %d\n", name, fd);
return (fd);
}

void fcroptrfile(char* name, FILE** fp) {
   if ((*fp = fopen(name, "w"))==NULL) {
	fprintf(stderr, "CAN'T FOPEN FILE = %s for writing.\n", name);
	exit(-1);
   }
}

int croptrfile(char *name)
{
int fd;
if ((fd = open (name, O_CREAT | O_TRUNC | O_RDWR, 0664)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
//printf("Opened %s as fp %d\n", name, fd);
return (fd);
}

int freed(FILE* fp, void* pntr, size_t record_size, size_t num_records) {
   size_t num_read;
   if ((num_read = fread(pntr, record_size, num_records, fp))!=num_records) {
	fprintf(stderr, "FREAD ERROR\n");
	fprintf(stderr, "%ld records attempted %ld read\n", num_records, num_read);
	fprintf(stderr, "%s\n", strerror(errno));
	exit(-1);
   }
   return num_read*record_size;
}

int reed(int fd, void *pntr, int length)
{
int temp;
if ((temp = read(fd, pntr, length)) < length)
   {
   fprintf (stderr, "READ ERROR\n");
   fprintf (stderr, "%d attempted  %d read\n", length, temp);
   fprintf (stderr, "%s\n", strerror(errno));
   exit(-1);
   }
//printf("Reading %d bytes from fd %d\n", length, fd);
return(temp);
}

int frite(FILE* fp, void *pntr, size_t record_size, size_t num_records) {
   size_t num_write;
   if ((num_write = fwrite(pntr, record_size, num_records, fp))!=num_records) {
        fprintf(stderr, "FWRITE ERROR\n");
        fprintf(stderr, "%ld records attempted %ld written\n", num_records, num_write);
        fprintf(stderr, "%s\n", strerror(errno));
        exit(-1);
   }
   return num_records*record_size;
}

int rite(int fd, void *pntr, int length)
{
int temp;
if ((temp = write(fd, pntr, length)) < length)
   {
   fprintf (stderr, "WRITE ERROR\n");
   fprintf (stderr, "%d attempted  %d written\n", length, temp);
   exit(-1);
   }
//printf("Wrote %d bytes to fd %d\n", length, fd);
return(temp);
}

void *check_realloc(void *ptr,size_t len)
{
if (debug) {
    char buf[256];
    sprintf(buf, "Reallocating %ld bytes.", len);
    write_log(buf);
}

ptr = realloc(ptr,len);
//fprintf(stderr,"Reallocing %ld bytes to pointer %ld.\n", len, ptr);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory reallocation error\n");
   fprintf(stderr,"*****  memory error in process %d\n", my_global_id);
   fprintf(stderr,"Tried to realloc %ld bytes.\n", len);
   void* array[20];
   size_t size;
   char** strings;
   backtrace(array, 20);
   strings = backtrace_symbols(array, size);
   int i;
   for (i=0; i<size; i++) {
        fprintf(stderr, "stacktrace %d) %s\n", i, strings[i]);
   }
   fflush(stderr);
   exit(-1);
   }

return(ptr);
}

void *check_malloc(size_t len)
{
char *ptr;

if (debug) {
	char buf[256];
	sprintf(buf, "Allocating %ld bytes.", len);
	write_log(buf);
}
//fprintf(stderr,"%d) Allocating %ld bytes.\n", my_global_id, len);

//ptr = (char *) malloc (len);
//
ptr = malloc(len);
 
if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory allocation error\n");
   fprintf(stderr,"*****  memory error in process %d\n", my_global_id);
   fprintf(stderr,"Tried to allocate %ld bytes.\n", len);
   fflush(stderr);
   if (debug) {
		char buf[256];
		sprintf(buf,"*****  memory error in process %d\n", my_global_id);
		write_log(buf);
		sprintf(buf,"Tried to allocate %ld bytes.\n", len);
		write_log(buf);
		close_log();
	}
   void* array[20];
   size_t size;
   char** strings;
   backtrace(array, 20);
   strings = backtrace_symbols(array, size);
   int i;
   for (i=0; i<size; i++) {
	fprintf(stderr, "stacktrace %d) %s\n", i, strings[i]);
   }
   fflush(stderr);
   MPI_Abort(MPI_COMM_WORLD, -1);
   exit(-1);
   }
 
return(ptr);
}

void fortran_rite(int fd,int nargs, ...)
{
va_list ap;
void *ptr[MAX_VAR_LIST];
int len[MAX_VAR_LIST];
int totlen = 0;
int i;

va_start(ap,nargs);
for(i=0;i<nargs;i++)
   {
   ptr[i] = va_arg(ap,void *);
   len[i] = va_arg(ap,int);
   totlen = totlen + len[i];
   }
va_end(ap);

rite(fd,&totlen,sizeof(int));

for(i=0;i<nargs;i++)
   rite(fd,ptr[i],len[i]);

rite(fd,&totlen,sizeof(int));
}

void write_seis_ascii(char* sname, char* comp, float *st, float *dt, int nt, float *ts, char* filename) {
  char stitle[128], header[128];
  int i, j, nt6;
  char* full_filename = malloc(sizeof(char) * (strlen(filename) + strlen(comp) + 2));  full_filename = strcpy(full_filename, filename);
  full_filename = strcat(full_filename, ".");
  full_filename = strcat(full_filename, comp);
  FILE* fpw = fopfile(full_filename,"w");

  sprintf(stitle,"%-10s%3s %s\n",sname,comp,"TITLE");
  fprintf(fpw,"%s",stitle);

  sprintf(header,"%d %12.5e 0 0 %12.5e 0.0 0.0 0.0\n",nt,*dt,*ts);
  fprintf(fpw,"%s",header);

  for (i=0; i<nt; i++) {
	fprintf(fpw,"%13.5e",st[i]);
	if ((i+1)%6==0) {
	  fprintf(fpw,"\n");
	}
  }
  if (nt%6!=0) {
	fprintf(fpw,"\n");
  }
  fflush(fpw);
  fclose(fpw);
  free(full_filename);
}

void write_seis_binary(struct seisheader* header, float* st, float* dt, int nt, float* ts, char* filename) {
  int written;
  char outfile[256];
  FILE* fpout;
  //Combine header and data into single buffer so we only have 1 write
  int header_size = sizeof(struct seisheader);
  int to_write = header_size + nt * sizeof(float);
  char* buffer = check_malloc(to_write);
  memcpy(buffer, header, header_size);
  memcpy(buffer+header_size, st, nt*sizeof(float));
  //sprintf(outfile, "%s.%s", filename, comp);
  //Write to file
  fpout = fopfile(filename, "wb");
  written = fwrite(buffer, 1, to_write, fpout);
  if (written!=to_write) {
	fprintf(stderr, "Error in writing to file %s; wrote %d instead of %d bytes.\n", filename, written, to_write);
        exit(-1);
  }
  fflush(fpout);
  fclose(fpout);
  free(buffer);
}

void write_seis(struct seisheader* header, char* filename, char *sname,char *comp,float *st,float *dt,int nt,float *ts, int output_binary) {
  FILE *fopfile(), *fpw;
  int i, j, nt6;
  //char outfile[256];
  
  //sprintf(outfile,"%s.%s",dir,stat,comp);
  
  if (!output_binary) {
	write_seis_ascii(sname, comp, st, dt, nt, ts, filename);
  } else {
	write_seis_binary(header, st, dt, nt, ts, filename);
  }
}


char *skipval(int j,char *str)
{
while(str[0] == ' ' || str[0] == '\t' || str[0] == '\b' || str[0] == '\n')
   str++;

while(j--)
   {
   while(str[0] != ' ' && str[0] != '\t' && str[0] != '\b' && str[0] != '\n')
      str++;

   while(str[0] == ' ' || str[0] == '\t' || str[0] == '\b' || str[0] == '\n')
      str++;
   }

return(str);
}

void makedir(char *ipath)
{
struct stat sbuf;
char stmp[256], str[128], path[1024];
int rtn, j;
mode_t mode = 00777;

strcpy(path,ipath);

j = 0;
while(path[j] != '\0')
   j++;

j--;
while(path[j] == '/')
   j--;
path[j+1] = '\0';
path[j+2] = '\0';

j = 0;
while(path[j] != '\0')
   {
   while(path[j] != '/' && path[j] != '\0')
      j++;

   if(j != 0)
      {
      strncpy(stmp,path,j);
      stmp[j] = '\0';

      rtn = stat(stmp,&sbuf); /* stat directory path to see if it already exists */

      if(rtn == -1 && errno == ENOENT) /* try to make the directory path */
         {
         rtn = mkdir(stmp,mode);

         if(rtn == -1)
            {
            if(errno != EEXIST)
               {
               sprintf(str,"makedir() cannot make directory %s, exiting",stmp);
               perror(str);
               exit(-1);
               }
            }
         }

      else if(rtn == -1 && errno != ENOENT) /* some other problem */
         {
         sprintf(str,"problem with stat() on %s, exiting",stmp);
         perror(str);
         exit(-1);
         }
      }
   j++;
   }

/*
   Double check to make sure directory exists.  This is a brute-force
   method, but I ran inot problems with automounted directories using the
   error-checking above.  RWG 9/20/99
*/

rtn = mkdir(stmp,mode);
if(rtn == -1 && errno != EEXIST)
   {
   sprintf(str,"makedir() cannot make directory %s, exiting",stmp);
   perror(str);
   exit(-1);
   }

}

int parse_rup_geom(char* rup_geom_file, struct rup_geom_point** rg_points) {
        FILE* fp_in;
        int i, num_rows, num_cols;
        int num_pts = -1;
        float rake, dip, strike;
        fopfile_ro(rup_geom_file, &fp_in);
        char line[256];
	//Skip first 3 rows
	for (i=0; i<3; i++) {
                fgets(line, 256, fp_in);
        }
	//Next row has rows
	fscanf(fp_in, "NumRows = %d\n", &num_rows);
	//Then cols
        fscanf(fp_in, "NumCols = %d\n", &num_cols);
        num_pts = num_rows * num_cols;
        //Then comment lines
        fgets(line, 256, fp_in);
	while (line[0] == '#') {
                fgets(line, 256, fp_in);
        }
	//Create buffers
	*rg_points = check_malloc(num_pts * sizeof(struct rup_geom_point));
        //Read lines and insert
        for (i=0; i<num_pts-1; i++) {
                sscanf(line, "%f %f %f %f %f %f\n", &((*rg_points)[i].lat), &((*rg_points)[i].lon), &((*rg_points)[i].dep), &rake, &dip, &strike);
                fgets(line, 256, fp_in);
        }
	//add last line
	sscanf(line, "%f %f %f %f %f %f\n", &((*rg_points)[num_pts-1].lat), &((*rg_points)[num_pts-1].lon), &((*rg_points)[num_pts-1].dep), &rake, &dip, &strike);
        fclose(fp_in);
        return num_pts;
}
