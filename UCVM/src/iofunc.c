#include        "include.h"
#include        "function.h"

#define RDONLY_FLAGS    O_RDONLY
#define RDWR_FLAGS      O_RDWR
#define CROPTR_FLAGS    O_CREAT | O_TRUNC | O_RDWR

#if _FILE_OFFSET_BITS == 64

#undef RDONLY_FLAGS
#undef RDWR_FLAGS
#undef CROPTR_FLAGS

#define RDONLY_FLAGS    O_RDONLY | O_LARGEFILE
#define RDWR_FLAGS      O_RDWR | O_LARGEFILE
#define CROPTR_FLAGS    O_CREAT | O_TRUNC | O_RDWR | O_LARGEFILE

#endif

FILE *fopfile(const char *name,const char *mode)
{
FILE *fp;

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s\n", name, mode);
   exit(-1);
   }
return(fp);
}

int opfile_ro(const char *name)
{
int fd;
if ((fd = open (name, RDONLY_FLAGS, 0444)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

int opfile(const char *name)
{
int fd;
if ((fd = open (name, O_SYNC | RDWR_FLAGS, 0664)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

int croptrfile_sync(const char *name)
{
int fd;
if ((fd = open (name, O_SYNC | CROPTR_FLAGS, 0664)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

int croptrfile(const char *name)
{
int fd;
if ((fd = open (name, CROPTR_FLAGS, 0664)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

int reed(int fd, void *pntr, int length)
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

int rite(int fd, void *pntr, int length)
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

void *check_realloc(void *ptr,size_t len)
{
ptr = (char *) realloc (ptr,len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory reallocation error\n");
   exit(-1);
   }

return(ptr);
}

void *check_malloc(size_t len)
{
char *ptr;

ptr = (char *) malloc (len);
 
if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory allocation error\n");
   exit(-1);
   }
 
return(ptr);
}

void makedir(char *path)
{
char stmp[256], str[128];
int rtn, j;
mode_t mode = 00777;

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
      rtn = mkdir(stmp,mode);

         /*
      if(rtn == -1)
         {

fprintf(stderr,"**** errno=%d path=%s\n",errno,stmp);
fprintf(stderr,"\tEACCES=%d\n \tEEXIST=%d\n \tEFAULT=%d\n \tEIO=%d\n \tELOOP=%d\
n \tEMLINK=%d\n \tEMULTIHOP=%d\n \tENAMETOOLONG=%d\n \tENOENT=%d\n \tENOLINK=%d\
n \tENOSPC=%d\n \tENOTDIR=%d\n \tEROFS=%d\n\n",EACCES,
EEXIST,
EFAULT,
EIO,
ELOOP,
EMLINK,
EMULTIHOP,
ENAMETOOLONG,
ENOENT,
ENOLINK,
ENOSPC,
ENOTDIR,
EROFS);

         if(errno != EEXIST)
            {
            sprintf(str,"makedir() cannot make directory %s, exiting",stmp);
            perror(str);
            exit(-1);
            }
         }
         */
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
