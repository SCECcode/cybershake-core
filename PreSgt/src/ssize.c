#include <features.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>
#include <sys/times.h>

#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#include <sys/file.h>
#include <sys/resource.h>
#include <sys/signal.h>
#include <sys/syscall.h>

int main(int ac,char **av)
{
int ix, iy, iz;
long long indx;

while(scanf("%d %d %d",&ix,&iy,&iz) == 3)
   {
   indx = (long long)(ix)*(long long)(100000000) + (long long)(iy)*(long long)(10000) + (long long)(iz);

   fprintf(stderr,"%.4d %.4d %.4d\n",ix,iy,iz);
   fprintf(stderr,"%.12Ld\n",indx);
   }
}
