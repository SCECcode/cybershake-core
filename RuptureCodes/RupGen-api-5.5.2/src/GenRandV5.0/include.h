#ifndef RUPGEN_INCLUDE_H
#define RUPGEN_INCLUDE_H

#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include <sys/file.h>
#include <sys/resource.h>
#include <sys/signal.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <sys/time.h>
#include <sys/types.h>

#include <string.h>
#include "fftw3.h"

#include "rupgen_api.h"
#include "rupgen_defs.h"

#ifdef _USE_MEMCACHED
#include <libmemcached/memcached.h>
#endif

#endif
