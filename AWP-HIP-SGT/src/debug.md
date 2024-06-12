# AWP DEBUG utilities
These notes describe a set of scripts that are helpful for both debugging and further
developing AWP. Before you use any of these scripts, please commit your work beforehand so
that you can easily revert the changes in case something goes wrong.

## Memory issues
We have found that certain bugs in AWP are due to memory errors. Tools such as `valgrind`
and `cuda-memcheck` are excellent for reporting many memory related issues. You are highly
encouraged to use these tools during development. So please use them to detect and fix
memory related issues before the code is used in production. As a precautionary measure,
we can avoid uninitialized memory errors by always initializing the allocated memory
immediately after it has been allocated.

### Uninitialized memory
The bash script `addmemset.sh` can be used to automatically add a statement that will
zero initialize memory after each device allocation call. **Warning:** Running this
script can modify all .c and .cu files in the source directory. Please commit any changes
and any untracked files before proceeding. No special instructions are required to run
this script, call
```bash
$ bash addmemset.sh
```

### Counting function calls
The number of calls to `malloc` should match the number of calls to `free`.  A simply way
to check that this is the case is to count number of matches from a grep search. For
convenience, the bash script `numcalls.sh` will count the number of occurrences for some
functions residing in either .c or .cu files. Example usage:
```bash
 bash numcalls.sh malloc free cudaMalloc cudaFree cudaMallocHost cudaFreeHost cudaMemset
malloc: 19
free: 19
cudaMalloc: 173
cudaFree: 188
cudaMallocHost: 40
cudaFreeHost: 40
cudaMemset: 204
```

## MPI and CUDA errors
In many cases it is a lot easier to identify problems if errors are reported as soon as
they occur. However, in practice it can be quite tedious to check all error codes. When
it comes to cuda and MPI error codes, it is quite easy to check them because each CUDA
function starts with `cuda...` and returns an error code. The same holds true for MPI
functions. With the help of the preprocessor macros `CUCHK` and MPICHK` we can intercept
any error from these libraries and report the error message, as well as file, function,
and line number at which it occured. These macros are defined in `pmcl3d_cons.h`. Both of
these macros can be disabled by defining `NDEBUG`. To avoid having to modify the source,
you can instead modify the Makefile by adding the following to it,
```bash
CFLAGS = ... -DNDEBUG
GFLAGS = ... -DNDEBUG
```
The dots `...` represent any already existing arguments.


### Automatic checking
The
makefile `Makefile.debug` can modify the source to enable or disable these additional
checks. Both commands will modify source files, but the disabling call listed below should return to the
source files to their original state. To be safe, make sure to commit any changes
(including untracked files) before proceeding.

To enable, call
```bash
$ make -f Makefile.debug enable
```
To disable, call
```bash
$ make -f Makefile.debug disable
```

### Manual checking
You can also add these checks manually. Here is an example that demonstrates how to check
if a call to cudaMalloc is successful. We add the following statements to `pmcl3d.c`
```
float *var;
CUCHK(cudaMalloc((void**)&var, -1));
```

In this example, the number of bytes is invalid and will therefore generate an error. The
error message is:
```
CUDA error in pmcl3d.c:776 main(): invalid argument.
```

Since MPI exits by default on any error, it is necessary to first enable error handling.
After `MPI_Init(...)` in `main()`, add the statement
```
MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
```
After doing that, you can check for MPI errors by wrapping the MPI function call of
interest in `MPICHK` (similar to `CUCHK`).





