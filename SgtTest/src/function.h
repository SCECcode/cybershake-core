void *check_malloc(int);
void *check_realloc(void *, int);
float *read_wccseis(char *, struct statdata *, float *, int);
void write_wccseis(char *, struct statdata *, float *, int);
FILE *fopfile(char*, char*);
int opfile_ro(char *);
int opfile(char *);
int croptrfile(char *);
int reed(int, void *, int);
int rite(int, void *, int);
void getheader(char *,struct statdata *);

void swap_in_place(int,char *);

