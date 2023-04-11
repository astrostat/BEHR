
#define size_t long

void **set_matrix(long nrow, long ncol, size_t size );
void f_matrix( void **m, long nrow, long ncol, size_t size);
double **matrix(long nrl, long nrh, long ncl, long nch);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_matrix2(double **m, long nrow);

void advnst(long k);
float genbet(float aa,float bb);
float genchi(float df);
float genexp(float av);
float genf(float dfn, float dfd);
float gengam(float a,float r);
void genmn(float *parm,float *x,float *work);
void genmul(long n,float *p,long ncat,long *ix);
float gennch(float df,float xnonc);
float gennf(float dfn, float dfd, float xnonc);
float gennor(float av,float sd);
void genprm(long *iarray,int larray);
float genunf(float low,float high);
void getsd(long *iseed1,long *iseed2);
void gscgn(long getset,long *g);
long ignbin(long n,float pp);
long ignnbn(long n,float p);
long ignlgi(void);
long ignpoi(float mu);
long ignuin(long low,long high);
void initgn(long isdtyp);
long mltmod(long a,long s,long m);
void phrtsd(char* phrase,long* seed1,long* seed2);
float ranf(void);
void setall(long iseed1,long iseed2);
void setant(long qvalue);
void setgmn(float *meanv,float *covm,long p,float *parm);
void setsd(long iseed1,long iseed2);
float sexpo(void);
float sgamma(float a);
float snorm(void);

int sameString(char *str1, char *str2, int length);
int *turning(double *dens, int length);
void ppoints(int n, double left, double right, double *out);
void segment(int n, double left, double right, double *out);
int compare_doubles(const void *a, const void *b);
void split_fn(double *array);
void HPD_int(double *dens, double curr, double degree, double level, int LOG, double *temp);
void summary_fn(double *array, int length, double level, int HPD, double *results);
void HPD_fn(double *dens, double *points, double level, int LOG, double *results);
double maximum(double left, double right, double src, double bkg, double ratio, double eff, double *prior);
double findprob(double lam, int nrow, double **priormat);
double post_(double *param);
double dens_S(double param);
double dens_H(double param);

