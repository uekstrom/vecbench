#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#ifdef HAS_MKL
#include "mkl.h"
#endif

/* Vector size in MB */
#ifndef VECMB
#define VECMB 1
#endif

/* Number of doubles in the vector */
#define VECLEN ((1 << 17)*VECMB)

/* Number of calls to the benchmark function */
#ifndef NCALLS
#define NCALLS 10000
#endif



/* Fill the c array with doubles from (0,1) */
void fill_01(double c[], size_t len)
{
  size_t i;
  double x = M_PI;
  for (i=0;i<len;i++)
    {
      x += 0.1;
      c[i] = fmod(x,1);
    }
}

void copy(double dst[], const double src[], size_t len)
{
  size_t i;
  for (i=0;i<len;i++)
    dst[i] = src[i];
}

void copy_memcpy(double dst[], const double src[], size_t len)
{
  size_t i;
  memcpy(dst,src,len*sizeof(double));
}

#ifdef HAS_MKL
void copy_dcopy(double dst[], const double src[], size_t len)
{
  int one = 1, l = len;
  dcopy(&l, src, &one, dst, &one);
}
#endif

void poly1(double dst[], const double src[], size_t len)
{
  size_t i;
  for (i=0;i<len;i++)
    dst[i] = 1.2 + src[i]*2.3;
}

void poly2(double dst[], const double src[], size_t len)
{
  size_t i;
  for (i=0;i<len;i++)
    dst[i] = 1.2 + src[i]*(2.3 + 3.4*src[i]);
}

void rat22(double dst[], const double src[], size_t len)
{
  size_t i;
  for (i=0;i<len;i++)
    dst[i] = (1.2 + src[i]*(2.3 + 3.4*src[i]))/(5.2 + src[i]*(6.3 + 7.4*src[i]));
}

void logx(double dst[], const double src[], size_t len)
{
  size_t i;
  for (i=0;i<len;i++)
    dst[i] = log(src[i]);
}

void sqrtx(double dst[], const double src[], size_t len)
{
  size_t i;
  for (i=0;i<len;i++)
    dst[i] = sqrt(src[i]);
}

void asinhx(double dst[], const double src[], size_t len)
{
  size_t i;
  for (i=0;i<len;i++)
    dst[i] = asinh(src[i]);
}

void asinh_log(double dst[], const double src[], size_t len)
{
  size_t i;
  for (i=0;i<len;i++)
    dst[i] = asinh(sqrt(src[i]) + sqrt(1+src[i]));
}

void expx(double dst[], const double src[], size_t len)
{
  size_t i;
  for (i=0;i<len;i++)
    dst[i] = exp(src[i]);
}



typedef struct  __attribute__((aligned(32))) vecstruct
{
  double c[1<<30];
} vec_t;

void expx_vec_helper(vec_t * __restrict__ dst,
		     const vec_t * __restrict__ src,
		     size_t len)
{
  size_t i;
  for (i=0;i<len;i++)
    dst->c[i] = exp(src->c[i]);
}

void expx_vec(double dst[], const double src[], size_t len)
{
  expx_vec_helper((vec_t *)dst,(const vec_t *)src,len);
}

#define RUNFN(FN) \
  printf(#FN);\
  tick = clock();\
  for (i=0;i<NCALLS;i++)\
    {\
      FN(out,inp,VECLEN);\
      s += out[i % VECLEN]; \
    }\
  printf(" %.3f\n",(clock() - tick)/(double)CLOCKS_PER_SEC);

int main()
{
  int i;
  double s = 0;
  double *inp, *out;
  clock_t tick;
  printf("Input data size (MB) %i\n",(sizeof(*inp)*VECLEN)/(1 << 20));
  printf("Repetitions %i\n",NCALLS);
  inp = malloc(sizeof(*inp)*VECLEN);
  assert(inp && "Out of memory");
  out = malloc(sizeof(*out)*VECLEN);
  assert(out && "Out of memory");
  fill_01(inp,VECLEN);
  RUNFN(poly2);
  RUNFN(copy);
  RUNFN(copy_memcpy);
#ifdef HAS_MKL
  RUNFN(copy_dcopy);
#endif
  RUNFN(poly1);
  RUNFN(poly2);
  RUNFN(sqrtx);
  RUNFN(logx);
  RUNFN(expx);
  RUNFN(asinhx);
  RUNFN(asinh_log);
  printf("End of testing\n");
  free(inp);
  free(out);
  if (s == 1.2345667)
    return 1;
  else
    return 0;
}

