/***********************************************************************
                      Continuous Wavelet Transform

 Author  : Stepan V.Karpenko
 Date    : 01-04-2004
 Comments:
 History :
   31-05-2004 Stepan V.Karpenko
    Added two optimized versions of cwt().
***********************************************************************/

#ifndef _CWT_H_
#define _CWT_H_

#include "cwtwlets.h"


/* CWT structure */
typedef struct {
        char wname[WNAMELEN];      /* Wavelet name */
        double amin, astep, amax;  /* Minimum and maximum scale */
        double bstep;              /* Step of wavelet translation */
        unsigned long siglen;      /* Signal length */
        unsigned long rows, cols;  /* Dimension of cwt */
        double **cwt;              /* Wavelet coefficients */
} cwt_t;

#ifdef __cplusplus
extern "C" {
#endif

/*
      Perform continuous wavelet transform.
      s - source signal of length n;
      amin, astep, amax - Minimum scale, step and maximum scale;
      bstep - Wavelet translation step size;
      ivalp - used for increasing discretization in ivalp times, by
              inserting additional samples. This parameter affects
              precision;
      psi - Wavelet function;
      cwt - result;

      Returns 0 on success and 1 on error.
*/
int cwt( double *s, unsigned long n, double amin, double astep, double amax,
         double bstep, unsigned long ivalp, psi_t *psi, cwt_t *cwt );

/*
      Perform continuous wavelet transform. (optimized version 1)
      s - source signal of length n;
      amin, astep, amax - Minimum scale, step and maximum scale;
      bstep - Wavelet translation step size;
      ivalp - used for increasing discretization in ivalp times, by
              inserting additional samples. This parameter affects
              precision;
      wavelet - Wavelet;
      part - Wavelet part (REAL/IMAG);
      cwt - result;

      Returns 0 on success and 1 on error.
*/
int cwto1( double *s, unsigned long n, double amin, double astep, double amax,
           double bstep, unsigned long ivalp, cwtwlet_t *wavelet, long part,
           cwt_t *cwt );

/*
      Perform continuous wavelet transform. (optimized version 2)
      s - source signal of length n;
      amin, astep, amax - Minimum scale, step and maximum scale;
      bstep - Wavelet translation step size;
      ivalp - used for increasing discretization in ivalp times, by
              inserting additional samples. This parameter affects
              precision;
      wavelet - Wavelet;
      part - Wavelet part (REAL/IMAG);
      npoints - Points number for wavelet precompution
                (greater value - higher precision);
      cwt - result;

      Returns 0 on success and 1 on error.
*/
int cwto2( double *s, unsigned long n, double amin, double astep, double amax,
           double bstep, unsigned long ivalp, cwtwlet_t *wavelet, long part,
           unsigned long npoints, cwt_t *cwt );

/*
     Frees cwt structure allocated by cwt*() routines
*/
void free_cwt(cwt_t *cwt);

#ifdef __cplusplus
}
#endif

#endif
