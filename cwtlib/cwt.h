/***********************************************************************
                      Continuous Wavelet Transform

 Author  : Stepan V.Karpenko
 Date    : 01-04-2004
 Comments:
 History :
***********************************************************************/

#ifndef _CWT_H_
#define _CWT_H_

/* CWT structure */
typedef struct {
        double amin, astep, amax;    /* Minimum and maximum scale */
        double bstep;                /* Step of wavelet */
        unsigned siglen;             /* Length of signal */
        unsigned rows, cols;         /* Dimension of cwt */
        double **cwt;                /* Wavelet coefficients */
} cwt_t;

/* Wavelet function */
typedef double (psi_t)(double, double, double);

#ifdef __cplusplus
extern "C" {
#endif

/*
      Perform continuous wavelet transform.
      s - source signal of length n;
      amin, astep, amax - Minimum scale, step and maximum scale;
      bstep - Wavelet offset step size;
      ivalp - number of parts of interval between two samples, this
              parameter needs to increase discretization;
      psi - Wavelet function;
      cwt - result;

      Returns 0 on success and 1 on error.
*/
int cwt( double *s, unsigned n, double amin, double astep, double amax,
         double bstep, unsigned ivalp, psi_t *psi, cwt_t *cwt );

/*
     Frees cwt structure allocated by cwt()
*/
void free_cwt(cwt_t *cwt);

#ifdef __cplusplus
}
#endif

#endif
