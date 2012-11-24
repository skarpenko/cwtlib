/***********************************************************************
          Boundaries processing in Continuous Wavelet Transform

 Author  : Stepan V.Karpenko
 Date    : 19-04-2004
 Comments:
 History :
***********************************************************************/

#ifndef _FIXCWT_H_
#define _FIXCWT_H_

#ifdef __cplusplus
extern "C" {
#endif


/*
      Fix boundaries in CWT
      wt - cwt_t structure;
      s - source signal;
      ivalp - number of parts of interval between two samples, this
              parameter needs to increase discretization;
      psi - Wavelet function;

      Returns 0 on success and 1 on error.
*/
int fixcwt(cwt_t *wt, double *s, unsigned long ivalp, psi_t *psi);

#ifdef __cplusplus
}
#endif

#endif
