/***********************************************************************
          Boundaries processing in Continuous Wavelet Transform

 Author  : Stepan V.Karpenko
 Date    : 19-04-2004
 Comments:
 History :
***********************************************************************/

#include <stdlib.h>
#include <math.h>
#include "cwt.h"
#include "cwtwlets.h"

/*
      Fix boundaries in CWT
      wt - cwt_t structure;
      s - source signal;
      ivalp - number of parts of interval between two samples, this
              parameter needs to increase discretization;
      psi - Wavelet function;

      Returns 0 on success and 1 on error.
*/
int fixcwt(cwt_t *wt, double *s, unsigned long ivalp, psi_t *psi)
{
     double avg = 0.0;
     unsigned long i,j;
     unsigned long n = wt->siglen;
     double *l;
     cwt_t l_wt;

     /* find average value */
     for(i=0; i<n; i++) avg += s[i]; avg/=(double)n;

     /* create series of average values */
     l = (double *)malloc(n*sizeof(double));
     if(!l) return 1;
     for(i=0; i<n; i++) l[i] = avg;

     /* perform continuous wavelet transform */
     if(cwt(l, n, wt->amin, wt->astep, wt->amax, wt->bstep, ivalp, psi, &l_wt))
          return 1;

     /* subtract decompositions */
     for(i=0; i < wt->rows; i++)
         for(j=0; j < wt->cols; j++)
               wt->cwt[i][j] -= l_wt.cwt[i][j];

     /* free allocated structures */
     free(l);
     free_cwt(&l_wt);

     return 0;
}
