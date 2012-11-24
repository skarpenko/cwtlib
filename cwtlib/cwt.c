/***********************************************************************
                      Continuous Wavelet Transform

 Author  : Stepan V.Karpenko
 Date    : 01-04-2004
 Comments:
 History :
   07-04-2004 Stepan V.Karpenko
    Added minor precautions against possible bugs and incompatibility
    in compilers.
    Now actual maximum scale stores in cwt structure.
***********************************************************************/

#include <math.h>
#include <stdlib.h>
#include "cwtwlets.h"
#include "cwt.h"

#define NR_END 1
#define FREE_ARG char*

/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
static double **Matrixd(long nrl, long nrh, long ncl, long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        double **m;

        /* allocate pointers to rows */
        m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
        if (!m) return NULL;
        m += NR_END;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
        if (!m[nrl]) { free((FREE_ARG) (m+nrl-NR_END)); return NULL; }
        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

/* free a double matrix allocated by dmatrix() */
static void free_Matrixd(double **m, long nrl, long nrh, long ncl, long nch)
{
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

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
int cwt( double *s, unsigned long n, double amin, double astep, double amax,
         double bstep, unsigned long ivalp, psi_t *psi, cwt_t *cwt )
{
        double a, b, i, istep;
        unsigned long dx, dy;

        if( (amin > amax) || (amin <= 0.0) || (astep <= 0.0) ||
            (amax <= 0.0) || (bstep <= 0.0) || !ivalp )
                return 1;

        /* initialize cwt dimension */
        cwt->cols = (unsigned long)( ceil((double)n / bstep) );
        cwt->rows = (unsigned long)((amax - amin)/astep + 1.0);

        cwt->cwt = Matrixd(0, cwt->rows-1, 0, cwt->cols-1);
        if(!cwt->cwt) return 1;

        istep = 1.0 / (double)ivalp;
        /* Scales */
        for (dy = 0, a = amin; dy < cwt->rows; dy++, a+=astep)
        {
           /* Offsets */
            for ( dx = 0, b = 0.0; dx < cwt->cols; dx++, b+=bstep)
            {
               /* Perform convolution */
                cwt->cwt[dy][dx] = 0.0;
                for (i=0.0; i < n; i+=istep)
                    cwt->cwt[dy][dx] += s[(unsigned long)i] * psi(i, a, b);
                cwt->cwt[dy][dx] *= 1/(sqrt(a) * (double)ivalp);
            }
        }

        /* store transform info */
        cwt->amin = amin;
        cwt->astep = astep;
        cwt->amax = a - astep; /* Last processed scale */
        cwt->bstep = bstep;
        cwt->siglen = n;

        return 0;
}

/*
     Frees cwt structure allocated by cwt()
*/
void free_cwt(cwt_t *cwt)
{
        free_Matrixd(cwt->cwt, 0, cwt->rows-1, 0, cwt->cols-1);
        cwt->cwt = NULL;
        cwt->rows = cwt->cols = 0;
        cwt->amin = cwt->astep = cwt->amax = cwt->bstep = 0.0;
        cwt->siglen = 0;
}
