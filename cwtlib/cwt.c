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
   30-05-2004 Stepan V.Karpenko
    Minor bug fix in cwt(). ("i<n" replaced with "i<=n-1")
    Added vector allocation routines.
   31-05-2004 Stepan V.Karpenko
    Added two optimized versions of cwt().
***********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "cwtwlets.h"
#include "cwt.h"

#define NR_END 1
#define FREE_ARG char*


/************************ PRIVATE ************************/


/*
 * Memory allocation routines were taken from Numerical Recipes (www.nr.com)
 */

/* allocate a double vector with subscript range v[nl..nh] */
static double *Vectord(long nl, long nh)
{
        double *v;

        v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
        if (!v) return NULL;
        return v-nl+NR_END;
}

/* free a double vector allocated with Vectord() */
static void free_Vectord(double *v, long nl, long nh)
{
        free((FREE_ARG) (v+nl-NR_END));
}

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

/* free a double matrix allocated by Matrixd() */
static void free_Matrixd(double **m, long nrl, long nrh, long ncl, long nch)
{
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}


/************************ PUBLIC ************************/


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
            /* Translations */
            for ( dx = 0, b = 0.0; dx < cwt->cols; dx++, b+=bstep)
            {
                /* Perform convolution */
                cwt->cwt[dy][dx] = 0.0;
                for (i=0.0; i <= n-1; i+=istep)
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
        cwt->wname[0] = 0; /* Wavelet names is not supported in this version */

        return 0;
}

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
           cwt_t *cwt )
{
        double a, b, i, istep;
        unsigned long dx, dy;
        double t1, t2;
        psi_t *psi;

        if( (amin > amax) || (amin <= 0.0) || (astep <= 0.0) ||
            (amax <= 0.0) || (bstep <= 0.0) || !ivalp )
                return 1;

        /* Select part of wavelet (real/imaginary) */
        if(part == REAL) {
           psi = wavelet->real;
        } else {
           psi = wavelet->imag;
           if(psi == NULL) return 1;
        }

        /* initialize cwt dimension */
        cwt->cols = (unsigned long)( ceil((double)n / bstep) );
        cwt->rows = (unsigned long)((amax - amin)/astep + 1.0);

        cwt->cwt = Matrixd(0, cwt->rows-1, 0, cwt->cols-1);
        if(!cwt->cwt) return 1;

        istep = 1.0 / (double)ivalp;
        /* Scales */
        for (dy = 0, a = amin; dy < cwt->rows; dy++, a+=astep)
        {
            /* Translations */
            for ( dx = 0, b = 0.0; dx < cwt->cols; dx++, b+=bstep)
            {
                /* Compute wavelet boundaries */
                t1 = a*wavelet->esl + b;  t2 = a*wavelet->esr + b;
                if(t1<0.0) t1=0.0;        if(t2>=n) t2=(n-1);

                /* Perform convolution */
                cwt->cwt[dy][dx] = 0.0;
                for (i=t1; i <= t2; i+=istep)
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
        strncpy(cwt->wname, wavelet->wname, WNAMELEN); /* Store wavelet name */

        return 0;
}

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
           unsigned long npoints, cwt_t *cwt )
{
        double a, b, i, istep;
        unsigned long dx, dy;
        double t1, t2;
        double *W;
        double wstep;
        long L, R, j;
        double d, ind;
        psi_t *psi;

        if( (amin > amax) || (amin <= 0.0) || (astep <= 0.0) ||
            (amax <= 0.0) || (bstep <= 0.0) || !ivalp || !npoints )
                return 1;

        /* Select part of wavelet (real/imaginary) */
        if(part == REAL) {
           psi = wavelet->real;
        } else {
           psi = wavelet->imag;
           if(psi == NULL) return 1;
        }

        /* initialize cwt dimension */
        cwt->cols = (unsigned long)( ceil((double)n / bstep) );
        cwt->rows = (unsigned long)((amax - amin)/astep + 1.0);

        cwt->cwt = Matrixd(0, cwt->rows-1, 0, cwt->cols-1);
        if(!cwt->cwt) return 1;

        /* Precompute wavelet values */
        wstep = (double)( (wavelet->esr - wavelet->esl) / npoints );
        L = wavelet->esl / wstep;  R = wavelet->esr / wstep;
        W = Vectord(L, R);
        if(!W) {
            free_Matrixd(cwt->cwt, 0, cwt->rows-1, 0, cwt->cols-1);
            cwt->cwt = NULL;
            return 1;
        }
        for(j=L, i=wavelet->esl; j<=R; j++, i+=wstep)
             W[j] = psi(i, 1.0, 0.0);
        d = (double)( (R - L) / (wavelet->esr - wavelet->esl) );

        istep = 1.0 / (double)ivalp;
        /* Scales */
        for (dy = 0, a = amin; dy < cwt->rows; dy++, a+=astep)
        {
            /* Translations */
            for ( dx = 0, b = 0.0; dx < cwt->cols; dx++, b+=bstep)
            {
                /* Compute wavelet boundaries */
                t1 = a*wavelet->esl + b;  t2 = a*wavelet->esr + b;
                if(t1<0.0) t1=0.0;        if(t2>=n) t2=(n-1);

                /* Perform convolution */
                cwt->cwt[dy][dx] = 0.0;
                for (i=t1, j=0; i <= t2; i+=istep,j++) {
                    ind = d*(i-b)/a;
                    cwt->cwt[dy][dx] += s[(unsigned long)i] * W[(long)ind];
                }
                cwt->cwt[dy][dx] *= 1/(sqrt(a) * (double)ivalp);
            }
        }
        /* Free allocated memory */
        free_Vectord(W, L, R);

        /* store transform info */
        cwt->amin = amin;
        cwt->astep = astep;
        cwt->amax = a - astep; /* Last processed scale */
        cwt->bstep = bstep;
        cwt->siglen = n;
        strncpy(cwt->wname, wavelet->wname, WNAMELEN); /* Store wavelet name */

        return 0;
}

/*
     Frees cwt structure allocated by cwt*() routines
*/
void free_cwt(cwt_t *cwt)
{
        if(cwt->cwt)
            free_Matrixd(cwt->cwt, 0, cwt->rows-1, 0, cwt->cols-1);
        cwt->cwt = NULL;
        cwt->rows = cwt->cols = 0;
        cwt->amin = cwt->astep = cwt->amax = cwt->bstep = 0.0;
        cwt->siglen = 0;
        cwt->wname[0] = 0;
}
