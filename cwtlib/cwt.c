/*
 *   cwt.c - Continuous Wavelet Transform Routines
 *
 *   Continuous Wavelet Transform Library
 *   Copyright (C) 2005 Stepan V.Karpenko <carp@mail.ru>
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the
 *   Free Software Foundation, Inc., 59 Temple Place, Suite 330,
 *   Boston, MA  02111-1307  USA
 */


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

/* compute fast Fourier transform */
static void fft(double *re, double *im, unsigned long n, int isign)
{
        unsigned long i, j, k, l, le, le1, ip, n2;
        double wpr, wpi, wr, wi, wtr, wti;

        n2 = n>>1;
        j = 1;
        for(i=0; i<n-1; i++) {
            if(i<j) {
               wtr     = re[j-1];
               wti     = im[j-1];
               re[j-1] = re[i];
               im[j-1] = im[i];
               re[i]   = wtr;
               im[i]   = wti;
            }
            k = n2;
            while(k<j) {
                j -= k;
                k >>= 1;
            }
            j += k;
        }
        l=1;
        k=n;
        while(k>>=1) {
            le1 = (le=1<<l++) >> 1;
            wtr = PI / (double)le1;
            wpr = cos(wtr); wpi = -isign*sin(wtr);
            wr = 1.0;       wi = 0.0;
            for(j=0; j<le1; j++) {
                for(i=j; i<n; i+=le) {
                    ip = i + le1;
                    wtr    = wr*re[ip] - wi*im[ip];
                    wti    = wi*re[ip] + wr*im[ip];
                    re[ip] = re[i] - wtr;
                    im[ip] = im[i] - wti;
                    re[i]  = re[i] + wtr;
                    im[i]  = im[i] + wti;
                }
                wr = (wtr=wr)*wpr - wi*wpi;
                wi = wi*wpr + wtr*wpi;
            }
        }
}


/************************ PUBLIC ************************/


/*
      Perform continuous wavelet transform. (Initial version)
*/
int cwt( double *s, unsigned long n, double amin, double astep, double amax,
         double bstep, unsigned long ivalp, cwtwlet_t *wavelet, long part,
         cwt_t *cwt )
{
        double a, b, i, istep;
        unsigned long dx, dy;
        psi_t *psi;
        double sign;

        if( (amin > amax) || (amin <= 0.0) || (astep <= 0.0) ||
            (amax <= 0.0) || (bstep <= 0.0) || !ivalp )
                return 1;

        /* Select part of wavelet (real/imaginary) */
        if(part == REAL) {
           psi = wavelet->real;
           sign = 1.0;
        } else {
           psi = wavelet->imag;
           sign = -1.0;
        }
        if(psi == NULL) return 1;

        /* initialize cwt dimensions */
        cwt->cols = (unsigned long)ceil( (double)n / bstep );
        cwt->rows = (unsigned long)floor( (amax - amin) / astep ) + 1;

        cwt->cwt = Matrixd(0, cwt->rows-1, 0, cwt->cols-1);
        if(!cwt->cwt) return 1;

        istep = 1.0 / (double)ivalp;
        /* Scales */
        for (dy = 0, a = amin; dy < cwt->rows; dy++, a+=astep)
        {
            /* Translations */
            for (dx = 0, b = 0.0; dx < cwt->cols; dx++, b+=bstep)
            {
                /* Perform convolution */
                cwt->cwt[dy][dx] = 0.0;
                for (i = 0.0; i <= n-1; i+=istep)
                    cwt->cwt[dy][dx] += s[(unsigned long)i] * psi(i, a, b);
                cwt->cwt[dy][dx] *= sign/(sqrt(a) * (double)ivalp);
            }
        }

        /* store transform info */
        cwt->amin = amin;
        cwt->astep = astep;
        cwt->amax = a - astep; /* Last processed scale */
        cwt->bstep = bstep;
        cwt->siglen = n;
        strncpy(cwt->wname, wavelet->wname, WNAMELEN); /* Store wavelet name */
        cwt->i = part; /* Complex part info */

        return 0;
}

/*
      Perform continuous wavelet transform. (optimized version 1)
*/
int cwto1( double *s, unsigned long n, double amin, double astep, double amax,
           double bstep, unsigned long ivalp, cwtwlet_t *wavelet, long part,
           cwt_t *cwt )
{
        double a, b, i, istep;
        unsigned long dx, dy;
        double t1, t2;
        psi_t *psi;
        double sign;

        if( (amin > amax) || (amin <= 0.0) || (astep <= 0.0) ||
            (amax <= 0.0) || (bstep <= 0.0) || !ivalp )
                return 1;

        /* Select part of wavelet (real/imaginary) */
        if(part == REAL) {
           psi = wavelet->real;
           sign = 1.0;
        } else {
           psi = wavelet->imag;
           sign = -1.0;
        }
        if(psi == NULL) return 1;

        /* initialize cwt dimensions */
        cwt->cols = (unsigned long)ceil( (double)n / bstep );
        cwt->rows = (unsigned long)floor( (amax - amin) / astep ) + 1;

        cwt->cwt = Matrixd(0, cwt->rows-1, 0, cwt->cols-1);
        if(!cwt->cwt) return 1;

        istep = 1.0 / (double)ivalp;
        /* Scales */
        for (dy = 0, a = amin; dy < cwt->rows; dy++, a+=astep)
        {
            /* Translations */
            for (dx = 0, b = 0.0; dx < cwt->cols; dx++, b+=bstep)
            {
                /* Compute wavelet boundaries */
                t1 = a*wavelet->esl + b;  t2 = a*wavelet->esr + b;
                if(t1<0.0) t1=0.0;        if(t2>=n) t2=(n-1);

                /* Perform convolution */
                cwt->cwt[dy][dx] = 0.0;
                for (i = t1; i <= t2; i+=istep)
                    cwt->cwt[dy][dx] += s[(unsigned long)i] * psi(i, a, b);
                cwt->cwt[dy][dx] *= sign/(sqrt(a) * (double)ivalp);
            }
        }

        /* store transform info */
        cwt->amin = amin;
        cwt->astep = astep;
        cwt->amax = a - astep; /* Last processed scale */
        cwt->bstep = bstep;
        cwt->siglen = n;
        strncpy(cwt->wname, wavelet->wname, WNAMELEN); /* Store wavelet name */
        cwt->i = part; /* Complex part info */

        return 0;
}

/*
      Perform continuous wavelet transform. (optimized version 2)
*/
int cwto2( double *s, unsigned long n, double amin, double astep, double amax,
           double bstep, unsigned long ivalp, cwtwlet_t *wavelet, long part,
           unsigned long npoints, cwt_t *cwt )
{
        /* variables for computing cwt */
        double a, b, i, istep;
        unsigned long dx, dy;
        psi_t *psi;
        double sign;
        /* variables for operating with precomputed wavelet */
        double t1, t2;
        double *W;
        double wstep;
        long L, R, j;
        double d, ind;

        if( (amin > amax) || (amin <= 0.0) || (astep <= 0.0) ||
            (amax <= 0.0) || (bstep <= 0.0) || !ivalp || !npoints )
                return 1;

        /* Select part of wavelet (real/imaginary) */
        if(part == REAL) {
           psi = wavelet->real;
           sign = 1.0;
        } else {
           psi = wavelet->imag;
           sign = -1.0;
        }
        if(psi == NULL) return 1;

        /* initialize cwt dimensions */
        cwt->cols = (unsigned long)ceil( (double)n / bstep );
        cwt->rows = (unsigned long)floor( (amax - amin) / astep ) + 1;

        cwt->cwt = Matrixd(0, cwt->rows-1, 0, cwt->cols-1);
        if(!cwt->cwt) return 1;

        /* Precompute wavelet values */
        wstep = (double)( (wavelet->esr - wavelet->esl) / npoints );
        L = (long)floor(wavelet->esl / wstep);  R = (long)ceil(wavelet->esr / wstep);
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
            for (dx = 0, b = 0.0; dx < cwt->cols; dx++, b+=bstep)
            {
                /* Compute wavelet boundaries */
                t1 = a*wavelet->esl + b;  t2 = a*wavelet->esr + b;
                if(t1<0.0) t1=0.0;        if(t2>=n) t2=(n-1);

                /* Perform convolution */
                cwt->cwt[dy][dx] = 0.0;
                for (i = t1; i <= t2; i+=istep) {
                    ind = d*(i-b)/a;
                    cwt->cwt[dy][dx] += s[(unsigned long)i] * W[(long)ind];
                }
                cwt->cwt[dy][dx] *= sign/(sqrt(a) * (double)ivalp);
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
        cwt->i = part; /* Complex part info */

        return 0;
}

/*
      Perform continuous wavelet transform. (optimized version 3)
*/
int cwto3( double *s, unsigned long n, double amin, double astep, double amax,
           double bstep, unsigned long ivalp, cwtwlet_t *wavelet, long part,
           unsigned long npoints, cwt_t *cwt )
{
        /* variables for computing cwt */
        double a, b, i, istep;
        unsigned long dx, dy;
        psi_t *psi;
        double sign;
        /* variables for operating with precomputed wavelet */
        double t1, t2;
        double *W;
        double wstep;
        long L, R, j;
        double d, ind;
        /* precomputed values */
        double ivlp_amin, a_esl, a_esr;


        if( (amin > amax) || (amin <= 0.0) || (astep <= 0.0) ||
            (amax <= 0.0) || (bstep <= 0.0) || !ivalp || !npoints )
                return 1;

        /* Select part of wavelet (real/imaginary) */
        if(part == REAL) {
           psi = wavelet->real;
           sign = 1.0;
        } else {
           psi = wavelet->imag;
           sign = -1.0;
        }
        if(psi == NULL) return 1;

        /* initialize cwt dimensions */
        cwt->cols = (unsigned long)ceil( (double)n / bstep );
        cwt->rows = (unsigned long)floor( (amax - amin) / astep ) + 1;

        cwt->cwt = Matrixd(0, cwt->rows-1, 0, cwt->cols-1);
        if(!cwt->cwt) return 1;

        /* Precompute wavelet values */
        wstep = (double)( (wavelet->esr - wavelet->esl) / npoints );
        L = (long)floor(wavelet->esl / wstep);  R = (long)ceil(wavelet->esr / wstep);
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
        ivlp_amin = (double)ivalp * amin;
        /* Scales */
        for (dy = 0, a = amin; dy < cwt->rows; dy++, a+=astep)
        {
            /* calculate ivalp for next a */
            if(dy && ivalp != 1) {
                ivalp = (unsigned long)ceil( ivlp_amin / a );
                istep = 1.0 / (double)ivalp;
            }
            /* small precompution */
            a_esl = a*wavelet->esl; a_esr = a*wavelet->esr;

            /* Translations */
            for (dx = 0, b = 0.0; dx < cwt->cols; dx++, b+=bstep)
            {
                /* Compute wavelet boundaries */
                t1 = a_esl + b;    t2 = a_esr + b;
                if(t1<0.0) t1=0.0; if(t2>=n) t2=(n-1);

                /* Perform convolution */
                cwt->cwt[dy][dx] = 0.0;
                for (i = t1; i <= t2; i+=istep) {
                    ind = d*(i-b)/a;
                    cwt->cwt[dy][dx] += s[(unsigned long)i] * W[(long)ind];
                }
                cwt->cwt[dy][dx] *= sign/(sqrt(a) * (double)ivalp);
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
        cwt->i = part; /* Complex part info */

        return 0;
}

/*
      Perform continuous wavelet transform. (FFT based version)
*/
int cwtft( double *s_re, double *s_im, unsigned long n, double amin,
           double astep, double amax, cwtwlet_t *wavelet,
           cwt_t *cwt_re, cwt_t *cwt_im )
{
        /* variables for computing cwt */
        double a;
        unsigned long dx, dy;
        unsigned long rows, cols;
        double *f_re, *f_im;
        double *r_re, *r_im;
        /* variables for wavelet computation */
        double w, W_re, W_im;
        psi_t *Psi_re, *Psi_im;
        /* precomputed values */
        double sqrt_a_n;
        double twoPIn = 2.0 * PI / (double)n;

        /* check that n is an integer power of two */
        dx=n; while(dx!=1) { if( dx&1 && dx>1 ) return 1; dx>>=1; }

        if( (amin > amax) || (amin <= 0.0) || (astep <= 0.0) ||
            (amax <= 0.0) || (!cwt_re && !cwt_im) )
                return 1;

        /* init pointers to wavelet in frequency domain */
        Psi_re = wavelet->realft;
        Psi_im = wavelet->imagft;
        if(!Psi_re && !Psi_im) return 1;

        /* allocate memory for source signal */
        f_re = Vectord(0, n-1);
        if(!f_re) return 1;
        f_im = Vectord(0, n-1);
        if(!f_im) {
            free_Vectord(f_re, 0, n-1);
            return 1;
        }

        /* initialize cwt dimensions */
        cols = n;
        rows = (unsigned long)floor( (amax - amin) / astep ) + 1;

        /* allocate memory for transform */
        r_re = NULL; r_im = NULL;
        if(cwt_re) {
            cwt_re->cols = cols;
            cwt_re->rows = rows;
            cwt_re->cwt = Matrixd(0, rows-1, 0, cols-1);
            if(cwt_re->cwt) r_re++; /* now r_re is not NULL */
         } else
            r_re = Vectord(0, cols-1);

         if(cwt_im) {
             cwt_im->cols = cols;
             cwt_im->rows = rows;
             cwt_im->cwt = Matrixd(0, rows-1, 0, cols-1);
             if(cwt_im->cwt) r_im++; /* now r_im is not NULL */
        } else
            r_im = Vectord(0, cols-1);

        /* check for memory allocation errors (r_re or r_im = NULL) */
        if(!r_re || !r_im) {
            free_Vectord(f_re, 0, n-1);
            free_Vectord(f_im, 0, n-1);

            if(cwt_re) {
               if(cwt_re->cwt) {
                   free_Matrixd(cwt_re->cwt, 0, rows-1, 0, cols-1);
                   cwt_re->cwt = NULL;
               }
            } else
               if(r_re) free_Vectord(r_re, 0, cols-1);

            if(cwt_im) {
               if(cwt_im->cwt) {
                   free_Matrixd(cwt_im->cwt, 0, rows-1, 0, cols-1);
                   cwt_im->cwt = NULL;
               }
            } else
               if(r_im) free_Vectord(r_im, 0, cols-1);

            return 1;
        }

        /* copy source signal */
        for (dx = 0; dx < n; dx++) {
            if(s_re) f_re[dx] = s_re[dx]; else f_re[dx] = 0.0;
            if(s_im) f_im[dx] = s_im[dx]; else f_im[dx] = 0.0;
        }

        /* forward Fourier transform */
        fft(f_re, f_im, n, 1);

        /* Scales */
        for (dy = 0, a = amin; dy < rows; dy++, a+=astep)
        {
            sqrt_a_n = sqrt(a) / (double)n; /* precompution */

            /* set pointers to transform result */
            if(cwt_re) r_re = cwt_re->cwt[dy];
            if(cwt_im) r_im = cwt_im->cwt[dy];

            /* Convolute */
            for (dx = 0; dx < cols; dx++)
            {
                /* calculate wave number w */
                if(dx<=n>>1)
                  w = twoPIn * (double)dx;
                else
                  w =-twoPIn * (double)(n-dx);

                /* calculate wavelet */
                if(!Psi_re) W_re = 0.0;
                else W_re = sqrt_a_n * Psi_re(w, a, 0.0);
                if(!Psi_im) W_im = 0.0;
                else W_im = sqrt_a_n * Psi_im(w, a, 0.0);

                /* compute result */
                r_re[dx] = f_re[dx] * W_re + f_im[dx] * W_im;
                r_im[dx] = f_im[dx] * W_re - f_re[dx] * W_im;
            }
            /* inverse Fourier transform */
            fft(r_re, r_im, n, -1);
        }

        /* store transform info */
        /* real part */
        if(cwt_re) {
            cwt_re->amin = amin;
            cwt_re->astep = astep;
            cwt_re->amax = a - astep; /* Last processed scale */
            cwt_re->bstep = 1.0;
            cwt_re->siglen = n;
            strncpy(cwt_re->wname, wavelet->wname, WNAMELEN);
            cwt_re->i = REAL;
        }
        /* imaginary part */
        if(cwt_im) {
            cwt_im->amin = amin;
            cwt_im->astep = astep;
            cwt_im->amax = a - astep;
            cwt_im->bstep = 1.0;
            cwt_im->siglen = n;
            strncpy(cwt_im->wname, wavelet->wname, WNAMELEN);
            cwt_im->i = IMAG;
        }

        /* Free allocated memory */
        free_Vectord(f_re, 0, n-1);
        free_Vectord(f_im, 0, n-1);
        if(!cwt_re)
            free_Vectord(r_re, 0, cols-1);
        if(!cwt_im)
            free_Vectord(r_im, 0, cols-1);

        return 0;
}

/*
    Copies one cwt structure to another
*/
int copy_cwt(cwt_t *src, cwt_t *dst)
{
        unsigned long dx,dy;

        /* allocate new area */
        dst->cwt = Matrixd(0, src->rows-1, 0, src->cols-1);
        if(!dst->cwt) return 1;

        /* copy cwt params */
        dst->cols   = src->cols;
        dst->rows   = src->rows;
        dst->amin   = src->amin;
        dst->astep  = src->astep;
        dst->amax   = src->amax;
        dst->bstep  = src->bstep;
        dst->siglen = src->siglen;
        dst->i      = src->i;
        strncpy(dst->wname, src->wname, WNAMELEN);

        /* copy cwt values */
        for (dy = 0; dy < dst->rows; dy++)
            for (dx = 0; dx < dst->cols; dx++)
                dst->cwt[dy][dx] = src->cwt[dy][dx];

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
        cwt->i = -1;
}
