/*
 *   Continuous Wavelet Transform Library
 *   Copyright (C) 2005 Stepan V.Karpenko
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

/***********************************************************************
 Title   : Continuous Wavelet Transform Routines
 Author  : Stepan V.Karpenko
 Date    : 01-04-2004
 Comments:
 History :
   31-05-2004 Stepan V.Karpenko
    Added two optimized versions of cwt().
   08-10-2004 Stepan V.Karpenko
    Added cwto3() and copy_cwt().
   11-08-2005 Stepan V.Karpenko
    New interface of initial version cwt() routine.
   14-08-2005 Stepan V.Karpenko
    "i" field inserted to cwt_t, that indicates complex part of WT.
***********************************************************************/

#ifndef _CWT_H_
#define _CWT_H_

#include "cwtwlets.h"


/* CWT structure */
typedef struct {
        char wname[WNAMELEN];      /* Wavelet name */
        long i;                    /* Complex part (REAL/IMAG) */
        double amin, astep, amax;  /* Minimum and maximum scale */
        double bstep;              /* Step of wavelet translation */
        unsigned long siglen;      /* Signal length */
        unsigned long rows, cols;  /* Dimensions of cwt */
        double **cwt;              /* Wavelet coefficients */
} cwt_t;

#ifdef __cplusplus
extern "C" {
#endif

/*
      Perform continuous wavelet transform. (Initial version)
      s - source signal of length n;
      amin, astep, amax - Minimum scale, step and maximum scale;
      bstep - Wavelet translation step size;
      ivalp - used for increasing discretization in ivalp times, by
              inserting additional samples. This parameter affects
              precision;
      wavelet - Wavelet;
      part - Complex part of wavelet (REAL/IMAG);
      cwt - result;

      Returns 0 on success and 1 on error.
*/
int cwt( double *s, unsigned long n, double amin, double astep, double amax,
         double bstep, unsigned long ivalp, cwtwlet_t *wavelet, long part,
         cwt_t *cwt );

/*
      Perform continuous wavelet transform. (optimized version 1)
      s - source signal of length n;
      amin, astep, amax - Minimum scale, step and maximum scale;
      bstep - Wavelet translation step size;
      ivalp - used for increasing discretization in ivalp times, by
              inserting additional samples. This parameter affects
              precision;
      wavelet - Wavelet;
      part - Complex part of wavelet (REAL/IMAG);
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
      part - Complex part of wavelet (REAL/IMAG);
      npoints - Points number for wavelet precompution
                (greater value - higher precision);
      cwt - result;

      Returns 0 on success and 1 on error.
*/
int cwto2( double *s, unsigned long n, double amin, double astep, double amax,
           double bstep, unsigned long ivalp, cwtwlet_t *wavelet, long part,
           unsigned long npoints, cwt_t *cwt );

/*
      Perform continuous wavelet transform. (optimized version 3)
      s - source signal of length n;
      amin, astep, amax - Minimum scale, step and maximum scale;
      bstep - Wavelet translation step size;
      ivalp - used for increasing discretization in ivalp times, by
              inserting additional samples. This parameter affects
              precision;
      wavelet - Wavelet;
      part - Complex part of wavelet (REAL/IMAG);
      npoints - Points number for wavelet precompution
                (greater value - higher precision);
      cwt - result;

      Returns 0 on success and 1 on error.
*/
int cwto3( double *s, unsigned long n, double amin, double astep, double amax,
           double bstep, unsigned long ivalp, cwtwlet_t *wavelet, long part,
           unsigned long npoints, cwt_t *cwt );

/*
    Copies one cwt structure to another
    (NOTE: dst must be empty or memory leak will occur)

    Returns 0 on success and 1 on error.
*/
int copy_cwt(cwt_t *src, cwt_t *dst);

/*
     Frees cwt structure allocated by cwt*() routines
*/
void free_cwt(cwt_t *cwt);

#ifdef __cplusplus
}
#endif

#endif
