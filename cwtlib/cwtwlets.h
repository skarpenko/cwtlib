/*
 *   Continuous Wavelet Transform Library
 *   Copyright (C) 2004 Stepan V.Karpenko
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
 Title   : Wavelet Functions for Continuous Wavelet Transform
 Author  : Stepan V.Karpenko
 Date    : 01-04-2004
 Comments:
 History :
   09-04-2004 Stepan V.Karpenko
    Added cwtwlets[] array and cwtwlet_t type.
   30-05-2004 Stepan V.Karpenko
    Added entries for effective support boundaries.
   31-05-2004 Stepan V.Karpenko
    Added Gauss1 and Gauss2 (link to MexHat).
   13-10-2004 Stepan V.Karpenko
    Added Haar and French Hat.
***********************************************************************/

#ifndef _CWTWLETS_H_
#define _CWTWLETS_H_


/* Definitions */
#ifndef PI
  #define PI  3.14159265358979323846
#endif
#define TINY 1E-20
#define WNAMELEN 10

/* Wavelet function */
typedef double (psi_t)(double, double, double);

/* Wavelet structure */
typedef struct {
        char wname[WNAMELEN];  /* Wavelet name */
        double esl, esr;       /* Effective support boundaries */
        psi_t *real;           /* Real part */
        psi_t *imag;           /* Imaginary part, if not exist must be NULL */
} cwtwlet_t;

#ifdef __cplusplus
extern "C" {
#endif

/* Possible parts of wavelet */
enum { REAL, IMAG };

/* Enumeration of wavelet entries */
enum { HAAR, FRHAT, MEXHAT, MORLET, GAUSS1, GAUSS2 };
extern cwtwlet_t cwtwlets[]; /* Array of wavelets */

#ifdef __cplusplus
}
#endif

#endif
