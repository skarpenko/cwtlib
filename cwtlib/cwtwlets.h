/*
 *   cwtwlets.h - Wavelets for Continuous Wavelet Transform
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


#ifndef _CWTWLETS_H_
#define _CWTWLETS_H_


/* Definitions */
#ifndef PI
  #define PI  3.14159265358979323846
#endif
#define TINY 1E-200
#define WNAMELEN 10

/* Wavelet function */
typedef double (psi_t)(double, double, double);

/* Wavelet structure */
typedef struct {
        char wname[WNAMELEN];  /* Wavelet name */
        double esl, esr;       /* Effective support boundaries */
        psi_t *real;           /* Real part in time domain */
        psi_t *imag;           /* Imaginary part in time domain */
        psi_t *realft;         /* Real part in frequency domain */
        psi_t *imagft;         /* Imaginary part in frequency domain */
} cwtwlet_t;

#ifdef __cplusplus
extern "C" {
#endif

/* Complex parts */
enum { REAL, IMAG };

/* Enumeration of wavelet entries */
enum {
   HAAR,    FRHAT,   SPLASH,  MEXHAT,   POISSON,
   MORLET,  CMORLET, SHANNON, ESHANNON, WAVE,
   GAUSS1,  GAUSS2,  GAUSS3,  GAUSS4,   GAUSS5,
   GAUSS6,  GAUSS7,  GAUSS8,  GAUSS9,   GAUSS10,
   CGAUSS1, CGAUSS2, CGAUSS3, CGAUSS4,  CGAUSS5,
   CGAUSS6, CGAUSS7, CGAUSS8, CGAUSS9,  CGAUSS10,
   PAUL1,   PAUL2,   PAUL3,   PAUL4,    PAUL5
};

extern cwtwlet_t cwtwlets[]; /* Array of wavelets */

#ifdef __cplusplus
}
#endif

#endif
