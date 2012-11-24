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
    Added cwtwlets[] array.
   30-05-2004 Stepan V.Karpenko
    Added entries for effective support boundaries.
   31-05-2004 Stepan V.Karpenko
    Added Gauss1 and Gauss2 (link to MexHat).
   13-10-2004 Stepan V.Karpenko
    Added Haar and French Hat.
***********************************************************************/

#include <math.h>
#include <stdlib.h>
#include "cwtwlets.h"

/************************ PRIVATE ************************/

/* Haar */
static double HAARreal(double x, double a, double b)
{
      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      if( x >= -1.0/2.0 && x <  0.0 )     return  1.0;
      if( x >=  0.0     && x <  1.0/2.0 ) return -1.0;
      else return 0.0;
}

/* French Hat */
static double FRHATreal(double x, double a, double b)
{
      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      if( fabs(x) <= 1.0/3.0 )                   return  1.0;
      if( fabs(x) >  1.0/3.0 && fabs(x) <= 1.0 ) return -1.0/2.0;
      else return 0.0;
}

/* Mexican Hat */
static double MEXHATreal(double x, double a, double b)
{
      const double c = 0.86732507058407771; /* 2 / (sqrt(3) * pi^(-1/4)) */

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x = x * x;
      return c * (1 - x) * exp(-x/2);
}

/* Morlet real part */
static double MORLETreal(double x, double a, double b)
{
      const double c = 2 * PI * 0.8;
      double x2;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      return exp(-x2/2) * cos(c*x);
}

/* Morlet imaginary part */
static double MORLETimag(double x, double a, double b)
{
      const double c = 2 * PI * 0.8;
      double x2;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      return exp(-x2/2) * sin(c*x);
}

/* Gaussian 1 */
static double GAUSS1real(double x, double a, double b)
{
      const double c = 1.7864876834760046; /* 2 * (2/pi)^(1/4) */
      double x2;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      return -c * x * exp(-x2);
}

/************************ PUBLIC ************************/

cwtwlet_t cwtwlets[] = {
    { "Haar",   -1.0, 1.0, &HAARreal,   NULL },
    { "FrHat",  -1.5, 1.5, &FRHATreal,  NULL },
    { "MexHat", -5.0, 5.0, &MEXHATreal, NULL },
    { "Morlet", -4.0, 4.0, &MORLETreal, &MORLETimag },
    { "Gauss1", -5.0, 5.0, &GAUSS1real, NULL },
    { "Gauss2", -5.0, 5.0, &MEXHATreal, NULL },

    /* Last entry */
    { "", 0.0, 0.0, NULL, NULL }
};
