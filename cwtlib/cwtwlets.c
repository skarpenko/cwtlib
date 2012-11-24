/***********************************************************************
           Wavelet Functions for Continuous Wavelet Transform

 Author  : Stepan V.Karpenko
 Date    : 01-04-2004
 Comments:
 History :
   09-04-2004 Stepan V.Karpenko
    Added cwtwlets[] array.
***********************************************************************/

#include <math.h>
#include <stdlib.h>
#include "cwtwlets.h"

/************************ PRIVATE ************************/

/* Mexican Hat */
static double MEXHATreal(double x, double a, double b)
{
      const double c = 0.86732507058407771; /* 2 / sqrt(3) * pi^(-1/4) */

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

/************************ PUBLIC ************************/

cwtwlet_t cwtwlets[] = {
    { "MexHat", &MEXHATreal, NULL },
    { "Morlet", &MORLETreal, &MORLETimag },

    /* Last entry */
    { "", NULL, NULL }
};
