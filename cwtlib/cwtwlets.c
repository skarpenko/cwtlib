/***********************************************************************
           Wavelet Functions for Continuous Wavelet Transform

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
***********************************************************************/

#include <math.h>
#include <stdlib.h>
#include "cwtwlets.h"

/************************ PRIVATE ************************/

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
    { "MexHat", -5.0, 5.0, &MEXHATreal, NULL },
    { "Morlet", -4.0, 4.0, &MORLETreal, &MORLETimag },
    { "Gauss1", -5.0, 5.0, &GAUSS1real, NULL },
    { "Gauss2", -5.0, 5.0, &MEXHATreal, NULL },

    /* Last entry */
    { "", 0.0, 0.0, NULL, NULL }
};
