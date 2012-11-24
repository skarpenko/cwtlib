/***********************************************************************
           Wavelet Functions for Continuous Wavelet Transform

 Author  : Stepan V.Karpenko
 Date    : 01-04-2004
 Comments:
 History :
***********************************************************************/

#include <math.h>
#include "cwtwlets.h"

/* Mexican Hat */
double MEXHAT(double x, double a, double b)
{
      const double c = 0.86732507058407771; /* 2 / sqrt(3) * pi^(-1/4) */

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x = x * x;
      return c * (1 - x) * exp(-x/2);
}

/* Morlet real part */
double MORLETreal(double x, double a, double b)
{
      const double c = 2 * PI * 0.8;
      double x2;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      return exp(-x2/2) * cos(c*x);
}

/* Morlet imaginary part */
double MORLETimag(double x, double a, double b)
{
      const double c = 2 * PI * 0.8;
      double x2;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      return exp(-x2/2) * sin(c*x);
}
