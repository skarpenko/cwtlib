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
   11-08-2005 Stepan V.Karpenko
    Added Wave (link to Gauss1).
   12-08-2005 Stepan V.Karpenko
    Gauss1 is a little bit modified (x2 replaced with x2/2).
    Added Gauss2 - Gauss10. Gauss2 now is not a link to MexHat.
   14-08-2005 Stepan V.Karpenko
    Added complex Gauss1 - Gauss10.
    Added complex Morlet (CMORLET). MORLET is now real-valued wavelet.
   15-08-2005 Stepan V.Karpenko
    Added complex Shannon wavelet and real-valued Poisson wavelet.
   16-08-2005 Stepan V.Karpenko
    Added Paul1 - Paul5 wavelets.
   22-08-2005 Stepan V.Karpenko
    Added real-valued Splash wavelet and complex Shannon wavelet with
    exponential decay.
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

/* Splash */
static double SPLASHreal(double x, double a, double b)
{
      const double Fb = 0.2; /* bandwidth parameter */
                /* c = sqrt(2/Fb^3) */
      const double c = sqrt( 2.0/(Fb*Fb*Fb) );

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      return c * x * exp( -fabs(x)/Fb );
}

/* Mexican Hat */
static double MEXHATreal(double x, double a, double b)
{
                /* c = 2 / ( sqrt(3) * pi^(1/4) ) */
      const double c = 0.8673250705840776;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x = x * x;
      return c * (1.0 - x) * exp(-x/2);
}

/* Poisson */
static double POISSONreal(double x, double a, double b)
{
                /* c = 2.0 / sqrt(pi) */
      const double c = 1.128379167095513;
      double x2;
      double dnm;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      dnm = 1.0 + x2;
      return c * (1.0 - x2) / ( dnm*dnm );
}

/* Morlet */
static double MORLETreal(double x, double a, double b)
{
      const double Fb = 2.0; /* bandwidth parameter */
      const double Fc = 0.8; /* wavelet center frequency */
      double x2;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      return exp(-x2/Fb) * cos(2*PI*Fc*x);
}

/* Complex Morlet. Real part */
static double CMORLETreal(double x, double a, double b)
{
      const double Fb = 2.0; /* bandwidth parameter */
      const double Fc = 0.8; /* wavelet center frequency */
      const double c = 1.0 / sqrt(PI*Fb);
      double x2;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      return c * exp(-x2/Fb) * cos(2*PI*Fc*x);
}

/* Complex Morlet. Imaginary part */
static double CMORLETimag(double x, double a, double b)
{
      const double Fb = 2.0; /* bandwidth parameter */
      const double Fc = 0.8; /* wavelet center frequency */
      const double c = 1.0 / sqrt(PI*Fb);
      double x2;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      return c * exp(-x2/Fb) * sin(2*PI*Fc*x);
}

/* Complex Shannon wavelet. Real part */
static double SHANNONreal(double x, double a, double b)
{
                /* c = 1 / Integral( f(x)^2, x=-Inf..+Inf )^(1/2) */
      const double c = 2.2360679774997897;
      const double Fb = 2.0; /* bandwidth parameter */
      const double Fc = 0.8; /* wavelet center frequency */

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      return c * sin(PI*Fb*x) * cos(2*PI*Fc*x) / ( sqrt(Fb) * PI * x );
}

/* Complex Shannon wavelet. Imaginary part */
static double SHANNONimag(double x, double a, double b)
{
                /* c = 1 / Integral( f(x)^2, x=-Inf..+Inf )^(1/2) */
      const double c = 2.2360679774997897;
      const double Fb = 2.0; /* bandwidth parameter */
      const double Fc = 0.8; /* wavelet center frequency */

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      return c * sin(PI*Fb*x) * sin(2*PI*Fc*x) / ( sqrt(Fb) * PI * x );
}

/* Complex Shannon wavelet with exponential decay. Real part */
static double ESHANNONreal(double x, double a, double b)
{
                /* c = 1 / Integral( f(x)^2, x=-Inf..+Inf )^(1/2) */
      const double c = 5.9605667650728377;
      const double Fb = 1.0; /* bandwidth parameter */
      const double Fc = 1.0; /* wavelet center frequency */

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      return c * sin(PI*Fb*x) * cos(2*PI*Fc*x) /
             ( sqrt(Fb) * PI * x * exp(Fb*abs(x)) );
}

/* Complex Shannon wavelet with exponential decay. Imaginary part */
static double ESHANNONimag(double x, double a, double b)
{
                /* c = 1 / Integral( f(x)^2, x=-Inf..+Inf )^(1/2) */
      const double c = 5.9605667650728377;
      const double Fb = 1.0; /* bandwidth parameter */
      const double Fc = 1.0; /* wavelet center frequency */

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      return c * sin(PI*Fb*x) * sin(2*PI*Fc*x) /
             ( sqrt(Fb) * PI * x * exp(Fb*abs(x)) );
}

/* Gaussian 1 */
static double GAUSS1real(double x, double a, double b)
{
                /* c = sqrt(2) / pi^(1/4) */
      const double c = 1.062251932027197;
      double x2;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      return -c * x * exp(-x2/2);
}

/* Gaussian 2 */
static double GAUSS2real(double x, double a, double b)
{
                /* c = 2 / ( sqrt(3) * pi^(1/4) ) */
      const double c = 0.8673250705840776;
      double x2;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      return c * exp(-x2/2) * (x2 - 1.0);
}

/* Gaussian 3 */
static double GAUSS3real(double x, double a, double b)
{
                /* c = 2 * sqrt(2) / ( sqrt(15) * pi^(1/4) ) */
      const double c = 0.5485445389623982;
      double x2;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      return -c * x * exp(-x2/2) * (x2 - 3.0);
}

/* Gaussian 4 */
static double GAUSS4real(double x, double a, double b)
{
                /* c = 4 / ( sqrt(105) * pi^(1/4) ) */
      const double c = 0.2932093894547376;
      double x2, x4;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      return c * exp(-x2/2) * (x4 - 6.0*x2 + 3.0);
}

/* Gaussian 5 */
static double GAUSS5real(double x, double a, double b)
{
                /* c = 4 * sqrt(2) / ( sqrt(945) * pi^(1/4) ) */
      const double c = 0.1382202317273416;
      double x2, x4;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      return -c * x * exp(-x2/2) * (x4 - 10.0*x2 + 15.0);
}

/* Gaussian 6 */
static double GAUSS6real(double x, double a, double b)
{
                /* c = 8 / ( sqrt(10395) * pi^(1/4) ) */
      const double c = 0.05893730483821539;
      double x2, x4, x6;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      x6 = x4 * x2;
      return c * exp(-x2/2) * (x6 - 15.0*x4 + 45.0*x2 - 15.0);
}

/* Gaussian 7 */
static double GAUSS7real(double x, double a, double b)
{
                /* c = 8 * sqrt(2) / ( sqrt(135135) * pi^(1/4) ) */
      const double c = 0.02311711288066360;
      double x2, x4, x6;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      x6 = x4 * x2;
      return -c * x * exp(-x2/2) * (x6 - 21.0*x4 + 105.0*x2 - 105.0);
}

/* Gaussian 8 */
static double GAUSS8real(double x, double a, double b)
{
                /* c = 16 / ( sqrt(2027025) * pi^(1/4) ) */
      const double c = 0.008441176126088454;
      double x2, x4, x6, x8;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      x6 = x4 * x2;
      x8 = x4 * x4;
      return c * exp(-x2/2) * (x8 - 28.0*x6 + 210.0*x4 - 420.0*x2 + 105.0);
}

/* Gaussian 9 */
static double GAUSS9real(double x, double a, double b)
{
                /* c = 16 * sqrt(2) / ( sqrt(34459425) * pi^(1/4) ) */
      const double c = 0.002895299525125788;
      double x2, x4, x6, x8;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      x6 = x4 * x2;
      x8 = x4 * x4;
      return -c * x * exp(-x2/2) * (x8 - 36.0*x6 + 378.0*x4 - 1260.0*x2 + 945.0);
}

/* Gaussian 10 */
static double GAUSS10real(double x, double a, double b)
{
                /* c = 32 / ( sqrt(654729075) * pi^(1/4) ) */
      const double c = 0.0009393592071302543;
      double x2, x4, x6, x8, x10;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      x6 = x4 * x2;
      x8 = x4 * x4;
      x10 = x8 * x2;
      return c * exp(-x2/2) * (x10 - 45.0*x8 + 630.0*x6 - 3150.0*x4 + 4725.0*x2 - 945.0);
}

/* Complex Gaussian 1. Real part */
static double CGAUSS1real(double x, double a, double b)
{
                /* c = sqrt(2) * exp(1/8) / pi^(1/4) */
      const double c = 1.203689133543866;
      double x2;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( -ce*x - (1.0/2.0)*se );
}

/* Complex Gaussian 1. Imaginary part */
static double CGAUSS1imag(double x, double a, double b)
{
                /* c = sqrt(2) * exp(1/8) / pi^(1/4) */
      const double c = 1.203689133543866;
      double x2;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( se*x - (1.0/2.0)*ce );
}

/* Complex Gaussian 2. Real part */
static double CGAUSS2real(double x, double a, double b)
{
                /* c = 2 * exp(1/8) / ( sqrt(3) * pi^(1/4) ) */
      const double c = 0.9828080620384235;
      double x2;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( ce*x2 + se*x - (5.0/4.0)*ce );
}

/* Complex Gaussian 2. Imaginary part */
static double CGAUSS2imag(double x, double a, double b)
{
                /* c = 2 * exp(1/8) / ( sqrt(3) * pi^(1/4) ) */
      const double c = 0.9828080620384235;
      double x2;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( -se*x2 + ce*x + (5.0/4.0)*se );
}

/* Complex Gaussian 3. Real part */
static double CGAUSS3real(double x, double a, double b)
{
                /* c = 2 * sqrt(2) * exp(1/8) / ( sqrt(15) * pi^(1/4) ) */
      const double c = 0.6215823957634969;
      double x2, x3;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x3 = x2 * x;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( -ce*x3 - (3.0/2.0)*se*x2 + (15.0/4.0)*ce*x + (13.0/8.0)*se );
}

/* Complex Gaussian 3. Imaginary part */
static double CGAUSS3imag(double x, double a, double b)
{
                /* c = 2 * sqrt(2) * exp(1/8) / ( sqrt(15) * pi^(1/4) ) */
      const double c = 0.6215823957634969;
      double x2, x3;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x3 = x2 * x;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( se*x3 - (3.0/2.0)*ce*x2 - (15.0/4.0)*se*x + (13.0/8.0)*ce );
}

/* Complex Gaussian 4. Real part */
static double CGAUSS4real(double x, double a, double b)
{
                /* c = 4 * exp(1/8) / ( sqrt(105) * pi^(1/4) ) */
      const double c = 0.3322497660853046;
      double x2, x3, x4;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x3 = x2 * x;
      x4 = x2 * x2;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( ce*x4 + 2.0*se*x3 - (15.0/2.0)*ce*x2 - (13.0/2.0)*se*x
                  + (73.0/16.0)*ce );
}

/* Complex Gaussian 4. Imaginary part */
static double CGAUSS4imag(double x, double a, double b)
{
                /* c = 4 * exp(1/8) / ( sqrt(105) * pi^(1/4) ) */
      const double c = 0.3322497660853046;
      double x2, x3, x4;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x3 = x2 * x;
      x4 = x2 * x2;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( -se*x4 + 2.0*ce*x3 + (15.0/2.0)*se*x2 - (13.0/2.0)*ce*x
                  - (73.0/16.0)*se );
}

/* Complex Gaussian 5. Real part */
static double CGAUSS5real(double x, double a, double b)
{
                /* c = 4 * sqrt(210) * exp(1/8) / ( 315 * pi^(1/4) ) */
      const double c = 0.1566240417643754;
      double x2, x3, x4, x5;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x3 = x2 * x;
      x4 = x2 * x2;
      x5 = x4 * x;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( -ce*x5 - (5.0/2.0)*se*x4 + (25.0/2.0)*ce*x3 + (65.0/4.0)*se*x2
                  - (365.0/16.0)*ce*x - (281.0/32.0)*se );
}

/* Complex Gaussian 5. Imaginary part */
static double CGAUSS5imag(double x, double a, double b)
{
                /* c = 4 * sqrt(210) * exp(1/8) / ( 315 * pi^(1/4) ) */
      const double c = 0.1566240417643754;
      double x2, x3, x4, x5;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x3 = x2 * x;
      x4 = x2 * x2;
      x5 = x4 * x;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( se*x5 - (5.0/2.0)*ce*x4 - (25.0/2.0)*se*x3 + (65.0/4.0)*ce*x2
                  + (365.0/16.0)*se*x - (281.0/32.0)*ce );
}

/* Complex Gaussian 6. Real part */
static double CGAUSS6real(double x, double a, double b)
{
                /* c = 8 * sqrt(1155) * exp(1/8) / ( 3465 * pi^(1/4) ) */
      const double c = 0.06678471580535175;
      double x2, x3, x4, x5, x6;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x3 = x2 * x;
      x4 = x2 * x2;
      x5 = x4 * x;
      x6 = x4 * x2;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( ce*x6 + 3.0*se*x5 - (75.0/4.0)*ce*x4 - (65.0/2.0)*se*x3
                  + (1095.0/16.0)*ce*x2 + (843.0/16.0)*se*x - (1741.0/64.0)*ce );
}

/* Complex Gaussian 6. Imaginary part */
static double CGAUSS6imag(double x, double a, double b)
{
                /* c = 8 * sqrt(1155) * exp(1/8) / ( 3465 * pi^(1/4) ) */
      const double c = 0.06678471580535175;
      double x2, x3, x4, x5, x6;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x3 = x2 * x;
      x4 = x2 * x2;
      x5 = x4 * x;
      x6 = x4 * x2;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( -se*x6 + 3.0*ce*x5 + (75.0/4.0)*se*x4 - (65.0/2.0)*ce*x3
                  - (1095.0/16.0)*se*x2 + (843.0/16.0)*ce*x + (1741.0/64.0)*se );
}

/* Complex Gaussian 7. Real part */
static double CGAUSS7real(double x, double a, double b)
{
                /* c = 8 * sqrt(30030) * exp(1/8) / ( 45045 * pi^(1/4) ) */
      const double c = 0.02619512070009515;
      double x2, x3, x4, x5, x6, x7;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x3 = x2 * x;
      x4 = x2 * x2;
      x5 = x4 * x;
      x6 = x4 * x2;
      x7 = x6 * x;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( -ce*x7 - (7.0/2.0)*se*x6 + (105.0/4.0)*ce*x5 + (455.0/8.0)*se*x4
                  - (2555.0/16.0)*ce*x3 - (5901.0/32.0)*se*x2 + (12187.0/64.0)*ce*x
                  + (8485.0/128.0)*se );
}

/* Complex Gaussian 7. Imaginary part */
static double CGAUSS7imag(double x, double a, double b)
{
                /* c = 8 * sqrt(30030) * exp(1/8) / ( 45045 * pi^(1/4) ) */
      const double c = 0.02619512070009515;
      double x2, x3, x4, x5, x6, x7;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x3 = x2 * x;
      x4 = x2 * x2;
      x5 = x4 * x;
      x6 = x4 * x2;
      x7 = x6 * x;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( se*x7 - (7.0/2.0)*ce*x6 - (105.0/4.0)*se*x5 + (455.0/8.0)*ce*x4
                  + (2555.0/16.0)*se*x3 - (5901.0/32.0)*ce*x2 - (12187.0/64.0)*se*x
                  + (8485.0/128.0)*ce );
}

/* Complex Gaussian 8. Real part */
static double CGAUSS8real(double x, double a, double b)
{
                /* c = 16 * sqrt(1001) * exp(1/8) / ( 45045 * pi^(1/4) ) */
      const double c = 0.009565105669341758;
      double x2, x3, x4, x5, x6, x7, x8;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x3 = x2 * x;
      x4 = x2 * x2;
      x5 = x4 * x;
      x6 = x4 * x2;
      x7 = x6 * x;
      x8 = x4 * x4;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( ce*x8 + 4.0*se*x7 - 35.0*ce*x6 - 91.0*se*x5 + (2555.0/8.0)*ce*x4
                  + (1967.0/4.0)*se*x3 - (12187.0/16.0)*ce*x2 - (8485.0/16.0)*se*x
                  + (57233.0/256.0)*ce );
}

/* Complex Gaussian 8. Imaginary part */
static double CGAUSS8imag(double x, double a, double b)
{
                /* c = 16 * sqrt(1001) * exp(1/8) / ( 45045 * pi^(1/4) ) */
      const double c = 0.009565105669341758;
      double x2, x3, x4, x5, x6, x7, x8;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x3 = x2 * x;
      x4 = x2 * x2;
      x5 = x4 * x;
      x6 = x4 * x2;
      x7 = x6 * x;
      x8 = x4 * x4;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( -se*x8 + 4.0*ce*x7 + 35.0*se*x6 - 91.0*ce*x5 - (2555.0/8.0)*se*x4
                  + (1967.0/4.0)*ce*x3 + (12187.0/16.0)*se*x2 - (8485.0/16.0)*ce*x
                  - (57233.0/256.0)*se );
}

/* Complex Gaussian 9. Real part */
static double CGAUSS9real(double x, double a, double b)
{
                /* c = 16 * sqrt(34034) * exp(1/8) / ( 765765 * pi^(1/4) ) */
      const double c = 0.003280804178061403;
      double x2, x3, x4, x5, x6, x7, x8, x9;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x3 = x2 * x;
      x4 = x2 * x2;
      x5 = x4 * x;
      x6 = x4 * x2;
      x7 = x6 * x;
      x8 = x4 * x4;
      x9 = x8 * x;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( -ce*x9 - (9.0/2.0)*se*x8 + 45.0*ce*x7 + (273.0/2.0)*se*x6
                  - (4599.0/8.0)*ce*x5 - (17703.0/16.0)*se*x4 + (36561.0/16.0)*ce*x3
                  + (76365.0/32.0)*se*x2 - (515097.0/256.0)*ce*x - (328753.0/512.0)*se );
}

/* Complex Gaussian 9. Imaginary part */
static double CGAUSS9imag(double x, double a, double b)
{
                /* c = 16 * sqrt(34034) * exp(1/8) / ( 765765 * pi^(1/4) ) */
      const double c = 0.003280804178061403;
      double x2, x3, x4, x5, x6, x7, x8, x9;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x3 = x2 * x;
      x4 = x2 * x2;
      x5 = x4 * x;
      x6 = x4 * x2;
      x7 = x6 * x;
      x8 = x4 * x4;
      x9 = x8 * x;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( se*x9 - (9.0/2.0)*ce*x8 - 45.0*se*x7 + (273.0/2.0)*ce*x6
                  + (4599.0/8.0)*se*x5 - (17703.0/16.0)*ce*x4 - (36561.0/16.0)*se*x3
                  + (76365.0/32.0)*ce*x2 + (515097.0/256.0)*se*x - (328753.0/512.0)*ce );
}

/* Complex Gaussian 10. Real part */
static double CGAUSS10real(double x, double a, double b)
{
                /* c = 32 * sqrt(323323) * exp(1/8) / ( 14549535 * pi^(1/4) ) */
      const double c = 0.001064433432433728;
      double x2, x3, x4, x5, x6, x7, x8, x9, x10;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x3 = x2 * x;
      x4 = x2 * x2;
      x5 = x4 * x;
      x6 = x4 * x2;
      x7 = x6 * x;
      x8 = x4 * x4;
      x9 = x8 * x;
      x10 = x8 * x2;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( ce*x10 + 5.0*se*x9 - (225.0/4.0)*ce*x8 - 195.0*se*x7 + (7665.0/8.0)*ce*x6
                  + (17703.0/8.0)*se*x5 - (182805.0/32.0)*ce*x4 - (127275.0/16.0)*se*x3
                  + (2575485.0/256.0)*ce*x2 + (1643765.0/256.0)*se*x - (2389141.0/1024.0)*ce );
}

/* Complex Gaussian 10. Imaginary part */
static double CGAUSS10imag(double x, double a, double b)
{
                /* c = 32 * sqrt(323323) * exp(1/8) / ( 14549535 * pi^(1/4) ) */
      const double c = 0.001064433432433728;
      double x2, x3, x4, x5, x6, x7, x8, x9, x10;
      double ce, se;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x3 = x2 * x;
      x4 = x2 * x2;
      x5 = x4 * x;
      x6 = x4 * x2;
      x7 = x6 * x;
      x8 = x4 * x4;
      x9 = x8 * x;
      x10 = x8 * x2;
      ce = cos(x/2) * exp(-x2/2);
      se = sin(x/2) * exp(-x2/2);
      return c * ( -se*x10 + 5.0*ce*x9 + (225.0/4.0)*se*x8 - 195.0*ce*x7 - (7665.0/8.0)*se*x6
                  + (17703.0/8.0)*ce*x5 + (182805.0/32.0)*se*x4 - (127275.0/16.0)*ce*x3
                  - (2575485.0/256.0)*se*x2 + (1643765.0/256.0)*ce*x + (2389141.0/1024.0)*se );
}

/* Paul 1. Real part */
static double PAUL1real(double x, double a, double b)
{
                /* c = -2 * sqrt(2) / sqrt(pi) */
      const double c = -1.595769121605731;
      double x2, x4;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      return c * x / (x4 + 2.0*x2 + 1.0);
}

/* Paul 1. Imaginary part */
static double PAUL1imag(double x, double a, double b)
{
                /* c = -sqrt(2) / sqrt(pi) */
      const double c = -0.7978845608028655;
      double x2, x4;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      return c * (x2-1.0) / (x4 + 2.0*x2 + 1.0);
}

/* Paul 2. Real part */
static double PAUL2real(double x, double a, double b)
{
                /* c = 2 * sqrt(2) / sqrt(3*pi) */
      const double c = 0.9213177319235614;
      double x2, x4, x6;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      x6 = x4 * x2;
      return c * (3.0*x2 - 1.0) / (x6 + 3.0*x4 + 3.0*x2 + 1.0);
}

/* Paul 2. Imaginary part */
static double PAUL2imag(double x, double a, double b)
{
                /* c = 2 * sqrt(2) / sqrt(3*pi) */
      const double c = 0.9213177319235614;
      double x2, x4, x6;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      x6 = x4 * x2;
      return c * x * (x2 - 3.0) / (x6 + 3.0*x4 + 3.0*x2 + 1.0);
}

/* Paul 3. Real part */
static double PAUL3real(double x, double a, double b)
{
                /* c = -16 / sqrt(5*pi) */
      const double c = -4.037012035232256;
      double x2, x4, x6, x8;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      x6 = x4 * x2;
      x8 = x4 * x4;
      return c * x * (x2 - 1.0) / (x8 + 4.0*x6 + 6.0*x4 + 4.0*x2 + 1.0);
}

/* Paul 3. Imaginary part */
static double PAUL3imag(double x, double a, double b)
{
                /* c = -4 / sqrt(5*pi) */
      const double c = -1.009253008808064;
      double x2, x4, x6, x8;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      x6 = x4 * x2;
      x8 = x4 * x4;
      return c * (x4 - 6.0*x2 + 1.0) / (x8 + 4.0*x6 + 6.0*x4 + 4.0*x2 + 1.0);
}

/* Paul 4. Real part */
static double PAUL4real(double x, double a, double b)
{
                /* c = 8 * sqrt(2) / sqrt(35*pi) */
      const double c = 1.078936850151577;
      double x2, x4, x6, x8, x10;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      x6 = x4 * x2;
      x8 = x4 * x4;
      x10 = x8 * x2;
      return c * (5.0*x4 - 10.0*x2 + 1.0) /
             (x10 + 5.0*x8 + 10.0*x6 + 10.0*x4 + 5.0*x2 + 1.0);
}

/* Paul 4. Imaginary part */
static double PAUL4imag(double x, double a, double b)
{
                /* c = 8 * sqrt(2) / sqrt(35*pi) */
      const double c = 1.078936850151577;
      double x2, x4, x6, x8, x10;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      x6 = x4 * x2;
      x8 = x4 * x4;
      x10 = x8 * x2;
      return c * x * (x4 - 10.0*x2 + 5.0) /
             (x10 + 5.0*x8 + 10.0*x6 + 10.0*x4 + 5.0*x2 + 1.0);
}

/* Paul 5. Real part */
static double PAUL5real(double x, double a, double b)
{
                /* c = -32 / ( 3 * sqrt(7*pi) ) */
      const double c = -2.274598598644513;
      double x2, x4, x6, x8, x10, x12;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      x6 = x4 * x2;
      x8 = x4 * x4;
      x10 = x8 * x2;
      x12 = x6 * x6;
      return c * x * (3.0*x4 - 10.0*x2 + 3.0) /
             (x12 + 6.0*x10 + 15.0*x8 + 20.0*x6 + 15.0*x4 + 6.0*x2 + 1.0);
}

/* Paul 5. Imaginary part */
static double PAUL5imag(double x, double a, double b)
{
                /* c = -16 / ( 3 * sqrt(7*pi) ) */
      const double c = -1.137299299322257;
      double x2, x4, x6, x8, x10, x12;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      x6 = x4 * x2;
      x8 = x4 * x4;
      x10 = x8 * x2;
      x12 = x6 * x6;
      return c * (x6 - 15.0*x4 + 15.0*x2 - 1.0) /
             (x12 + 6.0*x10 + 15.0*x8 + 20.0*x6 + 15.0*x4 + 6.0*x2 + 1.0);
}

/************************ PUBLIC ************************/

cwtwlet_t cwtwlets[] = {
    { "Haar",      -1.0,  1.0, &HAARreal,      NULL          },
    { "FrHat",     -1.5,  1.5, &FRHATreal,     NULL          },
    { "Splash",    -2.0,  2.0, &SPLASHreal,    NULL          },
    { "MexHat",    -5.0,  5.0, &MEXHATreal,    NULL          },
    { "Poisson",  -15.0, 15.0, &POISSONreal,   NULL          },
    { "Morlet",    -4.0,  4.0, &MORLETreal,    NULL          },
    { "CMorlet",   -4.0,  4.0, &CMORLETreal,   &CMORLETimag  },
    { "Shannon",  -10.0, 10.0, &SHANNONreal,   &SHANNONimag  },
    { "eShannon",  -4.0,  4.0, &ESHANNONreal,  &ESHANNONimag },
    { "Wave",      -5.0,  5.0, &GAUSS1real,    NULL          },
    { "Gauss1",    -5.0,  5.0, &GAUSS1real,    NULL          },
    { "Gauss2",    -5.0,  5.0, &GAUSS2real,    NULL          },
    { "Gauss3",    -5.0,  5.0, &GAUSS3real,    NULL          },
    { "Gauss4",    -5.0,  5.0, &GAUSS4real,    NULL          },
    { "Gauss5",    -5.0,  5.0, &GAUSS5real,    NULL          },
    { "Gauss6",    -5.0,  5.0, &GAUSS6real,    NULL          },
    { "Gauss7",    -5.0,  5.0, &GAUSS7real,    NULL          },
    { "Gauss8",    -5.0,  5.0, &GAUSS8real,    NULL          },
    { "Gauss9",    -5.0,  5.0, &GAUSS9real,    NULL          },
    { "Gauss10",   -5.0,  5.0, &GAUSS10real,   NULL          },
    { "CGauss1",   -5.0,  5.0, &CGAUSS1real,   &CGAUSS1imag  },
    { "CGauss2",   -5.0,  5.0, &CGAUSS2real,   &CGAUSS2imag  },
    { "CGauss3",   -5.0,  5.0, &CGAUSS3real,   &CGAUSS3imag  },
    { "CGauss4",   -5.0,  5.0, &CGAUSS4real,   &CGAUSS4imag  },
    { "CGauss5",   -5.0,  5.0, &CGAUSS5real,   &CGAUSS5imag  },
    { "CGauss6",   -5.0,  5.0, &CGAUSS6real,   &CGAUSS6imag  },
    { "CGauss7",   -5.0,  5.0, &CGAUSS7real,   &CGAUSS7imag  },
    { "CGauss8",   -5.0,  5.0, &CGAUSS8real,   &CGAUSS8imag  },
    { "CGauss9",   -5.0,  5.0, &CGAUSS9real,   &CGAUSS9imag  },
    { "CGauss10",  -5.0,  5.0, &CGAUSS10real,  &CGAUSS10imag },
    { "Paul1",    -12.0, 12.0, &PAUL1real,     &PAUL1imag    },
    { "Paul2",     -6.0,  6.0, &PAUL2real,     &PAUL2imag    },
    { "Paul3",     -4.0,  4.0, &PAUL3real,     &PAUL3imag    },
    { "Paul4",     -3.0,  3.0, &PAUL4real,     &PAUL4imag    },
    { "Paul5",     -3.0,  3.0, &PAUL5real,     &PAUL5imag    },

    /* Last entry */
    { "", 0.0, 0.0, NULL, NULL }
};
