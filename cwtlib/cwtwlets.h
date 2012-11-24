/***********************************************************************
           Wavelet Functions for Continuous Wavelet Transform

 Author  : Stepan V.Karpenko
 Date    : 01-04-2004
 Comments:
 History :
***********************************************************************/

#ifndef _CWTWLETS_H_
#define _CWTWLETS_H_

/* Definitions */
#ifndef PI
  #define PI  3.14159265358979323846
#endif
#define TINY 1E-20

#ifdef __cplusplus
extern "C" {
#endif

/* Mexican Hat */
double MEXHAT(double x, double a, double b);

/* Morlet real part */
double MORLETreal(double x, double a, double b);

/* Morlet imaginary part */
double MORLETimag(double x, double a, double b);

#ifdef __cplusplus
}
#endif

#endif
