/***********************************************************************
           Wavelet Functions for Continuous Wavelet Transform

 Author  : Stepan V.Karpenko
 Date    : 01-04-2004
 Comments:
 History :
   09-04-2004 Stepan V.Karpenko
    Added cwtwlets[] array and cwtwlet_t type.
***********************************************************************/

#ifndef _CWTWLETS_H_
#define _CWTWLETS_H_

#include "cwt.h"

/* Definitions */
#ifndef PI
  #define PI  3.14159265358979323846
#endif
#define TINY 1E-20
#define WNAMELEN 10

/* Wavelet structure */
typedef struct {
        char wname[WNAMELEN];  /* Wavelet name */
        psi_t *real;           /* Real part */
        psi_t *imag;           /* Imaginary part, if not exist must be NULL */
} cwtwlet_t;

#ifdef __cplusplus
extern "C" {
#endif

/* Enumeration of wavelet entries */
enum { MEXHAT, MORLET };
extern cwtwlet_t cwtwlets[]; /* Array of wavelets */

#ifdef __cplusplus
}
#endif

#endif
