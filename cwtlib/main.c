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
 Title   : Continuous Wavelet Transform Demo Program
 Author  : Stepan V.Karpenko
 Date    : 01-04-2004
 Comments:
 History :
   31-05-2004 Stepan V.Karpenko
    Added support for optimized versions of cwt().
   08-10-2004 Stepan V.Karpenko
    CLK_TCK replaced with CLOCKS_PER_SEC.
    Added support for cwto3().
***********************************************************************/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include "cwtwlets.h"
#include "cwt.h"

/* CWT version */
/** #define NOOPT **/ /* Without optimization */
/** #define OPT1 **/  /* Optimized version 1 */
/** #define OPT2 **/  /* Optimized version 2 */
#define OPT3          /* Optimized version 3 */

#if defined(OPT2) || defined(OPT3)
  #define NPOINTS 60000 /* Points number for wavelet precompution */
#endif


#define WAVELET cwtwlets[MEXHAT]
#define AMIN 1    /* A min  */
#define ASTP 1    /* A step */
#define AMAX 128  /* A max  */
#define BSTP 1    /* B step */
#define PREC 2    /* Precision */
#define N 1024    /* Signal length */

double s[N];

clock_t msecs;


int main(int argc, char *argv[])
{
     FILE *fh;
     long i, j;
     cwt_t wt;
     float temp;

     if(argc<2) {
          printf("filename missing\n");
          return 1;
     }
     fh = fopen(argv[1], "r");
     if(!fh) {
          printf("can't open file\n");
          return 1;
     }
     for(i=0; i<N; i++) { fscanf(fh, "%f", &temp); s[i] = temp; }
     fclose(fh);
     fh = fopen("output.txt", "w");

     printf("Performing wavelet transform...");
     msecs = clock();
    #ifdef NOOPT
     if( cwt(s, N, AMIN, ASTP, AMAX, BSTP, PREC, WAVELET.real, &wt) )
    #endif
    #ifdef OPT1
     if( cwto1(s, N, AMIN, ASTP, AMAX, BSTP, PREC, &WAVELET, REAL, &wt) )
    #endif
    #ifdef OPT2
     if( cwto2(s, N, AMIN, ASTP, AMAX, BSTP, PREC, &WAVELET, REAL, NPOINTS, &wt) )
    #endif
    #ifdef OPT3
     if( cwto3(s, N, AMIN, ASTP, AMAX, BSTP, PREC, &WAVELET, REAL, NPOINTS, &wt) )
    #endif
     {
         printf("error performing wavelet transform!\n");
         return 1;
     }
     printf("Elapsed time: %u sec.\n", (clock()-msecs)/CLOCKS_PER_SEC);

     printf("\nTransform info:\nAmin\t\t= %f\nAstep\t\t= %f\nAmax\t\t= %f\n",
        wt.amin, wt.astep, wt.amax);
     printf("Bstep\t\t= %f\nSignal length\t= %u\n", wt.bstep, wt.siglen);
     printf("CWT array dimension: %ux%u\n", wt.rows, wt.cols);
     printf("Wavelet: %s\n\n", wt.wname);

     printf("Storing data.\n");
     for(i=0; i<wt.rows; i++) {
        for(j=0; j<wt.cols; j++) {
            fprintf(fh, " %.15f", wt.cwt[i][j]);
        }
        fprintf(fh, "\n");
     }
     fclose(fh);

     printf("Freeing allocated memory.\n");
     free_cwt(&wt);

     return 0;
}
