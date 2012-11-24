/***********************************************************************
                     Continuous Wavelet Transform

 Author  : Stepan V.Karpenko
 Date    : 01-04-2004
 Comments:
 History :
***********************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "cwtwlets.h"
#include "cwt.h"
#include "fixcwt.h"
char __id__string[] = "\nCopyright(C) Carp Co.,\nBased on CWTLib v1.3 by Stepan Karpenko.\n";
char __version__string[] = "Version 1.4";

/* #define NOOUT */                /* Uncomment for silent mode */
/* #define WINMAIN */              /* Uncomment when console not needed */
/* #define CORAMPL */              /* Uncomment for CWT normalization */
#define EXENAME "cwtmorlet.exe"    /* Name of executable file */

#define AMIN 1     /* A min */
#define ASTP 0.1     /* A step */
#define AMAX 17 //16   /* A max */
#define BSTP 0.1     /* B step */
#define PREC 30    /* Precision */

cwtwlet_t *WAVELET = &cwtwlets[MORLET];  /* Default wavelet */

long N = 36;     /* Signal length */
double *s;

clock_t msecs;


int main(int argc, char *argv[])
{
     FILE *fh;
     long i, j, dx, dy;
     cwt_t wt;
     double a, b;
     float temp;

     if(argc<3) {
          #ifndef NOOUT
            printf("filename missing\n");
          #endif
          return 1;
     }
     fh = fopen(argv[1], "r");
     if(!fh) {
          #ifndef NOOUT
            printf("can't open file: %s\n", argv[1]);
          #endif
          return 1;
     }

     if(argc == 5) {
         if( atol(argv[4]) )
           N = atol(argv[4]);
     }
     s = (double *)malloc( N*sizeof(double) );
     if(!s) {
          #ifndef NOOUT
            printf("Not enought memory!\n");
          #endif
          return 1;
     }

     for(i=0; i<N; i++) { fscanf(fh, "%f", &temp); s[i] = temp; }
     fclose(fh);
     fh = fopen(argv[2], "w");
     if(!fh) {
          #ifndef NOOUT
            printf("can't open file: %s\n", argv[2]);
          #endif
          return 1;
     }
     /* Choose wavelet */
     if(argc == 4) {
          i = 0;
          while( cwtwlets[i].wname[0] != 0 ) {
            if( !strcmp(argv[3], cwtwlets[i].wname) )
               break;
            i++;
          }
          if( cwtwlets[i].wname[0] != 0 )
             WAVELET = &cwtwlets[i];
          else {
             #ifndef NOOUT
              printf("Unknown wavelet");
             #endif
             return 1;
          }
     }
    
   #ifndef NOOUT
     printf("Performing wavelet transform...");
   #endif
     msecs = clock();
     if( cwto3(s, N, AMIN, ASTP, AMAX, BSTP, PREC, WAVELET, REAL, 128000, &wt) )
     {
         #ifndef NOOUT
           printf("error performing wavelet transform!\n");
         #endif
         return 1;
     }

   #ifndef NOOUT
     printf("Time elapsed: %u sec.\n", (clock()-msecs)/CLK_TCK);
      
     printf("\nTransform info:\nAmin\t\t= %f\nAstep\t\t= %f\nAmax\t\t= %f\n",
        wt.amin, wt.astep, wt.amax);
     printf("Bstep\t\t= %f\nSignal length\t= %u\n", wt.bstep, wt.siglen);
     printf("CWT array dimension: %ux%u\n", wt.rows, wt.cols);
     printf("Wavelet name: %s\n\n", WAVELET->wname);
   #endif

   #ifndef NOOUT
     printf("Fixing cwt boundaries...\n");
   #endif
     fixcwt(&wt, s, PREC, WAVELET);

   #ifdef CORAMPL
   #ifndef NOOUT
      printf("Normalizing amplitudes...\n");
   #endif
     for (dy = 0, a = wt.amin; dy < wt.rows; dy++, a+=wt.astep)
     {
         for ( dx = 0, b = 0.0; dx < wt.cols; dx++, b+=wt.bstep)
         {
           wt.cwt[dy][dx] *= 2.0 / ( sqrt(2 * a * PI) );
         }
     }
   #endif

   #ifndef NOOUT
     printf("Storing data...\n");
   #endif
     for(i=0; i<wt.rows; i++) {
        for(j=0; j<wt.cols; j++) {
            fprintf(fh, " %.15f", wt.cwt[i][j]);
        }
        fprintf(fh, "\n");
     }
     fclose(fh);

     #ifndef NOOUT
       printf("Freeing allocated memory...\n");
     #endif
     free_cwt(&wt);
     free(s);

     return 0;
}

#ifdef WINMAIN
/* WinMain*/
int __stdcall WinMain( void *hInstance, void *hPrevInstance, char *lpCmdLine, int iCmdShow )
{
    char *args[4] = {""EXENAME""};
    int argc = 3;
    int i=0;

    /* Parse command line string */
    if(!lpCmdLine[0]) return 1;

    /* First argument */
    args[1] = lpCmdLine;

    /* Find second argument */
    while(lpCmdLine[i] != 32 && lpCmdLine[i] != 0)
       i++;

    if(!i || !lpCmdLine[i])
      return 1;
    else {
      lpCmdLine[i] = 0;
      args[2] = &lpCmdLine[++i];
    }

    /* Find third argument */
    while(lpCmdLine[i] != 32 && lpCmdLine[i] != 0)
       i++;

    if(lpCmdLine[i] == 32) {
      lpCmdLine[i] = 0;
      args[3] = &lpCmdLine[++i];
      argc++;
    }

    /* Find fourth argument */
    while(lpCmdLine[i] != 32 && lpCmdLine[i] != 0)
       i++;

    if(lpCmdLine[i] == 32) {
      lpCmdLine[i] = 0;
      args[4] = &lpCmdLine[++i];
      argc++;
    }

    /* All other arguments not needed */
    while(lpCmdLine[i] != 32 && lpCmdLine[i] != 0)
       i++;

    if(lpCmdLine[i] == 32) lpCmdLine[i] = 0;

    /* Call main */
    return main(argc, args);
}
#endif
