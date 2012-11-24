/***********************************************************************
              Continuous Wavelet Transform Demo Program

 Author  : Stepan V.Karpenko
 Date    : 01-04-2004
 Comments:
 History :
***********************************************************************/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include "cwtwlets.h"
#include "cwt.h"

#define WAVELET cwtwlets[MEXHAT].real

#define N 1024
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
     if( cwt(s, N, 1, 1, 128, 1, 2, WAVELET, &wt) )
     {
         printf("error performing wavelet transform!\n");
         return 1;
     }
     printf("Time elapsed: %u sec.\n", (clock()-msecs)/CLK_TCK);

     printf("\nTransform info:\nAmin\t\t= %f\nAstep\t\t= %f\nAmax\t\t= %f\n",
        wt.amin, wt.astep, wt.amax);
     printf("Bstep\t\t= %f\nSignal length\t= %u\n", wt.bstep, wt.siglen);
     printf("CWT array dimension: %ux%u\n\n", wt.rows, wt.cols);

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
