// Error function compute

#include <math.h>
#include <stdio.h>
#include <time.h>
#include "cwtwlets.h"
#include "cwt.h"

#define NPOINTS 1000000

#define N1 100
#define N2 10000
#define DN 100

#define WAVELET cwtwlets[MORLET]
#define AMIN 1    /* A min  */
#define ASTP 1    /* A step */
#define AMAX 256  /* A max  */
#define BSTP 1    /* B step */
#define PREC 4    /* Precision */
#define N 256     /* Signal length */

double s[N];

double mid(cwt_t wt) {
  long i, j;
  double r=0;

     for(i=0; i<wt.rows; i++)
        for(j=0; j<wt.cols; j++)
            r+=fabs(wt.cwt[i][j]);
 
   return r/(wt.rows*wt.cols);
}


int main(int argc, char *argv[])
{
     FILE *fh, *fh2;
     long i, j, n;
     cwt_t wt1, wt2;
     float temp;
     double mid1, mid2;

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

    printf("WT N1\n");
//    cwt(s, N, AMIN, ASTP, AMAX, BSTP, PREC, WAVELET.real, &wt1);
    cwto2(s, N, AMIN, ASTP, AMAX, BSTP, PREC, &WAVELET, REAL, NPOINTS, &wt1);
    mid1 = mid(wt1);
    printf("Mid1=%f\n", mid1);

    fh = fopen("err.txt", "w");
    fh2 =  fopen("n.txt", "w");
    printf("WT N2\n");
    for(n=N1; n<=N2; n+=DN) {
      cwto2(s, N, AMIN, ASTP, AMAX, BSTP, PREC, &WAVELET, REAL, n, &wt2);
      for(i=0; i<wt2.rows; i++)
         for(j=0; j<wt2.cols; j++)
            wt2.cwt[i][j]-=wt1.cwt[i][j];

      mid2 = mid(wt2);
      fprintf(fh, "%f\n", (mid2*100.0)/mid1);
      fprintf(fh2, "%d\n", n);
      free_cwt(&wt2);
      printf("n=%d, Mid2=%f, %%=%f\r", n, mid2, (mid2*100.0)/mid1);
    }
    fclose(fh); fclose(fh2);
    free_cwt(&wt1);

     return 0;
}
