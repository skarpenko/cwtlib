#include <stdio.h>
#include <math.h>
#define TINY 0.00000000000001;

// c:=sqrt(2*Pi): w:=2*Pi*2:
// psi:=c*x*exp(-x^2) * (cos(w*x) + I*sin(w*x));

/* Norlet real part */
static double NORLETreal(double x, double a, double b)
{
      const double w = 2 * PI * 0.8;
      const double c = 2 * PI;
      double x2;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      return c * x * exp(-x2) * cos(w*x);
}

/* Norlet imaginary part */
static double NORLETimag(double x, double a, double b)
{
      const double w = 2 * PI * 0.8;
      const double c = 2 * PI;
      double x2;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      return c * x * exp(-x2) * sin(w*x);
}

int main() {
FILE *T, *Fr, *Fi;
double t;
T=fopen("t.txt", "w");
Fr=fopen("fr.txt", "w");
Fi=fopen("fi.txt", "w");

for(t=-6; t<=6; t+=0.01) {
   fprintf(T, "%f\n", t);
   fprintf(Fr, "%f\n", NORLETreal(t,1,0));
   fprintf(Fi, "%f\n", NORLETimag(t,1,0));
}

fclose(T); fclose(Fr); fclose(Fi);
}


