#include <stdio.h> //Нужна для функции getch();
#include <conio.h>
#include <math.h>

//#define M_PI ((float)3.141592653589793)

double mytan(double x, int i);
double mysin(double x, int i);
double mycos(double x, int i);
//double prevedpi(double x, int i);
double prevedpi(double x);

void main()
    {  int i;
	   printf("vvedite ygol\n");
	   //scanf("%c",&i);
	   //printf("tangens %c\n sinus %c\n cosinus %c", mytan, mysin, mycos);
       scanf("%d",&i);
	printf("tangens %lf, sinus %lf, cosinus %lf\n", (float)(mytan(i, 17)), (float)(mycos(i, 17)), (float)(mysin(i, 17)));
	printf("tangens %lf, sinus %lf, cosinus %lf\n", tan(i), cos(i), sin(i));
       getch();   
    }

//Тангенс на +-1/4 пи, i-количество звеньев в дроби (выше 15-17 бессмысленно точночть близка к макс)
double mytan(double x, int i)
   {
   x=prevedpi(x);
   
   double ret=x*x;//То, что возвращаем    
       for(; i>=0; i--)
          {  
               ret=x/( (2*i+1)  - x*ret  );    
          }
   return ret;
   }

//Тангенс на +-1/4 пи, i-количество звеньев в дроби тангенса
double mysin(double x, int i)
   {
   x=prevedpi(x);
   //В худшем случае аргумент близок к +-пи. Поэтому делим его на 4.
   double t=mytan(x/4, i);
   //Теперь получаем синус и косинус углов х/2
   double s=2*t/(1+t*t), c=(1-t*t)/(1+t*t);
   //Переобъявляем тангенс. теперь он у нас тангенс x/2
   t=s/c;
   //Находим искомый синус
   
   return ((2*t)/(1+t*t));
   }
double mycos(double x, int i)
   {
   x=prevedpi(x);
   //В худшем случае аргумент близок к +-пи. Поэтому делим его на 4.
   double t=mytan(x/4, i);
   //Теперь получаем синус и косинус углов х/2
   double s=2*t/(1+t*t), c=(1-t*t)/(1+t*t);
   //Переобъявляем тангенс. теперь он у нас тангенс x/2
   t=s/c;
   //Находим искомый koсинус
   
   return ((1-t*t)/(1+t*t));
   }


//Преведение по пи
double prevedpi(double x)
   {
   //сколько раз укладывается пи    
   int a=((int)(x/(M_PI)));   
       
       if(x/M_PI > 1.0 || x/M_PI < - 1.0 )
          {
          //Если укладывается четное число раз       
          if(a%2 == 0)
           x = x - a*M_PI;
          else
           x =-(x - a*M_PI);
           
          }
   return x;
   }
