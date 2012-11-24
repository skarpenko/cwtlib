#include <windows.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <list>
#include <math.h>
#include <cstdlib>
#include "fixed.h"
using namespace std;

typedef fixed_point<int, 16> fixed_t;

//const fixed_t F_PI = frac(22, 7);
//const fixed_t F_PI = frac<long long, 16>(355, 113);
const fixed_t F_PI = fixed_t(355) / fixed_t(113);

fixed_t prevedpi(fixed_t x);
fixed_t mytan(fixed_t x, int i=10);
fixed_t mysin(fixed_t x, int i=10);
fixed_t mycos(fixed_t x, int i=10);

void test_speed();

int main()
{
    srand( GetTickCount() );
/*
    double ang = 3.14;
    fixed_t a = 10;
    fixed_t b = 11;
    fixed_t c;
    
    c = a + b / 2;
    a = 9;
    c = sqrt(a, 10);
    
    cout << "PI = " << F_PI << endl;
    cout << c << endl;
   
    cout << "D: Tan = " << tan(ang) << "  Sin = " << sin(ang) << "  Cos = " << cos(ang) << endl;
    cout << "X: Tan = " << mytan(ang) << "  Sin = " << mysin(ang) << "  Cos = " << mycos(ang) << endl;
*/
cout << sqrt(fixed_t(9.0)) << endl; return 0;
    cout << endl << endl << "Ololo Test" << endl << endl;
    test_speed();

    return 0;
}

inline double GenVal()
{
    return static_cast<double>( rand() ) / RAND_MAX * 10.0;
}

template<class T>
void PrintItem(T item)
{
    cout << item << " ";
}

template<class T>
T Transf1(T item)
{
    return sqrt(item);
}

template<class T>
T Transf2(T item)
{
    return sin(item);
}

template<>
fixed_t Transf2<fixed_t> (fixed_t item)
{
    return mysin(item);
}

template<class T>
T Transf3(T item)
{
    T r = item * 3 + 2;
    return (r-2)/3;
}

void test_speed()
{
    const int N = 1000;
    const int T = 10000;
    typedef vector<float> dlist;
    typedef vector<fixed_t> flist;
    long dt;
    long i;
    
    dlist d(N);
    flist f(N);
    
    generate(d.begin(), d.end(), GenVal);
    generate(f.begin(), f.end(), GenVal);

//    for_each(d.begin(), d.end(), PrintItem<double>);
//    for_each(f.begin(), f.end(), PrintItem<fixed_t>);

    dt = GetTickCount();
    i = 0;
    while( GetTickCount() - dt < T ) {
        transform(d.begin(), d.end(), d.begin(), Transf1<double>);
//        transform(d.begin(), d.end(), d.begin(), Transf2<double>);
//        transform(d.begin(), d.end(), d.begin(), Transf3<double>);
        i++;
    }
    cout << "Double I = " << i << endl;

    dt = GetTickCount();
    i = 0;
    while( GetTickCount() - dt < T ) {
        transform(f.begin(), f.end(), f.begin(), Transf1<fixed_t>);
//        transform(f.begin(), f.end(), f.begin(), Transf2<fixed_t>);
//        transform(f.begin(), f.end(), f.begin(), Transf3<fixed_t>);
        i++;
    }
    cout << "Fixed I = " << i << endl;
}

//Преведение по пи
fixed_t prevedpi(fixed_t x)
{
   //сколько раз укладывается пи    
   int a=fixed_cast<int>(x/(F_PI));
       
       if(x/F_PI > 1 || x/F_PI < -1)
          {
          //Если укладывается четное число раз       
          if(a%2 == 0)
           x = x - a*F_PI;
          else
           x =-(x - a*F_PI);
           
          }
   return x;
}

fixed_t mytan(fixed_t x, int i)
   {
   x=prevedpi(x);
   
   fixed_t ret=x*x;//То, что возвращаем    
       for(; i>=0; i--)
          {  
               ret=x/( (2*i+1)  - x*ret  );    
          }
   return ret;
   }
   
//Тангенс на +-1/4 пи, i-количество звеньев в дроби тангенса
fixed_t mysin(fixed_t x, int i)
   {
   x=prevedpi(x);
   //В худшем случае аргумент близок к +-пи. Поэтому делим его на 4.
   fixed_t t=mytan(x/4, i);
   //Теперь получаем синус и косинус углов х/2
   fixed_t s=2*t/(1+t*t), c=(1-t*t)/(1+t*t);
   //Переобъявляем тангенс. теперь он у нас тангенс x/2
   t=s/c;
   //Находим искомый синус
   
   return ((2*t)/(1+t*t));
   }
fixed_t mycos(fixed_t x, int i)
   {
   x=prevedpi(x);
   //В худшем случае аргумент близок к +-пи. Поэтому делим его на 4.
   fixed_t t=mytan(x/4, i);
   //Теперь получаем синус и косинус углов х/2
   fixed_t s=2*t/(1+t*t), c=(1-t*t)/(1+t*t);
   //Переобъявляем тангенс. теперь он у нас тангенс x/2
   t=s/c;
   //Находим искомый koсинус
   
   return ((1-t*t)/(1+t*t));
   }

