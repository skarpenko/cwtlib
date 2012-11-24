#include <iostream>
#include "fixed.h"

using namespace std;

void main()
{
	fixed a, b, c;
	int d;

	a = int_to_fixed(100);
	b = int_to_fixed(10);
	d = a / b; //d будет содержать 10 в формате int
	cout << d << endl;
	c = a / 10; // с будет содержать 10 в формате fixed
	cout << fixed_to_int(c) << endl;

	c = a + b; //c будет содержать 110.0 в формате fixed
	cout << fixed_to_int(c) << endl;

	cin >> a;
}
