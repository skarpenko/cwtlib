#include <iostream>
using namespace std;


template<class T, unsigned _bits>
class basic_fixed {
    T value;
  public:
    typedef T basic_type;
    static const unsigned bits_used = _bits;

    basic_fixed() { value = T(); }
    template<class U> basic_fixed(U v) { value = static_cast<T>(v*(1<<_bits)); }

// add ctr for conver from other bits and types
template<class U, unsigned _bits2>
basic_fixed(const basic_fixed<U, _bits2>& a) { int b = a.value; };

    basic_fixed(const basic_fixed& a) { value = a.value; }
    
    basic_fixed& operator=(const basic_fixed& a)
    {
        if(this != &a)
            value = a.value;
        return *this;
    }

    basic_fixed operator+(const basic_fixed& a)
    {
        value += a.value;
        return *this;
    }

    basic_fixed operator-(const basic_fixed& a)
    {
        value += a.value;
        return *this;
    }

    basic_fixed operator*(const basic_fixed& a)
    {
        value = (value*a.value) >> _bits;
        return *this;
    }

    basic_fixed operator/(const basic_fixed& a)
    {
        //value = (value/a.value) << _bits;
        value = (value<<_bits)/a.value;
        return *this;
    }

    // for DEBUG
    void dbg_print_val() { cout << value << endl; }
    void dbg_print_parts()
    {
        T up = value >> _bits;
        T down = value & ((1<<_bits) - 1);
        
        cout << up << " : " << down / (double)(1<<_bits)  << endl;
    }
    double toDouble()
    {
        return static_cast<double>(value) / (1<<_bits);
    }
    
    T raw() { return value; }
    
};

// Geron sqrt
// Xn+1 = 1/2 * (Xn + a/Xn)
template<class T, unsigned _bits>
basic_fixed<T, _bits> sqrt(basic_fixed<T, _bits> v)
{
    typedef basic_fixed<T, _bits> fixed_type;
    /*const*/ fixed_type v12 = fixed_type(1) / fixed_type(2);
    fixed_type r = 1.0;
    int i;
    double rr=1.0, vv=9.0;
    cout << "v12 = " << v12.toDouble() << endl << "v = " << v.toDouble() << " r= " << r.toDouble() << endl;
    for(i=0; i<10; ++i) {
        cout << i << "/1: r = " << r.toDouble() << "\t\tv/r = " << (v/r).toDouble() << "\t\tr + v/r = " << (r + v/r).toDouble() << endl;
        r = v12 * (r + v / r);
        cout << i << "/2: r = " << rr << "\t\tv/r = " << vv/rr << "\t\tr + v/r = " << rr+vv/rr << endl;
        rr = 0.5 * (rr + vv / rr);
     }
    
    return r;
}

typedef basic_fixed<long long, 12> fixed_t;
typedef basic_fixed<int, 12> fixed2_t;


int main()
{
    fixed_t a = 2.55;
    fixed_t b = 3.76;
//    fixed2_t c = 3.1;
//    fixed_t z = c;

    cout << sizeof(fixed_t) << endl;
    b.dbg_print_parts(); return 0;
    cout << (a+b).toDouble() << endl;
    
    a = a / b;
    
    cout << 2.55/3.76 << endl << a.toDouble() << endl;

    fixed_t v = 9;
    fixed_t sq;
    
    sq = sqrt(v);
    cout << "Sqrt = " << sq.toDouble() << endl;

    cout << (fixed_t(0.25) / fixed_t(0.05)).toDouble() << endl;

    return 0;
}
