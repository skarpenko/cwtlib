#ifndef __FIXED_POINT__
#define __FIXED_POINT__

#include <iostream> //DEBUG

template<class T, unsigned FRAC_BITS>
class fixed_point {
    T value;
  public:
    typedef T basic_type;
    static const unsigned frac_bits = FRAC_BITS;

    fixed_point() { value = T(); }
    template<class U> fixed_point(U v) { value = static_cast<T>(v*(1 << frac_bits)); }
    fixed_point(const fixed_point& a) { value = a.value; }
    
    fixed_point& operator=(const fixed_point& a)
    {
        if(this != &a)
            value = a.value;
        return *this;
    }

    fixed_point& operator+=(const fixed_point& a)
    {
        value += a.value;
        return *this;
    }

    fixed_point& operator-=(const fixed_point& a)
    {
        value -= a.value;
        return *this;
    }

    fixed_point& operator*=(const fixed_point& a)
    {
        value = (value*a.value) >> frac_bits;
        return *this;
    }

    fixed_point& operator/=(const fixed_point& a)
    {
        //value = (value/a.value) << _bits;
        value = (value<<frac_bits)/a.value;
//        value = (value << (frac_bits/2) ) / (a.value >> (frac_bits/2) );
        return *this;
    }

    fixed_point& operator-()
    {
        value = -value;
        return *this;
    }

    fixed_point& operator+()
    {
        return *this;
    }
    
    // for DEBUG
    void dbg_print_val() { std::cout << value << std::endl; }
    void dbg_print_parts()
    {
        T up = value >> frac_bits;
        T down = value & ((1<<frac_bits) - 1);
        
        std::cout << up << " : " << down / (double)(1<<frac_bits)  << std::endl;
    }
    double toDouble()
    {
        return static_cast<double>(value) / (1<<frac_bits);
    }
/* CAREFUL WITH THIS ONE!
    operator double() const
    {
        return static_cast<double>(value) / (1<<frac_bits);
    }
*/
    T raw() const { return value; }
};

/* correct block */
template<class T, unsigned FRAC_BITS> inline fixed_point<T, FRAC_BITS>
operator+(const fixed_point<T, FRAC_BITS>& a1, const fixed_point<T, FRAC_BITS>& a2)
{
    return fixed_point<T, FRAC_BITS>(a1) += a2;
}

template<class T, unsigned FRAC_BITS, class U> inline fixed_point<T, FRAC_BITS>
operator+(const fixed_point<T, FRAC_BITS>& a1, const U a2)
{
    return a1 + fixed_point<T, FRAC_BITS>(a2);
}

template<class T, unsigned FRAC_BITS, class U> inline fixed_point<T, FRAC_BITS>
operator+(const U a1, const fixed_point<T, FRAC_BITS>& a2)
{
    return fixed_point<T, FRAC_BITS>(a1) + a2;
}
/////

template<class T, unsigned FRAC_BITS> inline fixed_point<T, FRAC_BITS>
operator-(const fixed_point<T, FRAC_BITS>& a1, const fixed_point<T, FRAC_BITS>& a2)
{
    return fixed_point<T, FRAC_BITS>(a1) -= a2;
}

template<class T, unsigned FRAC_BITS, class U> inline fixed_point<T, FRAC_BITS>
operator-(const fixed_point<T, FRAC_BITS>& a1, const U a2)
{
    return a1 - fixed_point<T, FRAC_BITS>(a2);
}

template<class T, unsigned FRAC_BITS, class U> inline fixed_point<T, FRAC_BITS>
operator-(const U a1, const fixed_point<T, FRAC_BITS>& a2)
{
    return fixed_point<T, FRAC_BITS>(a1) - a2;
}

template<class T, unsigned FRAC_BITS> inline fixed_point<T, FRAC_BITS>
operator*(const fixed_point<T, FRAC_BITS>& a1, const fixed_point<T, FRAC_BITS>& a2)
{
    return fixed_point<T, FRAC_BITS>(a1) *= a2;
}

template<class T, unsigned FRAC_BITS, class U> inline fixed_point<T, FRAC_BITS>
operator*(const fixed_point<T, FRAC_BITS>& a1, const U a2)
{
    return fixed_point<T, FRAC_BITS>(a1) *= a2;
}

template<class T, unsigned FRAC_BITS, class U> inline fixed_point<T, FRAC_BITS>
operator*(const U a1, const fixed_point<T, FRAC_BITS>& a2)
{
    return fixed_point<T, FRAC_BITS>(a1) *= a2;
}

template<class T, unsigned FRAC_BITS> inline fixed_point<T, FRAC_BITS>
operator/(const fixed_point<T, FRAC_BITS>& a1, const fixed_point<T, FRAC_BITS> a2)
{
    return fixed_point<T, FRAC_BITS>(a1) /= a2;
}

template<class T, unsigned FRAC_BITS, class U> inline fixed_point<T, FRAC_BITS>
operator/(const fixed_point<T, FRAC_BITS>& a1, const U a2)
{
    return fixed_point<T, FRAC_BITS>(a1) /= a2;
}

template<class T, unsigned FRAC_BITS, class U> inline fixed_point<T, FRAC_BITS>
operator/(const U a1, const fixed_point<T, FRAC_BITS>& a2)
{
    return fixed_point<T, FRAC_BITS>(a1) /= a2;
}

template<class T, unsigned FRAC_BITS> inline
bool operator==(const fixed_point<T, FRAC_BITS>& a1, const fixed_point<T, FRAC_BITS>& a2)
{
    return a1.raw() == a2.raw();
}

template<class T, unsigned FRAC_BITS, class U> inline
bool operator==(const fixed_point<T, FRAC_BITS>& a1, const U a2)
{
    return a1.raw() == fixed_point<T, FRAC_BITS>(a2).raw();
}

template<class T, unsigned FRAC_BITS, class U> inline
bool operator==(const U a1, const fixed_point<T, FRAC_BITS>& a2)
{
    return fixed_point<T, FRAC_BITS>(a1).raw() == a2.raw();
}

template<class T, unsigned FRAC_BITS> inline
bool operator!=(const fixed_point<T, FRAC_BITS>& a1, const fixed_point<T, FRAC_BITS>& a2)
{
    return a1.raw() != a2.raw();
}

template<class T, unsigned FRAC_BITS, class U> inline
bool operator!=(const fixed_point<T, FRAC_BITS>& a1, const U a2)
{
    return a1.raw() != fixed_point<T, FRAC_BITS>(a2).raw();
}

template<class T, unsigned FRAC_BITS, class U> inline
bool operator!=(const U a1, const fixed_point<T, FRAC_BITS>& a2)
{
    return fixed_point<T, FRAC_BITS>(a1).raw() != a2.raw();
}

template<class T, unsigned FRAC_BITS> inline
bool operator<(const fixed_point<T, FRAC_BITS>& a1, const fixed_point<T, FRAC_BITS>& a2)
{
    return a1.raw() < a2.raw();
}

template<class T, unsigned FRAC_BITS, class U> inline
bool operator<(const fixed_point<T, FRAC_BITS>& a1, const U a2)
{
    return a1.raw() < fixed_point<T, FRAC_BITS>(a2).raw();
}

template<class T, unsigned FRAC_BITS, class U> inline
bool operator<(const U a1, const fixed_point<T, FRAC_BITS>& a2)
{
    return fixed_point<T, FRAC_BITS>(a1).raw() < a2.raw();
}

template<class T, unsigned FRAC_BITS> inline
bool operator>(const fixed_point<T, FRAC_BITS>& a1, const fixed_point<T, FRAC_BITS>& a2)
{
    return a1.raw() > a2.raw();
}

template<class T, unsigned FRAC_BITS, class U> inline
bool operator>(const fixed_point<T, FRAC_BITS>& a1, const U a2)
{
    return a1.raw() > fixed_point<T, FRAC_BITS>(a2).raw();
}

template<class T, unsigned FRAC_BITS, class U> inline
bool operator>(const U a1, const fixed_point<T, FRAC_BITS>& a2)
{
    return fixed_point<T, FRAC_BITS>(a1).raw() > a2.raw();
}

// Geron sqrt
// Xn+1 = 1/2 * (Xn + a/Xn)
template<class T, unsigned FRAC_BITS> inline
fixed_point<T, FRAC_BITS> sqrt(const fixed_point<T, FRAC_BITS>& v, unsigned ni = 10)
{
    typedef fixed_point<T, FRAC_BITS> fixed_point_type;
    fixed_point_type r = 1.0;

    while(ni--) r = (r + v / r) / fixed_point_type(2);

    return r;
}

template<class T, unsigned FRAC_BITS> inline
fixed_point<T, FRAC_BITS> frac(const fixed_point<T, FRAC_BITS>& a, const fixed_point<T, FRAC_BITS>& b)
{
    return a / b;
}

template<class U, class T, unsigned FRAC_BITS> inline
U fixed_cast(const fixed_point<T, FRAC_BITS>& a)
{
    return static_cast<U>(a.raw()) / static_cast<U>(1<<a.frac_bits);
}

template<class T, unsigned FRAC_BITS> inline
std::ostream& operator<<(std::ostream& os, const fixed_point<T, FRAC_BITS>& a)
{
    double v = fixed_cast<double>(a);
    os << v;
    return os;
}

template<class T, unsigned FRAC_BITS> inline
std::istream& operator>>(std::istream& is, fixed_point<T, FRAC_BITS>& a)
{
    return is;
}

#endif
