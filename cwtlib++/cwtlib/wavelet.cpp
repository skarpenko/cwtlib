/*
 *   Continuous Wavelet Transform Library
 *   Copyright (C) 2004-2009 Stepan V. Karpenko <carp@mail.ru>
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

#include <cmath>
#include <stdexcept>
using namespace std;
#include <wavelet>


__CWTLIB_BEGIN_NAMESPACE


// ==================== Abstract wavelet class ================================

Wavelet::Wavelet(const string& Name)
{
    _name = Name;
}

Wavelet::Wavelet(const Wavelet& Src)
{
    _name = Src._name;
}

Wavelet::~Wavelet()
{}

const string& Wavelet::name() const
{
    return _name;
}

// ==================== Complex Morlet wavelet ================================

// default params
static const char         CMorlet_Name[] = "ComplexMorlet";
static const cwt_float_t  CMorlet_Fc     = 0.8;
static const cwt_float_t  CMorlet_Fb     = 2.0;

CMorlet::CMorlet()
  : Wavelet(CMorlet_Name)
{
    _fc = CMorlet_Fc;
    _fb = CMorlet_Fb;
    // compute L2 norm ...
    _c = 1.0 / sqrt(PI * _fb);
    // ... and effective support boundary values
    _effl = -2.0*_fb;
    _effr = +2.0*_fb;
}

CMorlet::CMorlet(cwt_float_t Fc, cwt_float_t Fb)
  : Wavelet(CMorlet_Name)
{
    if (Fc <= 0.0 || Fb <= 0.0)
        throw CWTLIB_EXCEPTION_INVALID_ARG();

    _fc = Fc;
    _fb = Fb;
    _c = 1.0 / sqrt(PI * _fb);
    _effl = -2.0*_fb;
    _effr = +2.0*_fb;
}

CMorlet::CMorlet(const CMorlet& Src)
  : Wavelet(Src)
{
    _fc = Src._fc;
    _fb = Src._fb;
    _c = Src._c;
    _effl = Src._effl;
    _effr = Src._effr;
}

cwt_float_t CMorlet::reT(cwt_float_t t) const
{
    return _c * exp(-(t*t) / _fb) * cos(PI2 * _fc * t);
}

cwt_float_t CMorlet::imT(cwt_float_t t) const
{
    return _c * exp(-(t*t) / _fb) * sin(PI2 * _fc * t);
}

cwt_float_t CMorlet::reF(cwt_float_t w) const
{
    cwt_float_t br;

    br = (w - PI2 * _fc);
    return exp(-_fb * br * br / 4.0);
}

cwt_float_t CMorlet::imF(cwt_float_t w) const
{
    return 0.0;
}

cwt_float_t CMorlet::fBand() const
{
    return _fb;
}

cwt_float_t CMorlet::cFreq() const
{
    return _fc;
}

cwt_float_t CMorlet::effL() const
{
    return _effl;
}

cwt_float_t CMorlet::effR() const
{
    return _effr;
}

Wavelet* CMorlet::clone() const
{
    return new CMorlet(*this);
}


__CWTLIB_END_NAMESPACE
