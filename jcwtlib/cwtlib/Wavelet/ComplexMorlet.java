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

package cwtlib.Wavelet;


/**
 * Complex Morlet wavelet
*/
public class ComplexMorlet extends Wavelet {
    /**
     * Default wavelet name
    */
    public final static String NAME = "ComplexMorlet";
    /**
     * Default central frequency value
    */
    public final static double FC = 0.8;
    /**
     * Default bandwidth parameter
    */
    public final static double FB = 2.0;

    private double _fc;    // central frequency
    private double _fb;    // bandwidth parameter
    private double _c;     // L2 norm
    // effective support params
    private double _effl;
    private double _effr;

    /**
     * Construct Complex Morlet wavelet with default
     * parameters.
    */
    public ComplexMorlet()
    {
        super(NAME);

        _fc = FC;
        _fb = FB;
        // compute L2 norm ...
        _c = 1.0 / Math.sqrt(Math.PI * _fb);
        // ... and effective support boundary values
        _effl = -2.0*_fb;
        _effr = +2.0*_fb;
    }

    /**
     * Construct Complex Morlet wavelet with user-provided
     * parameters.
     *
     *  @param Fc central frequency;
     *  @param Fb bandwidth parameter.
    */
    public ComplexMorlet(double Fc, double Fb)
    {
        super(NAME);
        if (Fc <= 0.0 || Fb <= 0.0)
            throw new IllegalArgumentException("Invalid parameter passed!");

        _fc = Fc;
        _fb = Fb;
        _c = 1.0 / Math.sqrt(Math.PI * _fb);
        _effl = -2.0*_fb;
        _effr = +2.0*_fb;
    }

    public double reT(double t)
    {
        return _c * Math.exp(-(t*t) / _fb) * Math.cos(2.0 * Math.PI * _fc * t);
    }

    public double imT(double t)
    {
        return _c * Math.exp(-(t*t) / _fb) * Math.sin(2.0 * Math.PI * _fc * t);
    }

    public double reF(double w)
    {
        double br;

        br = (w - 2.0 * Math.PI * _fc);

        return Math.exp(-_fb * br * br / 4.0);
    }

    public double imF(double w)
    {
        return 0.0;
    }

    /**
     * Returns bandwidth parameter
    */
    public double fBand()
    {
        return _fb;
    }

    public double cFreq()
    {
        return _fc;
    }

    public double effL()
    {
        return _effl;
    }

    public double effR()
    {
        return _effr;
    }
}
