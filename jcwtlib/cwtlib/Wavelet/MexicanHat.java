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
 * Mexican Hat wavelet
*/
public class MexicanHat extends Wavelet {
    /**
     * Default wavelet name
    */
    public final static String NAME = "MexicanHat";
    /**
     * Default central frequency value
    */
    public final static double FC = 1.0 / Math.PI;
    /**
     * Default L2 norm. c = 2 / ( sqrt(3) * pi^(1/4) )
    */
    public final static double C = 0.8673250705840776;
    /**
     * Default radius
    */
    public final static double R = 5.0;

    /**
     * Construct Mexican Hat wavelet
    */
    public MexicanHat()
    {
        super(NAME);
    }

    public double reT(double t)
    {
        t = t * t;
        return C * (1.0 - t) * Math.exp(-t / 2.0);
    }

    public double imT(double t)
    {
        return 0.0;
    }

    public double reF(double w)
    {
        w = w * w;
        return C * Math.sqrt(2.0 * Math.PI) * w * Math.exp(-w / 2.0);
    }

    public double imF(double w)
    {
        return 0.0;
    }

    public double cFreq()
    {
        return FC;
    }

    public double effL()
    {
        return -R;
    }

    public double effR()
    {
        return +R;
    }
}
