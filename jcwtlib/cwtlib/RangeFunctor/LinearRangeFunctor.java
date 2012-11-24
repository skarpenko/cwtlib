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

package cwtlib.RangeFunctor;


/**
 * Linear range functor class
 * (provides linear range of values)
*/
public class LinearRangeFunctor extends RangeFunctor {
    /**
     * Default functor name
    */
    public final static String NAME = "LinearRange";

    private double _start; // starting value
    private double _step;  // step
    private double _end;   // ending value
    private int _steps;    // total steps count

    /**
     * LinearRangeFunctor constructor.
     * 
     *  @param Start starting value;
     *  @param Step  step size;
     *  @param End   ending value (this will be recomputed inside to provide
     *               ectual ending value according to given step size).
    */ 
    public LinearRangeFunctor(double Start, double Step, double End)
      throws IllegalArgumentException
    {
        super(NAME);

        if (Start > End || Step <= 0.0)
            throw new IllegalArgumentException("Invalid range passed!");

        _steps = (int)Math.floor( (End - Start) / Step ) + 1;
        _start = Start;
        _step = Step;
        // obtain actual ending value according to provided step size
        _end = _start + (_steps - 1) * _step;
    }

    public double start()
    {
        return _start;
    }

    public double step()
    {
        return _step;
    }

    public double end()
    {
        return _end;
    }

    public int steps()
    {
        return _steps;
    }

    public double evaluate(int i) throws ArrayIndexOutOfBoundsException
    {
        if(i < 0 || i >= _steps)
            throw new ArrayIndexOutOfBoundsException("Index out if bounds!");

        return _start + i * _step;
    }
}
