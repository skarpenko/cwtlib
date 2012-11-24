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
 * Abstract range functor class
*/
public abstract class RangeFunctor implements Cloneable {
    private String _name;  // Functor name

    /**
     * RangeFunctor constructor
     *
     *  @param Name functor name.
    */
    protected RangeFunctor(String Name)
    {
        _name = Name;
    }

    /**
     * Starting value in range
    */
    public abstract double start();
    /**
     * Ending value in range
    */
    public abstract double end();
    /**
     * Total steps
    */
    public abstract int steps();

    /**
     * Functor name
    */
    public String name()
    {
        return _name;
    }

    /**
     * Used to obtain object clone
    */
    public RangeFunctor clone()
    {
        try {
            return (RangeFunctor)super.clone();
        }
        catch (CloneNotSupportedException ex) {
            // Cloneable interface implemented....
            // ... so that is not possible.
            throw new Error("Unknown Error!");
        }
    }

    /**
     * Evaluates value for specified step
    */
    public abstract double evaluate(int i);
}
