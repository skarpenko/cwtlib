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

package cwtlib;


/**
 * Signal class
*/
public class Signal {
    private double _re[];       // real part
    private double _im[];       // imaginary part
    private int _length;        // signal length (total samples count)
    private double _fs = 1.0;   // sample frequency
    private String _name;       // signal name

    /**
     * Construct empty Signal object
    */
    public Signal()
    { }

    /**
     * Constructor of Signal object.
     *
     *  @param Length signal length;
     *  @param Real   real part;
     *  @param Imag   imaginary part;
     *  @param Fs     sampling frequency;
     *  @param Name   signal name.
    */
    public Signal(int Length, double Real[], double Imag[], double Fs, String Name)
      throws IllegalArgumentException
    {
        if (Length < 0 || Fs <= 0.0)
            throw new IllegalArgumentException("Invalid argument Length or Fs!");

        _fs = Fs;
        _length = Length;
        _name = Name;

        if (_length == 0)
            return;

        _re = new double[_length];
        _im = new double[_length];

        // copy source data
        if (Real != null)
            System.arraycopy(Real, 0, _re, 0, _length);
        if (Imag != null)
            System.arraycopy(Imag, 0, _im, 0, _length);
    }

    /**
     * Get real value for given index of a sample
     *
     *  @param i index.
    */
    public double re(int i) throws ArrayIndexOutOfBoundsException
    {
        if(i < 0 || i >= _length)
            throw new ArrayIndexOutOfBoundsException("Index out of bounds!");

        return _re[i];
    }

    /**
     * Get real value for given time position of a sample
     *
     *  @param t time.
    */
    public double re(double t) throws ArrayIndexOutOfBoundsException
    {
        int i = (int)(t * _fs);

        if(i < 0 || i >= _length)
            throw new ArrayIndexOutOfBoundsException("Time value is out of range!");

        return _re[i];
    }

    /**
     * Set real value for given index of a sample
     *
     *  @param i index;
     *  @param v value to set.
    */
    public void re(int i, double v) throws ArrayIndexOutOfBoundsException
    {
        if(i < 0 || i >= _length)
            throw new ArrayIndexOutOfBoundsException("Index out of bounds!");

        _re[i] = v;
    }

    /**
     * Set real value for given time position of a sample
     *
     *  @param t time;
     *  @param v value to set.
    */
    public void re(double t, double v) throws ArrayIndexOutOfBoundsException
    {
        int i = (int)(t * _fs);

        if(i < 0 || i >= _length)
            throw new ArrayIndexOutOfBoundsException("Time value is out of range!");

        _re[i] = v;
    }

    /**
     * Get imaginary value for given index of a sample
     *
     *  @param i index.
    */
    public double im(int i) throws ArrayIndexOutOfBoundsException
    {
        if(i < 0 || i >= _length)
            throw new ArrayIndexOutOfBoundsException("Index out of bounds!");

        return _im[i];
    }

    /**
     * Get imaginary value for given time position of a sample
     *
     *  @param t time.
    */
    public double im(double t) throws ArrayIndexOutOfBoundsException
    {
        int i = (int)(t * _fs);

        if(i < 0 || i >= _length)
            throw new ArrayIndexOutOfBoundsException("Time value is out of range!");

        return _im[i];
    }

    /**
     * Set imaginary value for given index of a sample
     *
     *  @param i index;
     *  @param v value to set.
    */
    public void im(int i, double v) throws ArrayIndexOutOfBoundsException
    {
        if(i < 0 || i >= _length)
            throw new ArrayIndexOutOfBoundsException("Index out of bounds!");

        _im[i] = v;
    }

    /**
     * Set imaginary value for given time position of a sample
     *
     *  @param t time;
     *  @param v value to set.
    */
    public void im(double t, double v) throws ArrayIndexOutOfBoundsException
    {
        int i = (int)(t * _fs);

        if(i < 0 || i >= _length)
            throw new ArrayIndexOutOfBoundsException("Time value is out of range!");

        _im[i] = v;
    }

    /**
     * Get magnitude for given index of a sample
     *
     *  @param i index.
    */
    public double mag(int i) throws ArrayIndexOutOfBoundsException
    {
        if(i < 0 || i >= _length)
            throw new ArrayIndexOutOfBoundsException("Index out of bounds!");

        return Math.sqrt(_re[i]*_re[i] + _im[i]*_im[i]);
    }

    /**
     * Get magnitude for given time position of a sample
     *
     *  @param t time.
    */
    public double mag(double t) throws ArrayIndexOutOfBoundsException
    {
        int i = (int)(t * _fs);

        if(i < 0 || i >= _length)
            throw new ArrayIndexOutOfBoundsException("Time value is out of range!");

        return Math.sqrt(_re[i]*_re[i] + _im[i]*_im[i]);
    }

    /**
     * Get angle for given index of a sample
     *
     *  @param i index.
    */
    public double ang(int i) throws ArrayIndexOutOfBoundsException
    {
        if(i < 0 || i >= _length)
            throw new ArrayIndexOutOfBoundsException("Index out of bounds!");

        return Math.atan2(_im[i], _re[i]);
    }

    /**
     * Get angle for given time position of a sample
     *
     *  @param t time.
    */
    public double ang(double t) throws ArrayIndexOutOfBoundsException
    {
        int i = (int)(t * _fs);

        if(i < 0 || i >= _length)
            throw new ArrayIndexOutOfBoundsException("Time value is out of range!");

        return Math.atan2(_im[i], _re[i]);
    }

    /**
     * Signal length in samples
    */
    public int length()
    {
        return _length;
    }

    /**
     * Signal length in time
    */
    public double time()
    {
        return (double)(_length - 1) / _fs;
    }

    /**
     * Get signal name
    */
    public String getName()
    {
        return _name;
    }

    /**
     * Set signal name
     *
     *  @param Name new name.
    */
    public void setName(String Name)
    {
        _name = Name;
    }

    /**
     * Get signal sampling frequency
    */
    public double getFs()
    {
        return _fs;
    }

    /**
     * Set signal sampling frequency
     *
     *  @param Fs new sampling frequency.
    */
    public void setFs(double Fs) throws IllegalArgumentException
    {
        if (Fs <= 0.0)
            throw new IllegalArgumentException("Invalid Fs value!");

        _fs = Fs;
    }

    /**
     * Get sampling time step value
     * (the same as sampling frequency)
    */
    public double getDt()
    {
        return 1.0 / _fs;
    }

    /**
     * Set sampling time step value
     * (the same as sampling frequency)
     *
     *  @param Dt new time step.
    */
    public void setDt(double Dt) throws IllegalArgumentException
    {
        if (Dt <= 0.0)
            throw new IllegalArgumentException("Invalid Dt value!");

        _fs = 1.0 / Dt;
    }

    /**
     * Assign another signal to this one
     *
     *  @param Src source signal.
    */
    public void assign(Signal Src)
    {
        if (this == Src)
            return;

        // assign source fields
        _name = Src._name;
        _length = Src._length;
        _fs = Src._fs;
        _re = _im = null;

        if (_length == 0)
            return;

        // allocate storage and copy source data
        _re = new double[_length];
        _im = new double[_length];

        System.arraycopy(Src._re, 0, _re, 0, _length);
        System.arraycopy(Src._im, 0, _im, 0, _length);
    }

    /**
     * Set new signal size.
     *
     * If new size is less than current size, signal will be truncated to new size.
     * If new size is greater than current size, signal will be padded with zeros.
     *
     *  @param NewSize new size of a signal.
    */
    public void resize(int NewSize) throws IllegalArgumentException
    {
        if (NewSize < 0)
            throw new IllegalArgumentException("Invalid NewSize value!");

        // if resizing is not needed
        if (_length == NewSize)
            return;

        // delete internal data if NewSize is zero
        if (NewSize == 0) {
            _length = 0;
            _re = _im = null;
        } else {
            double tmp_re[];
            double tmp_im[];
            // data size need to be copied from old storage
            int cpy_size = (NewSize < _length) ? NewSize : _length;

            // allocate storage and copy new data
            tmp_re = new double[NewSize];
            tmp_im = new double[NewSize];

            System.arraycopy(_re, 0, tmp_re, 0, cpy_size);
            System.arraycopy(_im, 0, tmp_im, 0, cpy_size);

            // set new fields
            _re = tmp_re;
            _im = tmp_im;
            _length = NewSize;
        }
    }

    /**
     * Used to obtain object clone
    */
    public Signal clone()
    {
        return new Signal(_length, _re, _im, _fs, _name);
    }

    /**
     * Get reference to internal real data
    */
    public double[] reData()
    {
        return _re;
    }

    /**
     * Get reference to internal imaginary data
    */
    public double[] imData()
    {
        return _im;
    }
}
