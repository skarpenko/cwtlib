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


import cwtlib.Wavelet.*;
import cwtlib.RangeFunctor.*;


/**
 * CWT result class
*/
public class WTransform {
    private double _re[];                // real part of transform
    private double _im[];                // imaginary part of transform
    private int  _cols;                  // columns number
    private int _rows;                   // rows number
    private String _name;                // result name
    private Wavelet _wavelet;            // mother wavelet used for transform
    private RangeFunctor _scales;        // scales functor used
    private RangeFunctor _translations;  // translations functor used

    /**
     * Constructor of WTransform object (in most cases used in CWT algorithms,
     * but can be constructed manually).
     *
     *  @param Scales        functor which provides scales sequence;
     *  @param Translations  functor which provides translations sequence;
     *  @param MotherWavelet mother wavelet;
     *  @param Name          object name.
    */
    public WTransform(RangeFunctor Scales, RangeFunctor Translations,
                       Wavelet MotherWavelet, String Name)
      throws IllegalArgumentException
    {
        _name = Name;
        _rows = Scales.steps();
        _cols = Translations.steps();

        // transform result cannot be empty
        if (_rows <= 0 || _cols <= 0)
            throw new IllegalArgumentException("Invalid dimensions provided!");

        // copy necessary objects
        _scales = Scales.clone();
        _translations = Translations.clone();
        _wavelet = MotherWavelet.clone();

        // allocate storage
        _re = new double[ _rows * _cols ];
        _im = new double[ _rows * _cols ];
    }

    /**
     * Get real value for given row and column.
    */
    public double re(int row, int col) throws ArrayIndexOutOfBoundsException
    {
        if (row < 0 || row >= _rows || col < 0 || col >= _cols)
            throw new ArrayIndexOutOfBoundsException("Dimensions out of bounds!");

        return _re[row*_cols + col];
    }

    /**
     * Set real value for given row and column.
    */
    public void re(int row, int col, double v) throws ArrayIndexOutOfBoundsException
    {
        if (row < 0 || row >= _rows || col < 0 || col >= _cols)
            throw new ArrayIndexOutOfBoundsException("Dimensions out of bounds!");

        _re[row*_cols + col] = v;
    }

    /**
     * Get imaginary value for given row and column.
    */
    public double im(int row, int col) throws ArrayIndexOutOfBoundsException
    {
        if (row < 0 || row >= _rows || col < 0 || col >= _cols)
            throw new ArrayIndexOutOfBoundsException("Dimensions out of bounds!");

        return _im[row*_cols + col];
    }

    /**
     * Set imaginary value for given row and column.
    */
    public void im(int row, int col, double v) throws ArrayIndexOutOfBoundsException
    {
        if (row < 0 || row >= _rows || col < 0 || col >= _cols)
            throw new ArrayIndexOutOfBoundsException("Dimensions out of bounds!");

        _im[row*_cols + col] = v;
    }

    /**
     * Get magnitude for given row and column.
    */
    public double mag(int row, int col) throws ArrayIndexOutOfBoundsException
    {
        if (row < 0 || row >= _rows || col < 0 || col >= _cols)
            throw new ArrayIndexOutOfBoundsException("Dimensions out of bounds!");

        int idx = row*_cols + col;
        return Math.sqrt(_re[idx]*_re[idx] + _im[idx]*_im[idx]);
    }

    /**
     * Get angle for given row and column.
    */
    public double ang(int row, int col) throws ArrayIndexOutOfBoundsException
    {
        if (row < 0 || row >= _rows || col < 0 || col >= _cols)
            throw new ArrayIndexOutOfBoundsException("Dimensions out of bounds!");

        int idx = row*_cols + col;
        return Math.atan2(_im[idx], _re[idx]);
    }

    /**
     * Get rows count
    */
    public int rows()
    {
        return _rows;
    }

    /**
     * Get columns count
    */
    public int cols()
    {
        return _cols;
    }

    /**
     * Get MotherWavelet
    */
    public Wavelet motherWavelet()
    {
        return _wavelet;
    }

    /**
     * Get scales range
    */
    public RangeFunctor scales()
    {
        return _scales;
    }

    /**
     * Get translations range
    */
    public RangeFunctor translations()
    {
        return _translations;
    }

    /**
     * Get result name
    */
    public String getName()
    {
        return _name;
    }

    /**
     * Set new result name
    */
    public void setName(String Name)
    {
        _name = Name;
    }

    /**
     * Assign another result object to this one
    */
    public void assign(WTransform Src)
    {
        double tmp_re[], tmp_im[];
        RangeFunctor tmp_scales, tmp_translations;
        Wavelet tmp_wavelet;
        int tmp_cols, tmp_rows;

        if (Src == this)
            return;

        tmp_rows = Src.scales().steps();
        tmp_cols = Src.translations().steps();

        // check new dimensions
        if (tmp_rows <= 0 || tmp_cols <= 0)
            throw new IllegalArgumentException("Invalid dimensions provided!");

        // copy necessary objects
        tmp_scales = Src.scales().clone();
        tmp_translations = Src.translations().clone();
        tmp_wavelet = Src.motherWavelet().clone();

        // allocate temporary storage
        tmp_re = new double[ tmp_rows * tmp_cols ];
        tmp_im = new double[ tmp_rows * tmp_cols ];

        // copy data
        System.arraycopy(Src._re, 0, tmp_re, 0, tmp_rows * tmp_cols);
        System.arraycopy(Src._im, 0, tmp_im, 0, tmp_rows * tmp_cols);

        // assign new data
        _scales = tmp_scales;

        // translations
        _translations = tmp_translations;

        // wavelet
        _wavelet = tmp_wavelet;

        // translation data
        _re = tmp_re;
        _im = tmp_im;

        // dimensions
        _rows = tmp_rows;
        _cols = tmp_cols;

        // name
        _name = Src._name;
    }

    /**
     * Used to obtain object clone
    */
    public WTransform clone()
    {
        return new WTransform(_scales, _translations, _wavelet, _name);
    }

    /**
     * Get data proxy for direct acces to internal real data
    */
    public double[] reData()
    {
        return _re;
    }

    /**
     * Get data proxy for direct acces to internal imaginary data
    */
    public double[] imData()
    {
        return _im;
    }
}
