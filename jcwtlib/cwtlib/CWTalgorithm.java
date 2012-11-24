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
 * Collection of CWT algorithms (and CWT helpers)
*/
public final class CWTalgorithm {
    private CWTalgorithm() {}

    /**
     * Returns real part of complex multiplication
    */
    private static double cmplx_mul_re(double r1, double i1, double r2, double i2)
    {
        return r1*r2 - i1*i2;
    }

    /**
     * Returns imaginary part of complex multiplication
    */
    private static double cmplx_mul_im(double r1, double i1, double r2, double i2)
    {
        return r1*i2 + i1*r2;
    }

    /**
     * Compute fast Fourier transform
    */
    private static void fft(double re[], double im[], int n, int off, int isign)
    {
        int  i, j, oi, oj, k, l, le, le1, ip, n2;
        double wpr, wpi, wr, wi, wtr, wti;

        n2 = n>>1;
        j = 1;
        for (i=0; i<n-1; i++) {
            oj = off + j;
            oi = off + i;
            if (i<j) {
                wtr      = re[oj-1];
                wti      = im[oj-1];
                re[oj-1] = re[oi];
                im[oj-1] = im[oi];
                re[oi]   = wtr;
                im[oi]   = wti;
            }
            k = n2;
            while (k<j) {
                j -= k;
                k >>= 1;
            }
            j += k;
        }
        l=1;
        k=n;
        while ((k>>=1) != 0) {
            le1 = (le=1<<l++) >> 1;
            wtr = Math.PI / (double)le1;
            wpr = Math.cos(wtr); wpi = -isign*Math.sin(wtr);
            wr = 1.0;            wi = 0.0;
            for (j=0; j<le1; j++) {
                for (i=j; i<n; i+=le) {
                    oi = off + i;
                    ip = oi + le1;
                    wtr    = wr*re[ip] - wi*im[ip];
                    wti    = wi*re[ip] + wr*im[ip];
                    re[ip] = re[oi] - wtr;
                    im[ip] = im[oi] - wti;
                    re[oi] = re[oi] + wtr;
                    im[oi] = im[oi] + wti;
                }
                wr = (wtr=wr)*wpr - wi*wpi;
                wi = wi*wpr + wtr*wpi;
            }
        }
    }

    /**
     * CWT computation AS IS.
     *
     *  @param s              source signal need to be transformed;
     *  @param Scales         functor which provides scales sequence for
     *                        transform;
     *  @param Translations   functor which provides translations sequence for
     *                        transform;
     *  @param MotherWavelet  mother wavelet used in computations;
     *  @param ivalp          used to add additional time steps when computing
     *                        wavelet (affects precision);
     *  @param Name           name which will be assigned to result object.
     *
     *  @return Returns WTransform obect as a result.
    */
    public static WTransform cwt(Signal s, RangeFunctor Scales,
        RangeFunctor Translations, Wavelet MotherWavelet, int ivalp, String Name)
    {
        // Result
        WTransform wt;
        // references to internal Signal/WTransform data
        double s_re[], s_im[];
        double wt_re[], wt_im[];
        // signal params
        int n = s.length();
        double fs = s.getFs();
        // WT params
        double a, b, T;
        double i, istep;
        // indexes and dimensions
        int dx, dy;
        int rows, cols;
        int row, row_dx;


        // check arguments
        if (Scales.steps() <= 0 || Translations.steps() <= 0 || n <= 0 ||
            fs <= 0.0 || ivalp <= 0)
            throw new IllegalArgumentException("Invalid argument specified!");

        // create result object
        wt = new WTransform(Scales, Translations, MotherWavelet, Name);

        // obtain result dimensions and references to data
        rows = wt.rows();
        cols = wt.cols();
        wt_re = wt.reData();
        wt_im = wt.imData();
        s_re = s.reData();
        s_im = s.imData();

        // index step (used in convolution stage)
        istep = 1.0 / (double)ivalp;

        // Scales
        for (dy = 0; dy < rows; dy++) {
            // obtain current scale
            a = Scales.evaluate(dy) * fs;
            if (a == 0.0) a = Double.MIN_VALUE;

            // set starting index of current row
            row = dy * cols;

            // Translations
            for (dx = 0; dx < cols; dx++) {
                // obtain current translation
                b = Translations.evaluate(dx) * fs;

                // index of convolution result
                row_dx = row + dx;

                // Perform convolution
                wt_re[row_dx] = 0.0;
                wt_im[row_dx] = 0.0;
                for (i = 0.0; i < n; i += istep) {
                    T = (i - b) / a;
                    wt_re[row_dx] += cmplx_mul_re(s_re[(int)i], s_im[(int)i],
                                        MotherWavelet.reT(T), -MotherWavelet.imT(T));
                    wt_im[row_dx] += cmplx_mul_im(s_re[(int)i], s_im[(int)i],
                                        MotherWavelet.reT(T), -MotherWavelet.imT(T));
                    // NOTE: "-" before Wavelet imaginary part indicates complex
                    // conjunction.
                }
                wt_re[row_dx] *= 1.0 / (Math.sqrt(a) * (double)ivalp);
                wt_im[row_dx] *= 1.0 / (Math.sqrt(a) * (double)ivalp);
            }
        }

        return wt;
    }

    /**
     * CWT computation AS IS considering short-term nature of a wavelet
     *
     *  @param s             source signal need to be transformed;
     *  @param Scales        functor which provides scales sequence for
     *                       transform;
     *  @param Translations  functor which provides translations sequence for
     *                       transform;
     *  @param MotherWavelet mother wavelet used in computations;
     *  @param ivalp         used to add additional time steps when computing
     *                       wavelet (affects precision);
     *  @param Name          name which will be assigned to result object.
     *
     *  @return Returns WTransform obect as a result.
    */
    public static WTransform cwto1(Signal s, RangeFunctor Scales,
        RangeFunctor Translations, Wavelet MotherWavelet, int ivalp, String Name)
    {
        // Result
        WTransform wt;
        // references to internal Signal/WTransform data
        double s_re[], s_im[];
        double wt_re[], wt_im[];
        // signal params
        int n = s.length();
        double fs = s.getFs();
        // WT params
        double a, b, T;
        double i, istep;
        // wavelet params
        double t1, t2;
        // indexes and dimensions
        int dx, dy;
        int rows, cols;
        int row, row_dx;
        // precomputed values
        double a_esl, a_esr;


        // check arguments
        if (Scales.steps() == 0 || Translations.steps() == 0 || n == 0 ||
            fs <= 0.0 || ivalp == 0)
            throw new IllegalArgumentException("Invalid argument specified!");

        // create result object
        wt = new WTransform(Scales, Translations, MotherWavelet, Name);

        // obtain result dimensions and references to data
        rows = wt.rows();
        cols = wt.cols();
        wt_re = wt.reData();
        wt_im = wt.imData();
        s_re = s.reData();
        s_im = s.imData();

        // index step (used in convolution stage)
        istep = 1.0 / (double)ivalp;

        // Scales
        for (dy = 0; dy < rows; dy++) {
            // obtain current scale
            a = Scales.evaluate(dy) * fs;
            if (a == 0.0) a = Double.MIN_VALUE;

            // set starting index of current row
            row = dy * cols;

            // obtain wavelet support width on that scale
            a_esl = a*MotherWavelet.effL();  a_esr = a*MotherWavelet.effR();

            // Translations
            for (dx = 0; dx < cols; dx++) {
                // obtain current translation
                b = Translations.evaluate(dx) * fs;

                // index of convolution result
                row_dx = row + dx;

                // compute time range where wavelet presents
                t1 = a_esl + b;          t2 = a_esr + b;
                if (t1 < 0.0) t1 = 0.0;  if (t2 >= n) t2 = (double)(n - 1);

                // Perform convolution
                wt_re[row_dx] = 0.0;
                wt_im[row_dx] = 0.0;
                for (i = t1; i <= t2; i += istep) {
                    T = (i - b) / a;
                    wt_re[row_dx] += cmplx_mul_re(s_re[(int)i], s_im[(int)i],
                                        MotherWavelet.reT(T), -MotherWavelet.imT(T));
                    wt_im[row_dx] += cmplx_mul_im(s_re[(int)i], s_im[(int)i],
                                        MotherWavelet.reT(T), -MotherWavelet.imT(T));
                    // NOTE: "-" before Wavelet imaginary part indicates complex
                    // conjunction.
                }
                wt_re[row_dx] *= 1.0 / (Math.sqrt(a) * (double)ivalp);
                wt_im[row_dx] *= 1.0 / (Math.sqrt(a) * (double)ivalp);
            }
        }

        return wt;
    }

    /**
     * CWT computation AS IS considering short-term nature of a wavelet
     * with its precompution.
     *
     *  @param s             source signal need to be transformed;
     *  @param Scales        functor which provides scales sequence for
     *                       transform;
     *  @param Translations  functor which provides translations sequence for
     *                       transform;
     *  @param MotherWavelet mother wavelet used in computations;
     *  @param ivalp         used to add additional time steps when computing
     *                       wavelet (affects precision);
     *  @param npoints       number of wavelet values for precompution
     *                       (greater value - higher precision);
     *  @param Name          name which will be assigned to result object.
     *
     *  @return Returns WTransform obect as a result.
    */
    public static WTransform cwto2(Signal s, RangeFunctor Scales,
        RangeFunctor Translations, Wavelet MotherWavelet, int ivalp,
        int npoints, String Name)
    {
        // Result
        WTransform wt;
        // references to internal Signal/WTransform data
        double s_re[], s_im[];
        double wt_re[], wt_im[];
        // signal params
        int n = s.length();
        double fs = s.getFs();
        // WT params
        double a, b;
        double i, istep;
        // indexes and dimensions
        int dx, dy;
        int rows, cols;
        int row, row_dx;
        // variables for operating with precomputed wavelet
        double t1, t2;
        double W_re[];
        double W_im[];
        double wstep;
        int L, R, j;
        double d;
        int ind;
        // precomputed values
        double a_esl, a_esr;


        // check arguments
        if (Scales.steps() <= 0 || Translations.steps() <= 0 || n <= 0 ||
            fs <= 0.0 || ivalp <= 0 || npoints <= 0)
            throw new IllegalArgumentException("Invalid argument specified!");

        // create result object
        wt = new WTransform(Scales, Translations, MotherWavelet, Name);

        // obtain result dimensions and references to data
        rows = wt.rows();
        cols = wt.cols();
        wt_re = wt.reData();
        wt_im = wt.imData();
        s_re = s.reData();
        s_im = s.imData();

        //// Precompute wavelet values ////
        // time step for wavelet compution
        wstep = (double)( (MotherWavelet.effR() - MotherWavelet.effL()) / npoints );
        // left and right indexes
        L = (int)Math.floor(MotherWavelet.effL() / wstep);
        R = (int)Math.ceil(MotherWavelet.effR() / wstep);
        // init nr_vectors
        W_re = new double[R - L + 1];
        W_im = new double[R - L + 1];
        // fill vectors with wavelet values
        for (j = 0, i = MotherWavelet.effL(); j <= R-L; j++, i += wstep) {
            W_re[j] = MotherWavelet.reT(i);
            W_im[j] = MotherWavelet.imT(i);
        }
        // scale factor for indexing
        d = (double)( (R - L) / (MotherWavelet.effR() - MotherWavelet.effL()) );

        // index step (used in convolution stage)
        istep = 1.0 / (double)ivalp;

        // Scales
        for (dy = 0; dy < rows; dy++) {
            // obtain current scale
            a = Scales.evaluate(dy) * fs;
            if (a == 0.0) a = Double.MIN_VALUE;

            // set starting index of current row
            row = dy * cols;

            // obtain wavelet support width on that scale
            a_esl = a*MotherWavelet.effL(); a_esr = a*MotherWavelet.effR();

            // Translations
            for (dx = 0; dx < cols; dx++) {
                // obtain current translation
                b = Translations.evaluate(dx) * fs;

                // index of convolution result
                row_dx = row + dx;

                // compute time range where wavelet presents
                t1 = a_esl + b;          t2 = a_esr + b;
                if (t1 < 0.0) t1 = 0.0;  if (t2 >= n) t2 = (double)(n - 1);

                // Perform convolution
                wt_re[row_dx] = 0.0;
                wt_im[row_dx] = 0.0;
                for (i = t1; i <= t2; i += istep) {
                    ind = (int)(d * (i - b) / a) - L;
                    wt_re[row_dx] += cmplx_mul_re(s_re[(int)i], s_im[(int)i],
                                        W_re[ind], -W_im[ind]);
                    wt_im[row_dx] += cmplx_mul_im(s_re[(int)i], s_im[(int)i],
                                        W_re[ind], -W_im[ind]);
                    // NOTE: "-" before Wavelet imaginary part indicates complex
                    // conjunction.
                }
                wt_re[row_dx] *= 1.0 / (Math.sqrt(a) * (double)ivalp);
                wt_im[row_dx] *= 1.0 / (Math.sqrt(a) * (double)ivalp);
            }
        }

        return wt;
    }

    /**
     * CWT computation AS IS considering short-term nature of a wavelet
     * with its precompution and dynamic adjusting ivalp parameter.
     *
     *  @param s             source signal need to be transformed;
     *  @param Scales        functor which provides scales sequence for
     *                       transform;
     *  @param Translations  functor which provides translations sequence for
     *                       transform;
     *  @param MotherWavelet mother wavelet used in computations;
     *  @param ivalp         used to add additional time steps when computing
     *                       wavelet (affects precision);
     *  @param npoints       number of wavelet values for precompution
     *                       (greater value - higher precision);
     *  @param Name          name which will be assigned to result object.
     *
     *  @return Returns WTransform obect as a result.
    */
    public static WTransform cwto3(Signal s, RangeFunctor Scales,
        RangeFunctor Translations, Wavelet MotherWavelet, int ivalp,
        int npoints, String Name)
    {
        // Result
        WTransform wt;
        // references to internal Signal/WTransform data
        double s_re[], s_im[];
        double wt_re[], wt_im[];
        // signal params
        int n = s.length();
        double fs = s.getFs();
        // WT params
        double a, b;
        double i, istep;
        // indexes and dimensions
        int dx, dy;
        int rows, cols;
        int row, row_dx;
        // variables for operating with precomputed wavelet
        double t1, t2;
        double W_re[];
        double W_im[];
        double wstep;
        int L, R, j;
        double d;
        int ind;
        // precomputed values
        double ivalp_amin, a_esl, a_esr;


        // check arguments
        if (Scales.steps() <= 0 || Translations.steps() <= 0 || n <= 0 ||
            fs <= 0.0 || ivalp <= 0 || npoints <= 0)
            throw new IllegalArgumentException("Invalid argument specified!");

        // create result object
        wt = new WTransform(Scales, Translations, MotherWavelet, Name);

        // obtain result dimensions and references to data
        rows = wt.rows();
        cols = wt.cols();
        wt_re = wt.reData();
        wt_im = wt.imData();
        s_re = s.reData();
        s_im = s.imData();

        //// Precompute wavelet values ////
        // time step for wavelet compution
        wstep = (double)( (MotherWavelet.effR() - MotherWavelet.effL()) / npoints );
        // left and right indexes
        L = (int)Math.floor(MotherWavelet.effL() / wstep);
        R = (int)Math.ceil(MotherWavelet.effR() / wstep);
        // init nr_vectors
        W_re = new double[R - L + 1];
        W_im = new double[R - L + 1];
        // fill vectors with wavelet values
        for (j = 0, i = MotherWavelet.effL(); j <= R-L; j++, i += wstep) {
            W_re[j] = MotherWavelet.reT(i);
            W_im[j] = MotherWavelet.imT(i);
        }
        // scale factor for indexing
        d = (double)( (R - L) / (MotherWavelet.effR() - MotherWavelet.effL()) );

        // index step (used in convolution stage)
        istep = 1.0 / (double)ivalp;
        // this used to get default ivalp parameter for minimal possible scale
        ivalp_amin = (double)ivalp * Math.min(Scales.start(), Scales.end());

        // Scales
        for (dy = 0; dy < rows; dy++) {
            // obtain current scale
            a = Scales.evaluate(dy) * fs;
            if (a == 0.0) a = Double.MIN_VALUE;

            // set starting index of current row
            row = dy * cols;

            // compute ivalp and istep for current a
            ivalp = (int)Math.ceil( ivalp_amin / a );
            if (ivalp == 0) ivalp = 1; // ivalp cannot be 0
            istep = 1.0 / (double)ivalp;

            // obtain wavelet support width on that scale
            a_esl = a*MotherWavelet.effL(); a_esr = a*MotherWavelet.effR();

            // Translations
            for (dx = 0; dx < cols; dx++) {
                // obtain current translation
                b = Translations.evaluate(dx) * fs;

                // index of convolution result
                row_dx = row + dx;

                // compute time range where wavelet presents
                t1 = a_esl + b;          t2 = a_esr + b;
                if (t1 < 0.0) t1 = 0.0;  if (t2 >= n) t2 = (double)(n - 1);

                // Perform convolution
                wt_re[row_dx] = 0.0;
                wt_im[row_dx] = 0.0;
                for (i = t1; i <= t2; i += istep) {
                    ind = (int)(d * (i - b) / a) - L;
                    wt_re[row_dx] += cmplx_mul_re(s_re[(int)i], s_im[(int)i],
                                        W_re[ind], -W_im[ind]);
                    wt_im[row_dx] += cmplx_mul_im(s_re[(int)i], s_im[(int)i],
                                        W_re[ind], -W_im[ind]);
                    // NOTE: "-" before Wavelet imaginary part indicates complex
                    // conjunction.
                }
                wt_re[row_dx] *= 1.0 / (Math.sqrt(a) * (double)ivalp);
                wt_im[row_dx] *= 1.0 / (Math.sqrt(a) * (double)ivalp);
            }
        }

        return wt;
    }

    /**
     * CWT computation using FFT (fast Fourier transform).
     *
     *  @param s             source signal need to be transformed;
     *  @param Scales        functor which provides scales sequence for
     *                       transform;
     *  @param MotherWavelet mother wavelet used in computations;
     *  @param Name          name which will be assigned to result object.
     *
     *  @return Returns WTransform obect as a result.
    */
    public static WTransform cwtft(Signal s, RangeFunctor Scales,
        Wavelet MotherWavelet, String Name)
    {
        // Result
        WTransform wt;
        // indexes and dimensions
        int dx, dy;
        int rows, cols;
        // local copy of a source signal
        Signal f = s.clone();
        // signal params
        int n = f.length();
        double fs = f.getFs();
        // references to internal data
        double wt_re[], wt_im[];
        double f_re[], f_im[];
        // variables for wavelet computation
        double a, w, W_re, W_im;
        // precomputed values
        int dy_cols, dy_cols_dx;
        double sqrt_a_n;
        double twoPIn = 2.0 * Math.PI / (double)n;


        // check arguments
        if (Scales.steps() <= 0 || fs <= 0.0)
            throw new IllegalArgumentException("Invalid argument specified!");

        // check that signal length is an integer power of two
        dx = n;
        while (dx != 1) {
            if ( (dx&1)!=0 && dx>1 )
                throw new IllegalArgumentException("Signal length is not an integer power of 2!");
            dx >>= 1;
        }

        // create result object
        wt = new WTransform(Scales, new LinearRangeFunctor(0.0, f.getDt(), f.time()),
                        MotherWavelet, Name);

        // obtain result dimensions and data references
        rows = wt.rows();
        cols = wt.cols();
        wt_re = wt.reData();
        wt_im = wt.imData();
        f_re = f.reData();
        f_im = f.imData();

        // forward Fourier transform of a signal copy
        fft(f_re, f_im, n, 0, 1);

        // Scales
        for (dy = 0; dy < rows; dy++)
        {
            // get current scale
            a = Scales.evaluate(dy) * fs;
            if (a == 0.0) a = Double.MIN_VALUE;

            sqrt_a_n = Math.sqrt(a) / (double)n; // precompution

            // precomputed starting index of result row
            dy_cols = dy * cols;

            // Convolute
            for (dx = 0; dx < cols; dx++)
            {
                // calculate wave number w
                if (dx<=n>>1)
                  w = twoPIn * a * (double)(dx);
                else
                  w =-twoPIn * a * (double)(n-dx);

                // calculate wavelet
                W_re = sqrt_a_n * MotherWavelet.reF(w);
                W_im = sqrt_a_n * MotherWavelet.imF(w);

                // index of result point
                dy_cols_dx = dy_cols + dx;

                // compute result
                wt_re[dy_cols_dx] = f_re[dx] * W_re + f_im[dx] * W_im;
                wt_im[dy_cols_dx] = f_im[dx] * W_re - f_re[dx] * W_im;
            }
            // inverse Fourier transform
            fft(wt_re, wt_im, n, dy * cols, -1);
        }

        return wt;
    }
}
