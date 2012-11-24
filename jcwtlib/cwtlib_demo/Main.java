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

package cwtlib_demo;

import java.io.*;
import cwtlib.Wavelet.*;
import cwtlib.RangeFunctor.*;
import cwtlib.*;


// WT transform parts
enum wtResultType {
    Real,  /* Real part      */
    Imag,  /* Imaginary part */
    Magn,  /* Magnitude      */
    Angl   /* Angle (Theta)  */
}

public class Main {
    /**
     * Defines CWT algorithm version to use.
     * Valid versions is:
     *   0 - Not optimized version;
     *   1 - Optimized version 1;
     *   2 - Optimized version 2;
     *   3 - Optimized version 3;
     *   4 - FFT based version.
    */
    public static final int CWTVER = 3;

    /* WT computation params */

    /**
     * Mother wavelet
    */
    public static final Wavelet WAVELET = new ComplexMorlet();
    /**
     * Complex part of the transform to save
    */
    public static final wtResultType PART = wtResultType.Magn;
    /**
     * A min (starting scale)
    */
    public static final double AMIN = 1.0;
    /**
     * A step (scale step)
    */
    public static final double ASTP = 1.0;
    /**
     * A max (ending scale)
    */
    public static final double AMAX = 128.0;
    /**
     * Precision
    */
    public static final int PREC = 4;
    /**
     * Points number for wavelet precompution
    */
    public static final int NPOINTS = 60000;
    /**
     * Signal length for reading from a file
    */
    public static final int N = 1024;

    // Main()
    public static void main(String args[])
    {
        double data[];
        long msecs;

        // check command line params
        if (args.length < 2) {
            System.out.println("filename missing");
            System.out.println("Usage: java -jar cwtlib_demo.jar input_file output_file");
            return;
        }

        // initial signal object
        Signal s = new Signal(N, null, null, 1.0, "Test_Signal");

        // load signal data from file
        try {
            loadSignal(args[0], s);
        }
        catch (IOException ex) {
            System.out.println("failed to load data from file " + args[0]);
            return;
        }

        // create special object for wavelet transform
        // that is scales and translations functors
        LinearRangeFunctor Scales = new LinearRangeFunctor(AMIN, ASTP, AMAX);
        LinearRangeFunctor Translations = new LinearRangeFunctor(0.0, s.getDt(), s.time());
        // WT result goes here
        WTransform wt;

        // ok. lets start.
        System.out.print("Performing wavelet transform... ");

        msecs = System.currentTimeMillis();

        try {
            switch(CWTVER) {
                case 0:
                  wt = CWTalgorithm.cwt(s, Scales, Translations, WAVELET, PREC, "WT_test");
                  break;
                case 1:
                  wt = CWTalgorithm.cwto1(s, Scales, Translations, WAVELET, PREC, "WT_test");
                  break;
                case 2:
                  wt = CWTalgorithm.cwto2(s, Scales, Translations, WAVELET, PREC, NPOINTS, "WT_test");
                  break;
                case 3:
                  wt = CWTalgorithm.cwto3(s, Scales, Translations, WAVELET, PREC, NPOINTS, "WT_test");
                  break;
                case 4:
                  wt = CWTalgorithm.cwtft(s, Scales, WAVELET, "WT_Test");
                  break;
                default:
                  throw new Error("CWTVER is invalid!");
            }
        }
        catch(Exception e) {
            System.out.println("error performing wavelet transform!");
            System.out.println("Exception message: " + e.getMessage());
            return;
        }

        System.out.println("Elapsed time: " + 
            Double.toString( (System.currentTimeMillis() - msecs) / 1000.0 ) +
            " sec.");

        // print transform details
        System.out.println("\nTransform details:\n");
        System.out.println("  WT result name.............. " + wt.getName());
        System.out.println("  Scales functor name......... " + wt.scales().name());
        System.out.println("  Amin........................ " + Double.toString(wt.scales().start()));
        System.out.println("  Amax........................ " + Double.toString(wt.scales().end()));
        System.out.println("  Translations functor name... " + wt.translations().name());
        System.out.println("  Bmin........................ " + Double.toString(wt.translations().start()));
        System.out.println("  Bmax........................ " + Double.toString(wt.translations().end()));
        System.out.println("  Mother wavelet name......... " + wt.motherWavelet().name());
        System.out.println("  CWT result dimensions....... " + Integer.toString(wt.rows())
                                                             + " x " + Integer.toString(wt.cols()));
        System.out.println("  Signal name................. " + s.getName());
        System.out.println("  Signal length............... " + Integer.toString(s.length()) + " samples");
        System.out.println("  Signal sampling frequency... " + Double.toString(s.getFs()) + " Hz");
        switch(CWTVER) {
            case 0:
              System.out.println("  Algorithm used.............. cwt()");
              break;
            case 1:
              System.out.println("  Algorithm used.............. cwto1()");
              break;
            case 2:
              System.out.println("  Algorithm used.............. cwto2()");
              break;
            case 3:
              System.out.println("  Algorithm used.............. cwto3()");
              break;
            case 4:
              System.out.println("  Algorithm used.............. cwtft()");
              break;
            default:
              throw new Error("CWTVER is invalid!");
        }      

        // store result into file
        System.out.print("\nSaving ");
        switch(PART) {
          case Real: System.out.print("real part");      break;
          case Imag: System.out.print("imaginary part"); break;
          case Magn: System.out.print("magnitude");      break;
          case Angl: System.out.print("angle");          break;
          default:
              throw new Error("invalid parameter!");
        }
        System.out.print(" of the transform...");

        try {
            saveTransform(args[1], wt, PART);
        }
        catch (IOException ex)
        {
            System.out.println("failed to save data to file " + args[1]);
            return;
        }

        // all done.
        System.out.println(" done.");
    }

    /**
     * Loads signal from file
    */
    public static void loadSignal(String filename, Signal s) throws IOException
    {
        Reader rd = new BufferedReader(new FileReader(filename));
        StreamTokenizer stok = new StreamTokenizer(rd);
        double data[] = s.reData();
        int len = s.length();
        int i = 0;

        stok.parseNumbers();
        stok.nextToken();
        while (stok.ttype != StreamTokenizer.TT_EOF && i < len) {
            if (stok.ttype == StreamTokenizer.TT_NUMBER)
                data[i] = stok.nval;
            else
                break;
            i++;
            stok.nextToken();
        }
        rd.close();
    }

    /**
     * Stores WTransform data into file
    */
    public static void saveTransform(String filename, WTransform wt,
        wtResultType type) throws IOException
    {
        FileWriter fw = new FileWriter(filename);
        int i, j;

        for (j = 0; j < wt.rows(); j++) {
            for (i = 0; i < wt.cols(); i++ )
            {
                switch(PART) {
                  case Real: 
                    fw.write( Double.toString(wt.re(j, i)) + " " );
                    break;
                  case Imag:
                    fw.write( Double.toString(wt.im(j, i)) + " " );
                    break;
                  case Magn:
                    fw.write( Double.toString(wt.mag(j, i)) + " " );
                    break;
                  case Angl:
                    fw.write( Double.toString(wt.ang(j, i)) + " " );
                    break;
                  default:
                    throw new Error("invalid parameter!");
                }
            }
            fw.write("\n");
        }
        fw.close();
    }
}
