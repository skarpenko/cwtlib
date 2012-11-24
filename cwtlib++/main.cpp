/*
 *   main.cpp - Continuous Wavelet Transform Demo Program
 *
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

#include <ctime>
#include <iostream>
#include <fstream>
#include <stdexcept>
using namespace std;
#include <cwtlib>
using namespace cwtlib;

enum { Real, Imag, Magn, Angl };  // WT transform parts

// Defines CWT algorithm version to use.
// Valid versions is:
//   0 - Not optimized version;
//   1 - Optimized version 1;
//   2 - Optimized version 2;
//   3 - Optimized version 3;
//   4 - FFT based version.
#define CWTVER 3

// WT computation params
#define WAVELET ComplexMorlet()  /* Mother wavelet */
#define PART Magn                /* Complex part of the transform to save */
#define AMIN 1                   /* A min  */
#define ASTP 1                   /* A step */
#define AMAX 128                 /* A max  */
#define PREC 4                   /* Precision */
#define NPOINTS 60000            /* Points number for wavelet precompution */
#define N 1024                   /* Signal length for reading from a file */

// used for time tracking
clock_t msecs;

// Main()
int main(int argc, char *argv[])
{
    fstream file;
    cwt_uint_t i, j;
    // initial signal object
    Signal s(N, NULL, NULL, 1.0, "Test_Signal");

    // check command line params
    if(argc<3) {
        cout << "filename missing\nUsage: "
             << argv[0]
             << " input_file output_file\n";
        return 1;
    }

    // open input file
    file.open(argv[1], ios::in);
    if(!file) {
        cout << "can't open file " << argv[1] << endl;
        return 1;
    }

    // obtain DataProxy for direct access to signal data
    ReDataProxy data = s.reData();

    // read data from file
    for (i = 0; i < N && file; i++) {
        file >> data[i];
    }
    file.close(); // and close after read

    // open output file
    file.open(argv[2], ios::out);
    if(!file) {
        cout << "can't open file " << argv[2] << endl;
        return 1;
    }

    // create special object for wavelet transform
    // that is scales and translations functors
    LinearRangeFunctor Scales(AMIN, ASTP, AMAX);
    LinearRangeFunctor Translations(0.0, s.getDt(), s.time());
    // WT result goes here
    WTransform *wt;

    // ok. lets start.
    cout << "Performing wavelet transform... "; cout.flush();
    msecs = clock();
    try {
#if (CWTVER == 0)
        wt = CWTalgorithm::cwt(s, Scales, Translations, WAVELET, PREC, "WT_test");
#endif
#if (CWTVER == 1)
        wt = CWTalgorithm::cwto1(s, Scales, Translations, WAVELET, PREC, "WT_test");
#endif
#if (CWTVER == 2)
        wt = CWTalgorithm::cwto2(s, Scales, Translations, WAVELET, PREC, NPOINTS, "WT_test");
#endif
#if (CWTVER == 3)
        wt = CWTalgorithm::cwto3(s, Scales, Translations, WAVELET, PREC, NPOINTS, "WT_test");
#endif
#if (CWTVER == 4)
        wt = CWTalgorithm::cwtft(s, Scales, WAVELET, "WT_Test");
#endif
    }
    catch (exception& e) {
        cout << "error performing wavelet transform!\n"
             << "Exception message: " << e.what() << endl;
        return 1;
    }
    cout << "Elapsed time: " << ((clock()-msecs) / CLOCKS_PER_SEC) << " sec.\n";

    // print transform details
    cout << "\nTransform details:\n"
         << "  WT result name.............. " << wt->getName() << endl
         << "  Scales functor name......... " << wt->scales().name() << endl
         << "  Amin........................ " << wt->scales().start() << endl
         << "  Amax........................ " << wt->scales().end() << endl
         << "  Translations functor name... " << wt->translations().name() << endl
         << "  Bmin........................ " << wt->translations().start() << endl
         << "  Bmax........................ " << wt->translations().end() << endl
         << "  Mother wavelet name......... " << wt->motherWavelet().name() << endl
         << "  CWT result dimensions....... " << wt->rows() << " x " << wt->cols() << endl
         << "  Signal name................. " << s.getName() << endl
         << "  Signal length............... " << s.length() << " samples\n"
         << "  Signal sampling frequency... " << s.getFs() << " Hz\n"
#if (CWTVER == 0)
         << "  Algorithm used.............. cwt()\n";
#endif
#if (CWTVER == 1)
         << "  Algorithm used.............. cwto1()\n";
#endif
#if (CWTVER == 2)
         << "  Algorithm used.............. cwto2()\n";
#endif
#if (CWTVER == 3)
         << "  Algorithm used.............. cwto3()\n";
#endif
#if (CWTVER == 4)
         << "  Algorithm used.............. cwtft()\n";
#endif

    // store result into file
    cout << "\nSaving ";
    switch (PART) { 
        case Real:  cout << "real part";       break;
        case Imag:  cout << "imaginary part";  break;
        case Magn:  cout << "magnitude";       break;
        case Angl:  cout << "angle";           break;
        default: throw "invalid parameter!";
    }
    cout << " of the transform...";

    for(i = 0; i < wt->rows(); i++) {
        for(j = 0; j < wt->cols(); j++) {
            switch (PART) { 
                case Real:  file << " " << wt->re(i, j);   break;
                case Imag:  file << " " << wt->im(i, j);   break;
                case Magn:  file << " " << wt->mag(i, j);  break;
                case Angl:  file << " " << wt->ang(i, j);  break;
            }
        }
        file << endl;
    }
    file.close();

    cout << " done.\n";

    // all done.
    return 0;
}
