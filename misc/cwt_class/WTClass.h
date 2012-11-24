//////////////////////////////////////////////////////////////////////////
//                  Continuous Wavelet Transform Class
//
// Author   : Stepan V. Karpenko (carp@mail.ru)
// Date     : 06-09-2004
// Comments :
// History  :
//
//////////////////////////////////////////////////////////////////////////

#ifndef WTClass_H
#define WTClass_H

#include <math.h>
#include <vector>
#include <fstream>
#include "cwtlib/cwtwlets.h"
#include "cwtlib/cwt.h"
using namespace std;

class WTClass {
  private:
//=== Continuous wavelet transform parameters ===
    unsigned long Wavelet;
    double Amin, Astep, Amax;
    double Bstep;
    unsigned long Precision;
    unsigned long WavPoints;
    long WavPart;
    bool FixEdges;

//=== Flags ===
    bool series_flag, cwt_flag, globale_flag, locale_flag, norm_flag,
         localesk_flag, cwtsk_flag, globalem_flag;

    bool cwtsk_loc;         // CWT Skeleton Localization
                            // true - vertical, false - horizontal
    bool localesk_loc;      // Local Energy Function Skeleton Localization
                            // true - vertical, false - horizontal
  
//=== Private routines  ===
    // Class initialization routine
    void ClassInit();

    // Fixes CWT boundaries
    bool FixCWTBounds();

    // Calculates skeleton
    // w - source, loc - maxima localization
    // loc = true - vertical
    // loc = false - horizontal
    bool Skeleton(cwt_t *w, bool loc);

  public:
    vector<double> series;       // Time series for analysis
    vector<double> scales;       // Calculated scales
    vector<double> freqs;        // Calculated frequencies
    vector<double> transls;      // Wavelet translations
    vector<double> globale;      // Global Energy Function
    vector<double> globale_maxs; // Global Energy Function maximas
    cwt_t locale;                // Local Energy Function
    cwt_t locale_skel;           // Local Energy Function Skeleton
    cwt_t cwt;                   // Wavelet Transform
    cwt_t cwt_skel;              // Wavelet Transform Skeleton

    // Constructor
    WTClass();

    // Constructor with paramteters
    // dat - source series, n - data length
    WTClass(double *dat, long n);

    // Destructor
    ~WTClass();

    // Sets Wavelet Transform parameters
    // for more details see cwtlib
    void set_wt_params(unsigned long Wavelet, double Amin, double Astep,
                       double Amax, double Bstep, unsigned long Precision,
                       unsigned long WavPoints, long WavPart, bool FixEdges);

    // Add time series for analysis
    // dat - source series, n - data length
    void series_add(double *dat, long n);

    // Performs Continuous Wavelet Transform
    bool do_cwt();

    // Calculates CWT Skeleton
    // loc - localization
    bool cwt_skelet(bool loc);

    // Normalizes amplitudes of wavelet coefficients
    // NOTE: This method clears globale, locale, locale_skel and cwt_skel
    void norm_cwt();

    // Restores amplitudes of wavelet coefficients normalized by norm_cwt()
    // NOTE: This method clears globale, locale, locale_skel and cwt_skel
    void denorm_cwt();

    // Calculates Global Energy Function
    // NOTE: This method clears globale_maxs
    void global_energy();

    // Calculates Global Energy Function maximas
    bool global_energy_maxs();

    // Calculates Local Energy Function
    // NOTE: This method clears locale_skel
    bool local_energy();

    // Calculates Local Energy Function Skeleton
    // loc - localization
    bool local_energy_skelet(bool loc);

    // Returns row index associated with given frequency.
    // rounding controls method's behavior when given
    // frequency doesn't exist in freqs[] list.
    // =true, frequency will be rounded up
    // =false, frequency will be rounded down
    unsigned long freq2idx(double freq, bool rounding);

    // Returns row index associated with given scale.
    // rounding controls method's behavior when given
    // scale doesn't exist in scales[] list.
    // =true, scale will be rounded up
    // =false, scale will be rounded down
    unsigned long scale2idx(double scale, bool rounding);

    // Returns row index associated with given wavelet translation.
    // rounding controls method's behavior when given
    // translation doesn't exist in transls[] list.
    // =true, translation will be rounded up
    // =false, translation will be rounded down
    unsigned long transl2idx(double transl, bool rounding);

    // Returns list of index pairs [row,col] of maximal values of CWT
    // between low_fr and hig_fr frequencies.
    // rnlf and rnhf rounding parameter for low_fr and hig_fr
    // (see freq2idx for more details).
    // Be careful it doesn't check input params, so the
    // interval [low_fr..hig_fr] must be wide enought.
    // Output vector size is zero if error occured.
    vector<unsigned long> cwt_mvalsf(double low_fr, double hig_fr,
                                     bool rnlf, bool rnhf);

    // Returns list of index pairs [row,col] of maximal values of CWT
    // between low_sc and hig_sc scales.
    // rnls and rnhs rounding parameter for low_sc and hig_sc
    // (see scale2idx for more details).
    // Be careful it doesn't check input params, so the
    // interval [low_sc..hig_sc] must be wide enought.
    // Output vector size is zero if error occured.
    vector<unsigned long> cwt_mvalss(double low_sc, double hig_sc,
                                     bool rnls, bool rnhs);

    // Returns list of index pairs [row,col] of maximal values of Local
    // Energy Function between low_fr and hig_fr frequencies.
    // rnlf and rnhf rounding parameter for low_fr and hig_fr
    // (see freq2idx for more details).
    // Be careful it doesn't check input params, so the
    // interval [low_fr..hig_fr] must be wide enought.
    // Output vector size is zero if error occured.
    vector<unsigned long> locale_mvalsf(double low_fr, double hig_fr,
                                        bool rnlf, bool rnhf);

    // Returns list of index pairs [row,col] of maximal values of Local
    // Energy Function between low_sc and hig_sc scales.
    // rnls and rnhs rounding parameter for low_sc and hig_sc
    // (see scale2idx for more details).
    // Be careful it doesn't check input params, so the
    // interval [low_sc..hig_sc] must be wide enought.
    // Output vector size is zero if error occured.
    vector<unsigned long> locale_mvalss(double low_sc, double hig_sc,
                                        bool rnls, bool rnhs);

    // Returns indeces list of maximal values of Global Energy
    // Function between low_fr and hig_fr frequencies.
    // rnlf and rnhf rounding parameter for low_fr and hig_fr
    // (see freq2idx for more details).
    // Be careful it doesn't check input params, so the
    // interval [low_fr..hig_fr] must be wide enought.
    // Output vector size is zero if error occured.
    vector<unsigned long> globale_mvalsf(double low_fr, double hig_fr,
                                         bool rnlf, bool rnhf);

    // Returns indeces list of maximal values of Global Energy
    // Function between low_sc and hig_sc scales.
    // rnls and rnhs rounding parameter for low_sc and hig_sc
    // (see scale2idx for more details).
    // Be careful it doesn't check input params, so the
    // interval [low_sc..hig_sc] must be wide enought.
    // Output vector size is zero if error occured.
    vector<unsigned long> globale_mvalss(double low_sc, double hig_sc,
                                         bool rnls, bool rnhs);

    // Returns list of index pairs [row,col] of peaks of CWT
    // between low_fr and hig_fr frequencies.
    // rnlf and rnhf rounding parameter for low_fr and hig_fr
    // (see freq2idx for more details).
    // Be careful it doesn't check input params, so the
    // interval [low_fr..hig_fr] must be wide enought.
    // Output vector size is zero if error occured or nothing found.
    vector<unsigned long> cwt_peaksf(double low_fr, double hig_fr,
                                     bool rnlf, bool rnhf);

    // Returns list of index pairs [row,col] of peaks of CWT
    // between low_sc and hig_sc scales.
    // rnls and rnhs rounding parameter for low_sc and hig_sc
    // (see scale2idx for more details).
    // Be careful it doesn't check input params, so the
    // interval [low_sc..hig_sc] must be wide enought.
    // Output vector size is zero if error occured or nothing found.
    vector<unsigned long> cwt_peakss(double low_sc, double hig_sc,
                                     bool rnls, bool rnhs);

    // Returns list of index pairs [row,col] of peaks of Local
    // Energy Function between low_fr and hig_fr frequencies.
    // rnlf and rnhf rounding parameter for low_fr and hig_fr
    // (see freq2idx for more details).
    // Be careful it doesn't check input params, so the
    // interval [low_fr..hig_fr] must be wide enought.
    // Output vector size is zero if error occured or nothing found.
    vector<unsigned long> locale_peaksf(double low_fr, double hig_fr,
                                        bool rnlf, bool rnhf);

    // Returns list of index pairs [row,col] of peaks of Local
    // Energy Function between low_sc and hig_sc scales.
    // rnls and rnhs rounding parameter for low_sc and hig_sc
    // (see scale2idx for more details).
    // Be careful it doesn't check input params, so the
    // interval [low_sc..hig_sc] must be wide enought.
    // Output vector size is zero if error occured or nothing found.
    vector<unsigned long> locale_peakss(double low_sc, double hig_sc,
                                        bool rnls, bool rnhs);

    // Returns indeces list of peaks of Global Energy
    // Function between low_fr and hig_fr frequencies.
    // rnlf and rnhf rounding parameter for low_fr and hig_fr
    // (see freq2idx for more details).
    // Be careful it doesn't check input params, so the
    // interval [low_fr..hig_fr] must be wide enought.
    // Output vector size is zero if error occured or nothing found.
    vector<unsigned long> globale_peaksf(double low_fr, double hig_fr,
                                         bool rnlf, bool rnhf);

    // Returns indeces list of peaks of Global Energy
    // Function between low_sc and hig_sc scales.
    // rnls and rnhs rounding parameter for low_sc and hig_sc
    // (see scale2idx for more details).
    // Be careful it doesn't check input params, so the
    // interval [low_sc..hig_sc] must be wide enought.
    // Output vector size is zero if error occured.
    vector<unsigned long> globale_peakss(double low_sc, double hig_sc,
                                         bool rnls, bool rnhs);

    // Save CWT data to file
    bool save_cwt(char *file);

    // Save CWT Skeleton data to file
    bool save_cwt_skel(char *file);

    // Save Local Energy Function data to file
    bool save_locale(char *file);

    // Save Local Energy Function Skeleton data to file
    bool save_locale_skel(char *file);

    // Save Global Energy Function data to file
    bool save_globale(char *file);

    // Save Global Energy Function maximas data to file
    bool save_globale_maxs(char *file);

    // Save scales data to file
    bool save_scales(char *file);

    // Save frequencies data to file
    bool save_freqs(char *file);

    // Save wavelet translations data to file
    bool save_transls(char *file);

    // Returns true if time series loaded
    bool is_series();

    // Returns true if CWT performed
    bool is_cwt();

    // Returns true if CWT Skeleton calculated
    bool is_cwtsk();

    // Returns true if Global Energy Function calculated
    bool is_globale();

    // Returns true if Global Energy Function maximas calculated
    bool is_globalemxs();

    // Returns true if Local Energy Function calculated
    bool is_locale();

    // Returns true if Local Energy Function Skeleton calculated
    bool is_localesk();

    // Returns true if amplitudes of wavelet coefficients normalized
    // by norm_cwt()
    bool is_norm();

//========== Operators ==========

    void WTClass::operator =(WTClass& a);
};

#endif
