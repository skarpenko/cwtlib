//////////////////////////////////////////////////////////////////////////
//                  Continuous Wavelet Transform Class
//
// Author   : Stepan V. Karpenko (carp@mail.ru)
// Date     : 06-09-2004
// Comments :
// History  :
//
//////////////////////////////////////////////////////////////////////////

#include "WTClass.h"


//========== Private ==========

//
// Class initialization routine
//
void WTClass::ClassInit()
{
    // Initialising WT params
    Wavelet = MORLET;
    Amin = 1;
    Astep = 1;
    Amax = 266;
    Bstep = 1;
    Precision = 8;
    WavPoints = 128000;
    WavPart = REAL;
    FixEdges = true;
    
    // Setting flags
    series_flag = false;
    cwt_flag = false;
    globale_flag = false;
    globalem_flag = false;
    locale_flag = false;
    norm_flag = false;
    localesk_flag = false;
    cwtsk_flag = false;

    // free_cwt may crash if cwt structure is not initialized
    cwt.cwt = NULL;
    locale.cwt = NULL;
    cwt_skel.cwt = NULL;
    locale_skel.cwt = NULL;
}

//
// Fixes CWT boundaries
//
bool WTClass::FixCWTBounds()
{
     double avg = 0.0;
     unsigned long i,j;
     vector<double> l;
     cwt_t l_wt;

     // find average value
     for(i=0; i<series.size(); i++) avg += series[i]; avg/=(double)series.size();

     // create series of average values
     for(i=0; i<series.size(); i++) l.push_back(avg);

     // perform continuous wavelet transform
     if(cwto3(&l[0], l.size(), Amin, Astep, Amax, Bstep, Precision,
              &cwtwlets[Wavelet], WavPart, WavPoints, &l_wt))
          return false;

     // subtract decompositions
     for(i=0; i < cwt.rows; i++)
         for(j=0; j < cwt.cols; j++)
               cwt.cwt[i][j] -= l_wt.cwt[i][j];

     // free allocated structures
     free_cwt(&l_wt);

     return true;
}

//
// Calculates skeleton
// w - source, loc - maxima localization
// loc = true - vertical
// loc = false - horizontal
//
bool WTClass::Skeleton(cwt_t *w, bool loc) {
    unsigned long i, j, i_max, j_max;
    cwt_t tmp;

    if(loc) { i_max = w->cols; j_max = w->rows; }
     else   { i_max = w->rows; j_max = w->cols; }

    if(j_max < 3) return false;

    // create temporary cwt structure
    copy_cwt(w, &tmp);

    // calculate skeleton
    for (i = 0; i < i_max; i++)
    {
        for (j = 1; j < j_max - 1; j++)
        {
            if(loc) {
                if( !(tmp.cwt[j][i] > tmp.cwt[j-1][i] &&
                      tmp.cwt[j][i] > tmp.cwt[j+1][i]) )
                        w->cwt[j][i] = 0.0;    
            }
            else {
                if( !(tmp.cwt[i][j] > tmp.cwt[i][j-1] &&
                      tmp.cwt[i][j] > tmp.cwt[i][j+1]) )
                        w->cwt[i][j] = 0.0;
            }
        }
    }

    free_cwt(&tmp); // free allocated data
    return true;
}


//========== Public ==========

//
// Constructor
//
WTClass::WTClass() {
    ClassInit();
}

//
// Constructor with paramteters
// dat - source series, n - data length
//
WTClass::WTClass(double *dat, long n) {
    ClassInit();
    series_add(dat, n);
}

//
// Destructor
//
WTClass::~WTClass() {
    if(series_flag)  { series_flag = false;  series.clear(); }
    if(cwt_flag)     { cwt_flag = false;     free_cwt(&cwt);
                       scales.clear();       freqs.clear();
                       transls.clear(); }
    if(cwtsk_flag)    { cwtsk_flag = false;    free_cwt(&cwt_skel); }
    if(globale_flag)  { globale_flag = false;  globale.clear(); }
    if(globalem_flag) { globalem_flag = false; globale_maxs.clear(); }
    if(locale_flag)   { locale_flag = false;   free_cwt(&locale); }
    if(localesk_flag) { localesk_flag = false; free_cwt(&locale_skel); }
	norm_flag = false;
}

//
// Sets Wavelet Transform parameters
// for more details see cwtlib
//
void WTClass::set_wt_params(unsigned long Wavelet, double Amin, double Astep,
                            double Amax, double Bstep, unsigned long Precision,
                            unsigned long WavPoints, long WavPart, bool FixEdges) {
    // set given parameters
    WTClass::Wavelet = Wavelet;
    WTClass::Amin = Amin;
    WTClass::Astep = Astep;
    WTClass::Amax = Amax;
    WTClass::Bstep = Bstep;
    WTClass::Precision = Precision;
    WTClass::WavPoints = WavPoints;
    WTClass::WavPart = WavPart;
    WTClass::FixEdges = FixEdges;

    // Clear old results
    if(cwt_flag)     { cwt_flag = false;     free_cwt(&cwt);
                       scales.clear();       freqs.clear();
                       transls.clear(); }
    if(cwtsk_flag)    { cwtsk_flag = false;    free_cwt(&cwt_skel); }
    if(globale_flag)  { globale_flag = false;  globale.clear(); }
    if(globalem_flag) { globalem_flag = false; globale_maxs.clear(); }
    if(locale_flag)   { locale_flag = false;   free_cwt(&locale); }
    if(localesk_flag) { localesk_flag = false; free_cwt(&locale_skel); }
	norm_flag = false;
}

//
// Add time series for analysis
// dat - source series, n - data length
//
void WTClass::series_add(double *dat, long n) {
    if(series_flag)  { series_flag = false;  series.clear(); }
    if(cwt_flag)     { cwt_flag = false;     free_cwt(&cwt);
                       scales.clear();       freqs.clear();
                       transls.clear(); }
    if(cwtsk_flag)    { cwtsk_flag = false;    free_cwt(&cwt_skel); }
    if(globale_flag)  { globale_flag = false;  globale.clear(); }
    if(globalem_flag) { globalem_flag = false; globale_maxs.clear(); }
    if(locale_flag)   { locale_flag = false;   free_cwt(&locale); }
    if(localesk_flag) { localesk_flag = false; free_cwt(&locale_skel); }
	norm_flag = false;
	
    for(long i = 0; i < n; i++)
       series.push_back(dat[i]);
    series_flag = true;
}

//
// Performs Continuous Wavelet Transform
//
bool WTClass::do_cwt() {
    double WavFreq;
    double a, b;
    unsigned long dx, dy;
	
    if(cwt_flag)     { cwt_flag = false;     free_cwt(&cwt);
                       scales.clear();       freqs.clear();
                       transls.clear(); }
    if(cwtsk_flag)    { cwtsk_flag = false;    free_cwt(&cwt_skel); }
    if(globale_flag)  { globale_flag = false;  globale.clear(); }
    if(globalem_flag) { globalem_flag = false; globale_maxs.clear(); }
    if(locale_flag)   { locale_flag = false;   free_cwt(&locale); }
    if(localesk_flag) { localesk_flag = false; free_cwt(&locale_skel); }
	norm_flag = false;
	
    // Perform wavelet transform
    if(cwto3(&series[0], series.size(), Amin, Astep, Amax, Bstep, Precision,
             &cwtwlets[Wavelet], WavPart, WavPoints, &cwt))
	return false;
    // Fix CWT boundaries
    if(FixEdges)
     if(!FixCWTBounds()) {
         free_cwt(&cwt);
         return false;
     }
	
    // Select wavelet frequency
    switch(Wavelet) {
	case MORLET: WavFreq = 1.0; break;
	case MEXHAT:
	case GAUSS1:
	case GAUSS2: WavFreq = 0.5; break;
	default:
	    free_cwt(&cwt);
	    return false;
	    break;
    }
	
    // Fill Scales and Freqs vectors
    for (dy = 0, a = cwt.amin; dy < cwt.rows; dy++, a+=cwt.astep) {
	scales.push_back(a);
	freqs.push_back(WavFreq/a);
    }

    // Fill Translations vector
    for (dx = 0, b = 0; dx < cwt.cols; dx++, b+=cwt.bstep) {
	transls.push_back(b);
    }
	
    cwt_flag = true;
    return true;
}

//
// Calculates CWT Skeleton
// loc - localization
//
bool WTClass::cwt_skelet(bool loc) {
    unsigned long dx,dy;
    cwt_t tmp;

    if(!cwt_flag) return false;
    if(cwtsk_flag) { cwtsk_flag = false; free_cwt(&cwt_skel); }

    // make temporary copy
    if(copy_cwt(&cwt, &tmp))
        return false;

    // only absolute values needed
    for (dy = 0; dy < tmp.rows; dy++)
    	for (dx = 0; dx < tmp.cols; dx++)
    	    tmp.cwt[dy][dx] = fabs(tmp.cwt[dy][dx]);

    // Calculate skeleton
    if( !Skeleton(&tmp, loc) ) {
        free_cwt(&tmp);
        return false;
    }

    if(copy_cwt(&cwt, &cwt_skel))
        return false;

    // build skeleton
    for (dy = 0; dy < tmp.rows; dy++)
    	for (dx = 0; dx < tmp.cols; dx++)
    	    if( tmp.cwt[dy][dx] == 0.0 )
    	        cwt_skel.cwt[dy][dx] = 0.0;

    free_cwt(&tmp); // free temporary copy

    cwtsk_loc  = loc;
    cwtsk_flag = true;
    return true;
}

//
// Normalizes amplitudes of wavelet coefficients
// NOTE: This method clears globale, locale, locale_skel and cwt_skel
//
void WTClass::norm_cwt() {
    double a,b;
    unsigned long dx,dy;
	
    if(norm_flag) return;
    
    // Clear previously calculated functions
    if(globale_flag)  { globale_flag = false;  globale.clear(); }
    if(globalem_flag) { globalem_flag = false; globale_maxs.clear(); }
    if(locale_flag)   { locale_flag = false;   free_cwt(&locale); }
    if(cwtsk_flag)    { cwtsk_flag = false;    free_cwt(&cwt_skel); }
    if(localesk_flag) { localesk_flag = false; free_cwt(&locale_skel); }

    // Normalize amplitudes of wavelet coefficients
    for (dy = 0, a = cwt.amin; dy < cwt.rows; dy++, a+=cwt.astep)
    {
        for (dx = 0, b = 0.0; dx < cwt.cols; dx++, b+=cwt.bstep)
        {
            cwt.cwt[dy][dx] *= 2.0 / sqrt(2 * a * PI);
        }
    } // See Foster. This may be incorrect for some wavelets!

    norm_flag = true;
}

//
// Restores amplitudes of wavelet coefficients normalized by norm_cwt()
// NOTE: This method clears globale, locale, locale_skel and cwt_skel
//
void WTClass::denorm_cwt() {
    double a,b;
    unsigned long dx,dy;
	
    if(!norm_flag) return;

    // Clear previously calculated functions
    if(globale_flag)  { globale_flag = false;  globale.clear(); }
    if(globalem_flag) { globalem_flag = false; globale_maxs.clear(); }
    if(locale_flag)   { locale_flag = false;   free_cwt(&locale); }
    if(cwtsk_flag)    { cwtsk_flag = false;    free_cwt(&cwt_skel); }
    if(localesk_flag) { localesk_flag = false; free_cwt(&locale_skel); }

    // Restoring amplitudes of wavelet coefficients
    for (dy = 0, a = cwt.amin; dy < cwt.rows; dy++, a+=cwt.astep)
    {
    	for (dx = 0, b = 0.0; dx < cwt.cols; dx++, b+=cwt.bstep)
    	{
    	    cwt.cwt[dy][dx] *= sqrt(2 * a * PI) / 2.0;
    	}
    } // See Foster. This may be incorrect for some wavelets!

    norm_flag = false;
}

//
// Calculates Global Energy Function
// NOTE: This method clears globale_maxs
//
void WTClass::global_energy() {
    unsigned long dx,dy;
    double temp;
    
    if(globale_flag)  { globale_flag = false; globale.clear(); }
    if(globalem_flag) { globalem_flag = false; globale_maxs.clear(); }

    // Calculate global energy function
    for(dy=0; dy<cwt.rows; dy++) {
        temp = 0.0;
        for(dx=0; dx<cwt.cols; dx++) {
            temp += cwt.cwt[dy][dx] * cwt.cwt[dy][dx];
        }
        globale.push_back(temp/(double)cwt.cols);
    }
    globale_flag = true;
}

//
// Calculates Global Energy Function maximas
//
bool WTClass::global_energy_maxs() {
    unsigned long i;
    
    if(!globale_flag) return false;
    if(globale.size() < 3) return false;
    if(globalem_flag) { globalem_flag = false; globale_maxs.clear(); }

    // Calculate global energy function maximas
    globale_maxs.push_back(0.0);
    for (i = 1; i < globale.size()-1; i++)
    {
        if( globale[i] > globale[i-1] &&
            globale[i] > globale[i+1] ) {
             globale_maxs.push_back( globale[i] );
        } else { globale_maxs.push_back( 0.0 ); }
    }
    globale_maxs.push_back(0.0);

    globalem_flag = true;
    return true;
}

//
// Calculates Local Energy Function
// NOTE: This method clears locale_skel
//
bool WTClass::local_energy() {
    unsigned long dx,dy;
    
    if(locale_flag)   { locale_flag = false;   free_cwt(&locale); }
    if(localesk_flag) { localesk_flag = false; free_cwt(&locale_skel); }

    if(copy_cwt(&cwt, &locale))
        return false;

    // Calculate local energy function
    for (dy = 0; dy < locale.rows; dy++)
    	for (dx = 0; dx < locale.cols; dx++)
    	    locale.cwt[dy][dx] *= locale.cwt[dy][dx];

    locale_flag = true;
    return true;
}

//
// Calculates Local Energy Function Skeleton
// loc - localization
//
bool WTClass::local_energy_skelet(bool loc) {
    unsigned long dx,dy;

    if(!locale_flag) return false;
    if(localesk_flag) { localesk_flag = false; free_cwt(&locale_skel); }

    if(copy_cwt(&locale, &locale_skel))
        return false;

    // Calculate skeleton
    if( !Skeleton(&locale_skel, loc) ) {
        free_cwt(&locale_skel);
        return false;
    }

    localesk_loc  = loc;
    localesk_flag = true;
    return true;
}

//
// Returns row index associated with given frequency.
// rounding controls method's behavior when given
// frequency doesn't exist in freqs[] list.
// =true, frequency will be rounded up
// =false, frequency will be rounded down
//
unsigned long WTClass::freq2idx(double freq, bool rounding) {
    // Check boundaries
    for(unsigned long i=0; i<freqs.size(); i++) {
        if( (!i && freq >= freqs[i]) ||
            (i == freqs.size()-1 && freq <= freqs[i]) )
            return i;

        // Search for frequency
        if( freqs[i] == freq ) return i;
        if( freqs[i] > freq && freqs[i+1] < freq )
            if(rounding)
                return i;
            else
                return (i+1);
    }
}

//
// Returns row index associated with given scale.
// rounding controls method's behavior when given
// scale doesn't exist in scales[] list.
// =true, scale will be rounded up
// =false, scale will be rounded down
//
unsigned long WTClass::scale2idx(double scale, bool rounding) {
    // Check boundaries
    for(unsigned long i=0; i<scales.size(); i++) {
        if( (!i && scale <= scales[i]) ||
            (i == scales.size()-1 && scale >= scales[i]) )
            return i;

        // Search for scale
        if( scales[i] == scale ) return i;
        if( scales[i] < scale && scales[i+1] > scale )
            if(rounding)
                return (i+1);
            else
                return i;
    }
}

//
// Returns row index associated with given wavelet translation.
// rounding controls method's behavior when given
// translation doesn't exist in transls[] list.
// =true, translation will be rounded up
// =false, translation will be rounded down
//
unsigned long WTClass::transl2idx(double transl, bool rounding) {
    // Check boundaries
    for(unsigned long i=0; i<transls.size(); i++) {
        if( (!i && transl <= transls[i]) ||
            (i == transls.size()-1 && transl >= transls[i]) )
            return i;

        // Search for translation
        if( transls[i] == transl ) return i;
        if( transls[i] < transl && transls[i+1] > transl )
            if(rounding)
                return (i+1);
            else
                return i;
    }
}

//
// Returns list of index pairs [row,col] of maximal values of CWT
// between low_fr and hig_fr frequencies.
// rnlf and rnhf rounding parameter for low_fr and hig_fr
// (see freq2idx for more details).
// Be careful it doesn't check input params, so the
// interval [low_fr..hig_fr] must be wide enought.
// Output vector size is zero if error occured.
//
vector<unsigned long> WTClass::cwt_mvalsf(double low_fr, double hig_fr,
                                          bool rnlf, bool rnhf) {
    unsigned long dx, dy, low_y, hig_y;
    vector<unsigned long> ri, ril;

    if(!cwt_flag) return ril;

    low_y = freq2idx(hig_fr, rnhf);
    hig_y = freq2idx(low_fr, rnlf);

    // Search for maximal value inside given interval
    ri.push_back(low_y); ri.push_back(0);
    for (dy = low_y; dy <= hig_y; dy++)
    	for (dx = 0; dx < cwt.cols; dx++)
    	    if( fabs( cwt.cwt[dy][dx] ) > fabs( cwt.cwt[ ri[0] ][ ri[1] ] ) )
            { ri[0] = dy; ri[1] = dx; }

    // Search for other values that is equal to already discovered
    for (dy = low_y; dy <= hig_y; dy++)
    	for (dx = 0; dx < cwt.cols; dx++)
    	    if( fabs( cwt.cwt[dy][dx] ) == fabs( cwt.cwt[ ri[0] ][ ri[1] ] ) )
            { ril.push_back(dy); ril.push_back(dx); }

    return ril;
}

//
// Returns list of index pairs [row,col] of maximal values of CWT
// between low_sc and hig_sc scales.
// rnls and rnhs rounding parameter for low_sc and hig_sc
// (see scale2idx for more details).
// Be careful it doesn't check input params, so the
// interval [low_sc..hig_sc] must be wide enought.
// Output vector size is zero if error occured.
//
vector<unsigned long> WTClass::cwt_mvalss(double low_sc, double hig_sc,
                                          bool rnls, bool rnhs) {
    unsigned long dx, dy, low_y, hig_y;
    vector<unsigned long> ri, ril;

    if(!cwt_flag) return ril;

    low_y = scale2idx(low_sc, rnls);
    hig_y = scale2idx(hig_sc, rnhs);

    // Search for maximal value inside given interval
    ri.push_back(low_y); ri.push_back(0);
    for (dy = low_y; dy <= hig_y; dy++)
    	for (dx = 0; dx < cwt.cols; dx++)
    	    if( fabs( cwt.cwt[dy][dx] ) > fabs( cwt.cwt[ ri[0] ][ ri[1] ] ) )
            { ri[0] = dy; ri[1] = dx; }

    // Search for other values that is equal to already discovered
    for (dy = low_y; dy <= hig_y; dy++)
    	for (dx = 0; dx < cwt.cols; dx++)
    	    if( fabs( cwt.cwt[dy][dx] ) == fabs( cwt.cwt[ ri[0] ][ ri[1] ] ) )
            { ril.push_back(dy); ril.push_back(dx); }

    return ril;
}

//
// Returns list of index pairs [row,col] of maximal values of Local
// Energy Function between low_fr and hig_fr frequencies.
// rnlf and rnhf rounding parameter for low_fr and hig_fr
// (see freq2idx for more details).
// Be careful it doesn't check input params, so the
// interval [low_fr..hig_fr] must be wide enought.
// Output vector size is zero if error occured.
//
vector<unsigned long> WTClass::locale_mvalsf(double low_fr, double hig_fr,
                                             bool rnlf, bool rnhf) {
    unsigned long dx, dy, low_y, hig_y;
    vector<unsigned long> ri, ril;

    if(!locale_flag) return ril;

    low_y = freq2idx(hig_fr, rnhf);
    hig_y = freq2idx(low_fr, rnlf);

    // Search for maximal value inside given interval
    ri.push_back(low_y); ri.push_back(0);
    for (dy = low_y; dy <= hig_y; dy++)
    	for (dx = 0; dx < locale.cols; dx++)
    	    if( locale.cwt[dy][dx] < locale.cwt[ ri[0] ][ ri[1] ] )
    	    { ri[0] = dy; ri[1] = dx; }

    // Search for other values that is equal to already discovered
    for (dy = low_y; dy <= hig_y; dy++)
    	for (dx = 0; dx < locale.cols; dx++)
    	    if( locale.cwt[dy][dx] == locale.cwt[ ri[0] ][ ri[1] ] )
            { ril.push_back(dy); ril.push_back(dx); }

    return ril;
}

//
// Returns list of index pairs [row,col] of maximal values of Local
// Energy Function between low_sc and hig_sc scales.
// rnls and rnhs rounding parameter for low_sc and hig_sc
// (see scale2idx for more details).
// Be careful it doesn't check input params, so the
// interval [low_sc..hig_sc] must be wide enought.
// Output vector size is zero if error occured.
//
vector<unsigned long> WTClass::locale_mvalss(double low_sc, double hig_sc,
                                             bool rnls, bool rnhs) {
    unsigned long dx, dy, low_y, hig_y;
    vector<unsigned long> ri, ril;

    if(!locale_flag) return ril;

    low_y = scale2idx(low_sc, rnls);
    hig_y = scale2idx(hig_sc, rnhs);

    // Search for maximal value inside given interval
    ri.push_back(low_y); ri.push_back(0);
    for (dy = low_y; dy <= hig_y; dy++)
    	for (dx = 0; dx < locale.cols; dx++)
    	    if( locale.cwt[dy][dx] < locale.cwt[ ri[0] ][ ri[1] ] )
    	    { ri[0] = dy; ri[1] = dx; }

    // Search for other values that is equal to already discovered
    for (dy = low_y; dy <= hig_y; dy++)
    	for (dx = 0; dx < locale.cols; dx++)
    	    if( locale.cwt[dy][dx] == locale.cwt[ ri[0] ][ ri[1] ] )
            { ril.push_back(dy); ril.push_back(dx); }

    return ril;
}

//
// Returns indeces list of maximal values of Global Energy
// Function between low_fr and hig_fr frequencies.
// rnlf and rnhf rounding parameter for low_fr and hig_fr
// (see freq2idx for more details).
// Be careful it doesn't check input params, so the
// interval [low_fr..hig_fr] must be wide enought.
// Output vector size is zero if error occured.
//
vector<unsigned long> WTClass::globale_mvalsf(double low_fr, double hig_fr,
                                              bool rnlf, bool rnhf) {
    unsigned long i, ri, low_i, hig_i;
    vector<unsigned long> ril;

    if(!globale_flag) return ril;

    low_i = freq2idx(hig_fr, rnhf);
    hig_i = freq2idx(low_fr, rnlf);

    if( !(low_i - hig_i) ) { ril.push_back(low_i); return ril; }

    // Search for maximal value inside given interval
    ri = low_i;
    for(i=low_i+1; i<=hig_i; i++)
        if( globale[i] > globale[ri] ) ri = i;

    // Search for other values that is equal to globale[ri]
    for(i=low_i; i<=hig_i; i++)
       if( globale[i] == globale[ri] ) ril.push_back(i);

    return ril;
}

//
// Returns indeces list of maximal values of Global Energy
// Function between low_sc and hig_sc scales.
// rnls and rnhs rounding parameter for low_sc and hig_sc
// (see scale2idx for more details).
// Be careful it doesn't check input params, so the
// interval [low_sc..hig_sc] must be wide enought.
// Output vector size is zero if error occured.
//
vector<unsigned long> WTClass::globale_mvalss(double low_sc, double hig_sc,
                                              bool rnls, bool rnhs) {
    unsigned long i, ri, low_i, hig_i;
    vector<unsigned long> ril;

    if(!globale_flag) return ril;

    low_i = scale2idx(low_sc, rnls);
    hig_i = scale2idx(hig_sc, rnhs);

    if( !(low_i - hig_i) ) { ril.push_back(low_i); return ril; }

    // Search for maximal value inside given interval
    ri = low_i;
    for(i=low_i+1; i<=hig_i; i++)
        if( globale[i] > globale[ri] ) ri = i;

    // Search for other values that is equal to globale[ri]
    for(i=low_i; i<=hig_i; i++)
       if( globale[i] == globale[ri] ) ril.push_back(i);

    return ril;
}

//
// Returns list of index pairs [row,col] of peaks of CWT
// between low_fr and hig_fr frequencies.
// rnlf and rnhf rounding parameter for low_fr and hig_fr
// (see freq2idx for more details).
// Be careful it doesn't check input params, so the
// interval [low_fr..hig_fr] must be wide enought.
// Output vector size is zero if error occured or nothing found.
//
vector<unsigned long> WTClass::cwt_peaksf(double low_fr, double hig_fr,
                                          bool rnlf, bool rnhf) {
    unsigned long dx, dy, low_y, hig_y;
    vector<unsigned long> ril;

    if(!cwtsk_flag) return ril;

    low_y = freq2idx(hig_fr, rnhf);
    hig_y = freq2idx(low_fr, rnlf);

    // Search for peaks inside given interval
    for (dy = low_y; dy <= hig_y; dy++)
        for (dx = 0; dx < cwt_skel.cols; dx++)
            if( cwt_skel.cwt[dy][dx] != 0.0 )
            { ril.push_back(dy); ril.push_back(dx); }

    return ril;
}

//
// Returns list of index pairs [row,col] of peaks of CWT
// between low_sc and hig_sc scales.
// rnls and rnhs rounding parameter for low_sc and hig_sc
// (see scale2idx for more details).
// Be careful it doesn't check input params, so the
// interval [low_sc..hig_sc] must be wide enought.
// Output vector size is zero if error occured or nothing found.
//
vector<unsigned long> WTClass::cwt_peakss(double low_sc, double hig_sc,
                                          bool rnls, bool rnhs) {
    unsigned long dx, dy, low_y, hig_y;
    vector<unsigned long> ril;

    if(!cwtsk_flag) return ril;

    low_y = scale2idx(low_sc, rnls);
    hig_y = scale2idx(hig_sc, rnhs);

    // Search for peaks inside given interval
    for (dy = low_y; dy <= hig_y; dy++)
        for (dx = 0; dx < cwt_skel.cols; dx++)
            if( cwt_skel.cwt[dy][dx] != 0.0 )
            { ril.push_back(dy); ril.push_back(dx); }

    return ril;
}

//
// Returns list of index pairs [row,col] of peaks of Local
// Energy Function between low_fr and hig_fr frequencies.
// rnlf and rnhf rounding parameter for low_fr and hig_fr
// (see freq2idx for more details).
// Be careful it doesn't check input params, so the
// interval [low_fr..hig_fr] must be wide enought.
// Output vector size is zero if error occured or nothing found.
//
vector<unsigned long> WTClass::locale_peaksf(double low_fr, double hig_fr,
                                             bool rnlf, bool rnhf) {
    unsigned long dx, dy, low_y, hig_y;
    vector<unsigned long> ril;

    if(!localesk_flag) return ril;

    low_y = freq2idx(hig_fr, rnhf);
    hig_y = freq2idx(low_fr, rnlf);

    // Search for peaks inside given interval
    for (dy = low_y; dy <= hig_y; dy++)
        for (dx = 0; dx < locale_skel.cols; dx++)
            if( locale_skel.cwt[dy][dx] != 0.0 )
            { ril.push_back(dy); ril.push_back(dx); }

    return ril;
}

//
// Returns list of index pairs [row,col] of peaks of Local
// Energy Function between low_sc and hig_sc scales.
// rnls and rnhs rounding parameter for low_sc and hig_sc
// (see scale2idx for more details).
// Be careful it doesn't check input params, so the
// interval [low_sc..hig_sc] must be wide enought.
// Output vector size is zero if error occured or nothing found.
//
vector<unsigned long> WTClass::locale_peakss(double low_sc, double hig_sc,
                                             bool rnls, bool rnhs) {
    unsigned long dx, dy, low_y, hig_y;
    vector<unsigned long> ril;

    if(!localesk_flag) return ril;

    low_y = scale2idx(low_sc, rnls);
    hig_y = scale2idx(hig_sc, rnhs);

    // Search for peaks inside given interval
    for (dy = low_y; dy <= hig_y; dy++)
        for (dx = 0; dx < locale_skel.cols; dx++)
            if( locale_skel.cwt[dy][dx] != 0.0 )
            { ril.push_back(dy); ril.push_back(dx); }

    return ril;
}

//
// Returns indeces list of peaks of Global Energy
// Function between low_fr and hig_fr frequencies.
// rnlf and rnhf rounding parameter for low_fr and hig_fr
// (see freq2idx for more details).
// Be careful it doesn't check input params, so the
// interval [low_fr..hig_fr] must be wide enought.
// Output vector size is zero if error occured or nothing found.
//
vector<unsigned long> WTClass::globale_peaksf(double low_fr, double hig_fr,
                                              bool rnlf, bool rnhf) {
    unsigned long i, low_i, hig_i;
    vector<unsigned long> ril;

    if(!globalem_flag) return ril;

    low_i = freq2idx(hig_fr, rnhf);
    hig_i = freq2idx(low_fr, rnlf);

    if( !(low_i - hig_i) ) { ril.push_back(low_i); return ril; }

    // Search for peaks inside given interval
    for(i=low_i; i<=hig_i; i++)
        if( globale_maxs[i] != 0.0 ) ril.push_back(i);

    return ril;
}

//
// Returns indeces list of peaks of Global Energy
// Function between low_sc and hig_sc scales.
// rnls and rnhs rounding parameter for low_sc and hig_sc
// (see scale2idx for more details).
// Be careful it doesn't check input params, so the
// interval [low_sc..hig_sc] must be wide enought.
// Output vector size is zero if error occured.
//
vector<unsigned long> WTClass::globale_peakss(double low_sc, double hig_sc,
                                              bool rnls, bool rnhs) {
    unsigned long i, low_i, hig_i;
    vector<unsigned long> ril;

    if(!globalem_flag) return ril;

    low_i = scale2idx(low_sc, rnls);
    hig_i = scale2idx(hig_sc, rnhs);

    if( !(low_i - hig_i) ) { ril.push_back(low_i); return ril; }

    // Search for peaks inside given interval
    for(i=low_i; i<=hig_i; i++)
        if( globale_maxs[i] != 0.0 ) ril.push_back(i);

    return ril;
}

//
// Save CWT data to file
//
bool WTClass::save_cwt(char *file) {
    ofstream fstrm;

    if(!cwt_flag) return false;

    // Open file
    fstrm.open(file);
    if(!fstrm) return false;

    fstrm.precision(16);

    // Store data
    for(unsigned long dy=0; dy<cwt.rows; dy++) {
        for(unsigned long dx=0; dx<cwt.cols; dx++) {
            fstrm << " " << cwt.cwt[dy][dx];
        }
        fstrm << endl;
    }
    fstrm.close();

    return true;
}

//
// Save CWT Skeleton data to file
//
bool WTClass::save_cwt_skel(char *file) {
    ofstream fstrm;

    if(!cwtsk_flag) return false;

    // Open file
    fstrm.open(file);
    if(!fstrm) return false;

    fstrm.precision(16);

    // Store data
    for(unsigned long dy=0; dy<cwt_skel.rows; dy++) {
        for(unsigned long dx=0; dx<cwt_skel.cols; dx++) {
            fstrm << " " << cwt_skel.cwt[dy][dx];
        }
        fstrm << endl;
    }
    fstrm.close();

    return true;
}

//
// Save Local Energy Function data to file
//
bool WTClass::save_locale(char *file) {
    ofstream fstrm;

    if(!locale_flag) return false;

    // Open file
    fstrm.open(file);
    if(!fstrm) return false;

    fstrm.precision(16);

    // Store data
    for(unsigned long dy=0; dy<locale.rows; dy++) {
        for(unsigned long dx=0; dx<locale.cols; dx++) {
            fstrm << " " << locale.cwt[dy][dx];
        }
        fstrm << endl;
    }
    fstrm.close();

    return true;
}

//
// Save Local Energy Function Skeleton data to file
//
bool WTClass::save_locale_skel(char *file) {
    ofstream fstrm;

    if(!localesk_flag) return false;

    // Open file
    fstrm.open(file);
    if(!fstrm) return false;

    fstrm.precision(16);

    // Store data
    for(unsigned long dy=0; dy<locale_skel.rows; dy++) {
        for(unsigned long dx=0; dx<locale_skel.cols; dx++) {
            fstrm << " " << locale_skel.cwt[dy][dx];
        }
        fstrm << endl;
    }
    fstrm.close();

    return true;
}

//
// Save Global Energy Function data to file
//
bool WTClass::save_globale(char *file) {
    ofstream fstrm;

    if(!globale_flag) return false;

    // Open file
    fstrm.open(file);
    if(!fstrm) return false;

    fstrm.precision(16);

    // Store data
    for(unsigned long i=0; i<globale.size(); i++) {
        fstrm << globale[i] << endl;
    }
    fstrm.close();

    return true;
}

//
// Save Global Energy Function maximas data to file
//
bool WTClass::save_globale_maxs(char *file) {
    ofstream fstrm;

    if(!globalem_flag) return false;

    // Open file
    fstrm.open(file);
    if(!fstrm) return false;

    fstrm.precision(16);

    // Store data
    for(unsigned long i=0; i<globale_maxs.size(); i++) {
        fstrm << globale_maxs[i] << endl;
    }
    fstrm.close();

    return true;
}

//
// Save scales data to file
//
bool WTClass::save_scales(char *file) {
    ofstream fstrm;

    if(!cwt_flag) return false;

    // Open file
    fstrm.open(file);
    if(!fstrm) return false;

    fstrm.precision(16);

    // Store data
    for(unsigned long i=0; i<scales.size(); i++) {
        fstrm << scales[i] << endl;
    }
    fstrm.close();

    return true;
}

//
// Save frequencies data to file
//
bool WTClass::save_freqs(char *file) {
    ofstream fstrm;

    if(!cwt_flag) return false;

    // Open file
    fstrm.open(file);
    if(!fstrm) return false;

    fstrm.precision(16);

    // Store data
    for(unsigned long i=0; i<freqs.size(); i++) {
        fstrm << freqs[i] << endl;
    }
    fstrm.close();

    return true;
}

//
// Save wavelet translations data to file
//
bool WTClass::save_transls(char *file) {
    ofstream fstrm;

    if(!cwt_flag) return false;

    // Open file
    fstrm.open(file);
    if(!fstrm) return false;

    fstrm.precision(16);

    // Store data
    for(unsigned long i=0; i<transls.size(); i++) {
        fstrm << transls[i] << endl;
    }
    fstrm.close();

    return true;
}

//
// Returns true if time series loaded
//
bool WTClass::is_series() {
    return series_flag;
}

//
// Returns true if CWT performed
//
bool WTClass::is_cwt() {
    return cwt_flag;
}

//
// Returns true if CWT Skeleton calculated
//
bool WTClass::is_cwtsk() {
    return cwtsk_flag;
}

//
// Returns true if Global Energy Function calculated
//
bool WTClass::is_globale() {
    return globale_flag;
}

//
// Returns true if Global Energy Function maximas calculated
//
bool WTClass::is_globalemxs() {
    return globalem_flag;
}

//
// Returns true if Local Energy Function calculated
//
bool WTClass::is_locale() {
    return locale_flag;
}

//
// Returns true if Local Energy Function Skeleton calculated
//
bool WTClass::is_localesk() {
    return localesk_flag;
}

//
// Returns true if amplitudes of wavelet coefficients normalized
// by norm_cwt()
//
bool WTClass::is_norm() {
    return norm_flag;
}

//========== Operators ==========

void WTClass::operator =(WTClass& a) {
    // check against a=a
    if(this == &a) return;

    // clear local data
    if(series_flag)  { series_flag = false;  series.clear(); }
    if(cwt_flag)     { cwt_flag = false;     free_cwt(&cwt);
                       scales.clear();       freqs.clear();
                       transls.clear(); }
    if(cwtsk_flag)    { cwtsk_flag = false;    free_cwt(&cwt_skel); }
    if(globale_flag)  { globale_flag = false;  globale.clear(); }
    if(globalem_flag) { globalem_flag = false; globale_maxs.clear(); }
    if(locale_flag)   { locale_flag = false;   free_cwt(&locale); }
    if(localesk_flag) { localesk_flag = false; free_cwt(&locale_skel); }
    norm_flag = false;

    // copy wavelet transform parameters
    Wavelet   = a.Wavelet;
    Amin      = a.Amin;
    Astep     = a.Astep;
    Amax      = a.Amax;
    Bstep     = a.Bstep;
    Precision = a.Precision;
    WavPoints = a.WavPoints;
    WavPart   = a.WavPart;
    FixEdges  = a.FixEdges;

    // copy flags
    series_flag   = a.series_flag;
    cwt_flag      = a.cwt_flag;
    globale_flag  = a.globale_flag;
    globalem_flag = a.globalem_flag;
    locale_flag   = a.locale_flag;
    norm_flag     = a.norm_flag;
    localesk_flag = a.localesk_flag;
    cwtsk_flag    = a.cwtsk_flag;
    cwtsk_loc     = a.cwtsk_loc;
    localesk_loc  = a.localesk_loc;

    // copy data
    series       = a.series;
    scales       = a.scales;
    transls      = a.transls;
    freqs        = a.freqs;
    globale      = a.globale;
    globale_maxs = a.globale_maxs;

    if(cwt_flag)      copy_cwt(&a.cwt, &cwt);
    if(locale_flag)   copy_cwt(&a.locale, &locale);
    if(cwtsk_flag)    copy_cwt(&a.cwt_skel, &cwt_skel);
    if(localesk_flag) copy_cwt(&a.locale_skel, &locale_skel);
}
