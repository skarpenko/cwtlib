#include <stdio.h>
#include <math.h>
#include <iostream>
#include "WTClass.h"

#define N 256
double z[N];


int main(int argc, char **argv) {
    float temp;
    FILE *fh;
    long t, i;
    vector<unsigned long> r;
    WTClass *w, *d;
    fstream filo;
   

    for(t=0; t<N; t++)
       z[t] = sin(2*PI*0.01*t);

/*    fh = fopen(argv[1], "r");
    if(!fh) { printf("Fuck!\n"); return 1; }
    for(i=0; i<N; i++) { fscanf(fh, "%f", &temp); z[i] = temp; }
    fclose(fh);
*/
    w = new WTClass(z, N);

    r = w->cwt_mvalsf(0.0, 1.0, false, true);
    if(!r.size()) cout << "R clear" << endl;
    else cout << "R full" << endl;


    cout << "Performing CWT..." << endl;
    w->do_cwt();
    
    cout << "Normalizing amplitudes..." << endl;
    w->norm_cwt();

    cout << "Calculating global energy..." << endl;
    w->global_energy();

    cout << "Calculating local energy..." << endl;
    w->local_energy();
    
/*    cout << "Scales: ";
    for(t=0; t<w->scales.size(); t++)
       cout << w->scales[t] << " ";

    cout << "Freqs: ";
    for(t=0; t<w->freqs.size(); t++)
       cout << w->freqs[t] << " ";
*/
    cout << endl << endl << "Saving CWT..." << endl;
    w->save_cwt("CWT.txt");

    cout << endl << endl << "Saving LE..." << endl;
    w->save_locale("locale.txt");

    cout << endl << endl << "Saving GE..." << endl;
    w->save_globale("globale.txt");

    cout << endl << endl << "Preparation for saving denormalized CWT..." << endl;    
    d = new WTClass();
    *d = *w;
    d->denorm_cwt();

    cout << endl << endl << "Saving CWT..." << endl;
    d->save_cwt("dCWT.txt");

    r = w->cwt_mvalsf(0.0, 1.0, false, true);
    for(int i=0; i<r.size(); i+=2) {
       cout << r[i] << ", " << r[i+1] << " Amp: " << w->cwt.cwt[r[i]][r[i+1]] << endl;
    }

    r = w->globale_mvalsf(0.0, 1.0, false, true);
    for(int i=0; i<r.size(); i++) {
       cout << r[i] << " Amp: " << w->globale[r[i]] << endl;
    }
    w->cwt_skelet(false);
    w->save_cwt_skel("sk_f.txt");
    w->cwt_skelet(true);
    w->save_cwt_skel("sk_t.txt");

cout << endl;
    w->global_energy_maxs();
    r = w->globale_peaksf(0.0, 1.0, false, true);
    for(int i=0; i<r.size(); i++) {
       cout << r[i] << " Amp: " << w->globale[r[i]] << endl;
    }

    w->save_globale_maxs("globmaxs.txt");

    delete w; delete d;
}
