#include <string.h>
#include "math.h"
#include "array.h"
#include "inout.h"
#include "numpy.h"
#include "wavelet.h"
using namespace std;

Wavelet::Wavelet (float fp, float td, float dt, int nt, int flag) {
   this->fp = fp;
   this->td = td;
   this->dt = dt;
   this->nt = nt;
   this->nshift = 0;
   this->iwvl = flag;
   this->wv = new float[nt];
   memset(wv, 0, nt * sizeof(float));
}

void Wavelet::shift() {
   float tb, tw = 0;
   if (iwvl == 1) tw = 1.0 / fp;
   tb = tw - td; 
   if (tb > 0) {
       this->nshift = tb / dt;
       this->td += nshift * dt;
   }
}

void Wavelet::time_derivative (float p) {
   float *wv1 = new float[nt];
   memcpy(wv1, this->wv, nt * sizeof(float));
   for (int i = 1; i < nt - 1; i++)
        this->wv[i] = 0.5 * (wv1[i+1] - wv1[i-1]);
   this->wv[0] = 0;
   this->wv[nt-1] = 0;
}

void Wavelet::mkwvl(float amp) {
   if (iwvl == -1)
      readdata("wavelet.dat", wv, nt);
   else if (iwvl == 0) {                 // pulse
      int itd = td/dt;
      wv[itd] = 1.0;
   }
   else if (iwvl == 1) {            // ricker
      float t = 1.0 / fp;
      float t0 = t - td;
      int ntsh = int(t / dt + 0.5) + 5.1;
      int nts = ntsh * 2 + 1;
      int ntc = (t - t0) / dt + 1;
      int nt1 = ntc - ntsh;
      float b = PI * fp;
      b = b * b;
      float a = 2.0 * b;
      int it1 = max(nt1, 0);
      int it2 = min(nt1 + nts, nt);
      float tt;
      for (int it = it1; it < it2; it ++) {
           tt = it * dt + t0 - t;
           tt = tt * tt;
           wv[it] = amp*(1.0 - a * tt) * exp(-b * tt);
      }
   }
}

void Wavelet::delwvl () {
   delete [] wv;
}
