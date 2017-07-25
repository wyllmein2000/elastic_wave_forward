#ifndef WAVELET_H
#define WAVELET_H
#include "para.h"

class Wavelet {
   public:
     int iwvl;         // wavelet type
     int nt;           // length
     int nshift;       // shift backwards
     float fp;         // peak frequency
     float td;         // peak delay
     float dt;         // 
     float *wv; 
   public:
     Wavelet (float, float, float, int, int);
     void mkwvl(float);
     void shift();
     void time_derivative(float);
     void delwvl();
};

#endif
