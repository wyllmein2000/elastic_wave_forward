#ifndef NUMPY_H
#define NUMPY_H

#include <math.h>
#define PI 3.1415926
//#define PI (4*atan(1.0))

  int max(int, int);
float max(float, float);
float max(float *, int);
  int min(int, int);
float min(float, float);

//void interp (float *y, float *x, float fy, float fx, float dy, float dx, int ny, int nx);
void interp (float *y, float fy, float dy, int ny, float *x, float fx, float dx, int nx);

void interp2d (float *u, float fxo, float fzo, float dxo, float dzo, int nxo, int nzo, float *v, float fxi, float fzi, float dxi, float dzi, int nxi, int nzi);

void hamfunc (float *x, float p, int nw, int k1, int k2, int k3, int k4, int flag); 

#endif
