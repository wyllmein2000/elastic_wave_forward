#include "numpy.h"


/***********************************
 * compare two variables
 * ********************************/
int max(int a, int b) {
   if (a > b)
       return a;
   else
       return b;
}

float max(float a, float b) {
   if (a > b)
       return a;
   else
       return b;
}

float max(float *x, int n) {
   float p = x[0];
   for (int i = 1; i < n; i ++)
       if (p < x[i]) p = x[i];
   return p;
}

int min(int a, int b) {
   if (a < b)
       return a;
   else
       return b;
}

float min(float a, float b) {
   if (a < b)
       return a;
   else
       return b;
}

/***********************************
 * interpolation
 * ********************************/

/* 1d interpolation, where both input x and output y have uniform sample */
//void interp (float *y, float *x, float fy, float fx, float dy, float dx, int ny, int nx) {
void interp (float *y, float fy, float dy, int ny, float *x, float fx, float dx, int nx) {
   float a, ys;
   int ix;
   int ky_xb = (fx - fy) / dy;
   int ky_xe = (fx + nx * dx - dx - fy) / dy;

   for (int iy = 0; iy <= ky_xb; iy ++)
        y[iy] = x[0];

   for (int iy = ky_xe + 1; iy < ny; iy ++)
        y[iy] = x[nx-1];

   for (int iy = ky_xb + 1; iy <= ky_xe; iy ++) {
        ys = fy + iy * dy;
        a = (ys - fx) / dx;
        ix = int(a);
        a = a - ix;
        y[iy] = x[ix] * (1 - a) + x[ix+1] * a;
   }
}

/* 2d interpolation, input v[nxi*nzi] and output u[nxo*nzo] */
void interp2d (float *u, float fxo, float fzo, float dxo, float dzo, int nxo, int nzo, float *v, float fxi, float fzi, float dxi, float dzi, int nxi, int nzi) {

   int jx, jz;
   int iu, iv1, iv2, iv3, iv4;
   float x, z, a, b;

   for (int ix = 0; ix < nxo; ix ++) {
        x = fxo + ix * dxo;
        a = (x - fxi) / dxi;
        if (a < 0) a = 0;
        if (a > nxi - 1) a = nxi - 1.0;
        jx = int(a);
        a = a - jx;

        for (int iz = 0; iz < nzo; iz ++) {
             z = fzo + iz * dzo;
             b = (z - fzi) / dzi;
             if (b < 0) b = 0;
             if (b > nzi - 1) b = nzi - 1.0;
             jz = int(b);
             b = b - jz;
             
             iu = ix * nzo + iz;
             iv1 = jx * nzi + jz;
             iv2 = (jx + 1) * nzi + jz;
             iv3 = jx * nzi + jz + 1;
             iv4 = (jx + 1) * nzi + jz + 1;

             u[iu] = (1-a)*(1-b)*v[iv1] + a*(1-b)*v[iv2] + (1-a)*b*v[iv3] + a*b*v[iv4];
        }
   }
}

/***********************************
 * taper
 * ********************************/
void hamfunc (float *x, float p, int nw, int k1, int k2, int k3, int k4, int flag) {
   float q, p1=p, p2=1.0-p;
   if (flag==0) {
      for (int iw=0; iw<nw; iw++) {
          if (iw<k1 || iw>k4)
             x[iw]=0.0;
          else if (iw>=k1 && iw<k2) {
             q=(float)(k2-iw)/(k2-k1);
             x[iw]=p1+p2*cos(q*PI);
          }
          else if (iw>k3 && iw<=k4) {
             q=(float)(iw-k3)/(k4-k3);
             x[iw]=p1+p2*cos(q*PI);
          }
          else
             x[iw]=1.0;
      }
   }
}
