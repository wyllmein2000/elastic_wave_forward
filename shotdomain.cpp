#include "iostream"
#include "stdlib.h"
#include "numpy.h"
#include "para.h"
#define DIM 2
using namespace std;


void GetShotDomain (Shot *srcfd, Shot src, float fp, float fmax, float vmin, float vmax, int nxor, int ntor) {
   int dim = DIM;
   float dx = src.dx;
   float dz = src.dz;
   float dt = src.dt;

   float num_per_wavelength = 8.0;     // dispersion
   float stablility = 0.5;             // stablility
   
   if (dim == 2) {
       if (ntor == 2 && nxor == 2) {
           num_per_wavelength = 8.2;
           stablility = 0.707;
       }
       else if (ntor == 2 && nxor == 4) {
           num_per_wavelength = 8.0;
           stablility = 0.606;
       }
       else if (ntor == 2 && nxor == 6) {
           //num_per_wavelength = 4.0;
           num_per_wavelength = 5.0;
       }
       else if (ntor == 2 && nxor == 10) {
           num_per_wavelength = 5.0;
       }
   }
   else if (dim == 3) {
       if (ntor == 2 && nxor == 2) {
       }
       else if (ntor == 2 && nxor == 4) {
           stablility = 0.495;
       }
       else if (ntor == 2 && nxor == 6) {
       }
   }

   float p = 1.2;
   float wavelength = vmin / fmax;
   //float wavelength = vmin / fp;
   float dxc = wavelength / num_per_wavelength;

   int gx = dx / (p * dxc) + 1;
   int gz = dz / (p * dxc) + 1;
   dx = dx / gx;
   dz = dz / gz;
   //cout << " dx =" << dx << " dz =" << dz << endl; exit(0);
   
   float dtc = stablility * dx / vmax;
   if (dt > dtc) dt = dtc;

   float px = src.dx / dx;
   float pz = src.dz / dz;
   float pt = src.dt / dt;

   int nx = (src.nx - 1) * px + 1;
   int nz = (src.nz - 1) * pz + 1;
   int nt = (src.nt - 1) * pt + 1;

   srcfd->nx = nx;
   srcfd->nz = nz;
   srcfd->nt = nt;
   srcfd->dx = dx;
   srcfd->dz = dz;
   srcfd->dt = dt;

   srcfd->xs = src.xs;
   srcfd->zs = src.zs;
   srcfd->xg = src.xg;
   srcfd->zg = src.zg;

   srcfd->ixs = src.ixs * px;
   srcfd->izs = src.izs * pz;
   srcfd->ixg = src.ixg * px; 
   srcfd->izg = src.izg * pz; 
  
   srcfd->zr = new float[nx]; 
   srcfd->izr = new int[nx]; 
 
   //interp (srcfd->zr, src.zr, 0.0, 0.0, dx, src.dx, nx, src.nx);
   interp (srcfd->zr, 0.0, dx, nx, src.zr, 0.0, src.dx, src.nx);
   for (int ix = 0; ix < nx; ix ++)
        srcfd->izr[ix] = srcfd->zr[ix] / dz;

}
