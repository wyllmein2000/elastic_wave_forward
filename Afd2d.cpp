#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include "omp.h"
#include "math.h"
#include "numpy.h"
#include "inout.h"
#include "Afd2d.h"
using namespace std;

/******************************************
 *     input size: nx = nxi + 2 * nxpml
 * wavefield size: nxp = nx + 2 * nor;
 *  pml w.f. size: nxpml
 * ***************************************/
fd2d::fd2d (int nx, int nz, int nt, float dx, float dz, float dt, int nxb, int nzb, int order) {
   this->nor = order / 2;              // computation grid
   this->nx = nx;
   this->nz = nz;
   this->nt = nt;
   this->dx = dx;
   this->dz = dz;
   this->dt = dt;
   this->nxb = nxb;
   this->nzb = nzb;
   this->nxp = nx + 2 * nor;
   this->nzp = nz + 2 * nor;

   this->snap_fx = nor * dx;           // snapshot
   this->snap_fy = nor * dy;
   this->snap_fz = nor * dz;
   this->snap_dx = dx;
   this->snap_dy = dy;
   this->snap_dz = dz;
   this->snap_dt = dt;
   this->snap_nx = nx;
   this->snap_ny = ny;
   this->snap_nz = nz;
   this->snap_nt = nt;

   this->mwf = nxp * nzp;             // total wavefield size
   this->mxp = 2 * nxb * nz;          // pml wavefield size
   this->mzp = 2 * nx * nzb;
   
   this->p1 = new float[mwf];
   this->p2 = new float[mwf];
   this->p3 = new float[mwf];

   //this->DPxx = new float[mxp];
   //this->DPzz = new float[mzp];
   this->DPxx = new float[nx*nz];
   this->DPzz = new float[nx*nz];

   this->dpdx = new float[mxp];
   this->dedx = new float[mxp];
   this->dpdz = new float[mzp];
   this->dedz = new float[mzp];

   this->Expml = new float[mxp];
   this->Hxpml = new float[mxp];
   this->Ezpml = new float[mzp];
   this->Hzpml = new float[mzp];

   this->axpml = new float[2*nxb];
   this->azpml = new float[2*nzb];
   
}


void fd2d::initvar() {
   memset(p1, 0, mwf * sizeof(float));
   memset(p2, 0, mwf * sizeof(float));
   memset(p3, 0, mwf * sizeof(float));

   //memset(DPxx, 0, mxp * sizeof(float));
   //memset(DPzz, 0, mzp * sizeof(float));
   memset(DPxx, 0, nx * nz * sizeof(float));
   memset(DPzz, 0, nx * nz * sizeof(float));

   memset(dpdx, 0, mxp * sizeof(float));
   memset(dedx, 0, mxp * sizeof(float));
   memset(dpdz, 0, mzp * sizeof(float));
   memset(dedz, 0, mzp * sizeof(float));

   memset(Expml, 0, mxp * sizeof(float));
   memset(Hxpml, 0, mxp * sizeof(float));
   memset(Ezpml, 0, mzp * sizeof(float));
   memset(Hzpml, 0, mzp * sizeof(float));
}

void fd2d::initabc (int iz_free_surface, float fp, float *v) {
  this->ifs = 0;
  this->izf = iz_free_surface;
  this->iztop = 0;
  this->nbtop = nzb;

  if (this->izf >= 0) {              // free surface
      this->ifs == 1;
      this->iztop = izf + 1;
      this->nbtop = 0;
  }

  getpmlcoef (v);     // pml
}

void fd2d::snapset (float dxo, float dyo, float dzo, float dto, int nxo, int nyo, int nzo, int nto) {
  this->snap_fx = (nxb + nor) * dx;
  this->snap_fy = (nyb + nor) * dy;
  this->snap_fz = (nzb + nor) * dz;
  this->snap_dx = dxo;
  this->snap_dy = dyo;
  this->snap_dz = dzo;
  this->snap_dt = dto;
  this->snap_nx = nxo;
  this->snap_ny = nyo;
  this->snap_nz = nzo;
  this->snap_nt = nto;
}

void fd2d::forward(float *u, float *v, int ixs, int izs, int *izr, int iw) {
  int jx, jz;
  int ktm = snap_dt / dt;
  float dtdx = dt / dx;
  float dtdz = dt / dz;
  float dtxz = dt / (dx * dz);
  float *v2, *ptmp, *snapshot;

  getfdcoef (); 
  initvar ();

  v2 = new float[nx * nz];
  for (int i = 0; i < nx * nz; i ++)
       v2[i] = v[i] * v[i];

  /*
  #pragma omp parallel
  {
     int nthread = omp_get_num_threads();
     int ithread = omp_get_thread_num();
     #pragma omp critical
     cout << " OMP ithread = " << ithread << "/" << nthread << endl;
  }*/

  ofstream fsnap;
  if (iw == 1) {
      snapshot = new float[snap_nx * snap_nz];
      fsnap.open("temp/wavefield.bin", ios::out|ios::binary);
      cout << " model size: nx=" << nx << " nz=" << nz << endl;
      cout << "             dx=" << dx << " dz=" << dz << " dt=" << dt << endl;
      cout << " save u(t,z,x): nx=" << snap_nx << " nz=" << snap_nz << " nt=" << nt/ktm << endl;
      cout << " FD x-order=" << nor << endl;
      cout << " FD COEF: c0=" << coef[0] << " c1=" << coef[1] << " c2=" << coef[2] << endl;
  }

  for (int it = 0; it < nt; it ++) {

       if (iw == 1 && it % 500 == 0)
           cout << " *** it = " << it << "/" << nt << endl;

       update (v2, iztop);
       update_pml (it, v2);
       if (ifs == 1) update_fs ();

       // load source
       for (int ix = 0; ix < nx; ix ++)
            p1[IU(izs,ix)] += u[ix*nt+it] * pow(v[IV(izs,ix)]*dt,2)/(dx*dz);

       // output trace
       for (int ix = 0; ix < nx; ix ++) {
            jz = izr[ix];            
            u[ix*nt+it] = p1[IU(jz,ix)];
       }

       // output snapshot
       if (iw == 1 && it % ktm == 0) {
           interp2d (snapshot, snap_fx, snap_fz, snap_dx, snap_dz,
                     snap_nx, snap_nz, p1, 0.0, 0.0, dx, dz, nxp, nzp);
           fsnap.write((char *)snapshot, snap_nx * snap_nz * sizeof(float));
       }

       ptmp = p3;
       p3 = p2;
       p2 = p1;
       p1 = ptmp; 
  }

  if (iw == 1) {
      delete [] snapshot;
      fsnap.close();
  }

  delete [] v2;
  freevar ();
}


/****************************************************/
void fd2d::update(float *v2, int iztop) {

  int ix, iz, jx, jz;
   float t2 = dt * dt;
   float xs = 1.0 / (dx * dx);
   float zs = 1.0 / (dz * dz);

   #pragma omp parallel for private (ix, iz)
   for (ix = 0; ix < nx; ix ++)
   for (iz = iztop; iz < nz; iz++) {
        int ip = IU(iz, ix);
        int iv = IV(iz, ix);
        float px = coef[0] * p2[ip];
        float pz = coef[0] * p2[ip];
        for (int ic = 1; ic <= nor; ic++) {
            px += coef[ic]*(p2[IU(iz,ix+ic)]+p2[IU(iz,ix-ic)]);
            pz += coef[ic]*(p2[IU(iz+ic,ix)]+p2[IU(iz-ic,ix)]);
        }
        px *= xs;
        pz *= zs;

        p1[ip] = 2.0 * p2[ip] - p3[ip] + t2 * v2[iv] * (px+pz);

        DPxx[iv] = px;
        DPzz[iv] = pz;

        //if (ix < nxb) DPxx[iv] = px;
        //if (ix >= nx - nxb)  
        //if (iz < nzb) DPzz[iv] = pz;
   }
}

void fd2d::update_fs () {

   int ix, iz, jx, jz;

   #pragma omp parallel for private (ix, iz, jz)
   for (ix = 0; ix < nx; ix++) {
        p1[IU(izf,ix)] = 0;
        for (iz = 0; iz < izf; iz++) {
             jz = 2 * izf - iz;
             p1[IU(iz,ix)] = -p1[IU(jz,ix)];
   }
   }
}


void fd2d::update_pml (int istep, float *v2) {

   int nxx = 2 * nxb;
   int nzz = 2 * nzb;
   int ix, iz, jx, jz, kx, kz;
   float t2 = dt * dt;

   memset(dedx, 0, mxp*sizeof(float));
   memset(dedz, 0, mzp*sizeof(float));
   memset(dpdx, 0, mxp*sizeof(float));
   memset(dpdz, 0, mzp*sizeof(float));

   // x direction
   #pragma omp parallel for private (ix, jx, iz, jz, kx)
   for (iz=0; iz<nz ; iz++) {
   for (ix=1; ix<nxb; ix++) {
       jx=ix+nxb;
       kx=ix+nx-nxb;
       dpdx[IV(iz,ix)]=0.5*(p2[IU(iz,kx+1)]-p2[IU(iz,kx-1)])/dx;
       dpdx[IV(iz,jx)]=0.5*(p2[IU(iz,ix+1)]-p2[IU(iz,ix-1)])/dx;
   }
   }

   #pragma omp parallel for private (ix, iz)
   for (iz=0; iz<nz;  iz++)
   for (ix=1; ix<nxx-1; ix++)
       Expml[IV(iz,ix)]=axpml[ix]*Expml[IV(iz,ix)]+(1.0-axpml[ix])*dpdx[IV(iz,ix)];

   #pragma omp parallel for private (ix, iz)
   for (iz=0; iz<nz;  iz++)
   for (ix=1; ix<nxx-1; ix++)
       dedx[IV(iz,ix)]=0.5*(Expml[IV(iz,ix+1)]-Expml[IV(iz,ix-1)])/dx;

   #pragma omp parallel for private (ix, iz, jx, kx)
   for (iz = 0; iz < nz ; iz++) {
   for (ix = 1; ix < nxb; ix++) {
       jx = ix + nxb;
       kx = ix + nx - nxb;
       int iv = IV(iz, ix);
       int jv = IV(iz, jx);
       Hxpml[iv] = axpml[ix] * Hxpml[iv] + (1.0 - axpml[ix]) * (DPxx[IV(iz,kx)] - dedx[iv]);
       Hxpml[jv] = axpml[jx] * Hxpml[jv] + (1.0 - axpml[jx]) * (DPxx[IV(iz,ix)] - dedx[jv]);
       p1[IU(iz,kx)] -= t2*v2[IV(iz,kx)] * (dedx[iv] + Hxpml[iv]);
       p1[IU(iz,ix)] -= t2*v2[IV(iz,ix)] * (dedx[jv] + Hxpml[jv]);
   }
   }

   // z direction
   #pragma omp parallel for private (ix, iz, jz, kz)
   for (ix=0; ix<nx ; ix++)
   for (iz=1; iz<nzb; iz++) {
       jz=iz+nzb;
       kz=iz+nz-nzb;
       dpdz[IX(iz,ix)]=0.5*(p2[IU(kz+1,ix)]-p2[IU(kz-1,ix)])/dz;
       dpdz[IX(jz,ix)]=0.5*(p2[IU(iz+1,ix)]-p2[IU(iz-1,ix)])/dz;
   }

   #pragma omp parallel for private (ix, iz)
   for (ix=0; ix<nx;    ix++)
   for (iz=1; iz<nzz-1; iz++)
       Ezpml[IX(iz,ix)]=azpml[iz]*Ezpml[IX(iz,ix)]+(1.0-azpml[iz])*dpdz[IX(iz,ix)];

   #pragma omp parallel for private (ix, iz)
   for (ix=0; ix<nx;    ix++)
   for (iz=1; iz<nzz-1; iz++)
       dedz[IX(iz,ix)]=0.5*(Ezpml[IX(iz+1,ix)]-Ezpml[IX(iz-1,ix)])/dz;

   #pragma omp parallel for private (ix, iz, jz, kz)
   for (ix = 0; ix < nx ; ix++) {
   for (iz = 1; iz < nzb; iz++) {
       jz = iz + nzb;
       kz = iz + nz - nzb;
       int iv = IX(iz,ix);
       int jv = IX(jz,ix);
       Hzpml[iv] = azpml[iz] * Hzpml[iv] + (1.0 - azpml[iz]) * (DPzz[IV(kz,ix)] - dedz[iv]);
       Hzpml[jv] = azpml[jz] * Hzpml[jv] + (1.0 - azpml[jz]) * (DPzz[IV(iz,ix)] - dedz[jv]);
       p1[IU(kz, ix)] -= t2*v2[IV(kz,ix)] * (dedz[iv] + Hzpml[iv]);
       if (ifs == 0)
       p1[IU(iz, ix)] -= t2*v2[IV(iz,ix)] * (dedz[jv] + Hzpml[jv]);
   }
   }

}

void fd2d::getpmlcoef (float *v) {
   int flag=1;
   int ix,iz;
   int jx=nx/2;
   int jz=nz/2;
   int kxr=nx-nxb;
   int kxl=nxb-1;
   int kzb=nz-nzb;
   int kzt=nzb-1;
   int nxx=2*nxb;
   int nzz=2*nzb;
   float p,q,refl=0.000001;

   float disx=1.0*(nxb-1);
   float ux=-1.5*log(refl)/(disx*dx)*dt;
   float disz=1.0*(nzb-1);
   float uz=-1.5*log(refl)/(disz*dz)*dt;

   if (flag==1) {
      for (ix=0; ix<nxb; ix++) {
          p=ix/disx;
          axpml[      ix]=4.0/(2.0+ux*v[IV(jz, kxr+ix)]*p*p)-1.0;
          axpml[nxx-1-ix]=4.0/(2.0+ux*v[IV(jz, kxl-ix)]*p*p)-1.0;
      }
      for (iz=0; iz<nzb; iz++) {
          p=iz/disz;
          azpml[      iz]=4.0/(2.0+uz*v[IV(kzb+iz, jx)]*p*p)-1.0;
          azpml[nzz-1-iz]=4.0/(2.0+uz*v[IV(kzt-iz, jx)]*p*p)-1.0;
      }
   }
   else if (flag==2) {
      for (ix=0; ix<nxb; ix++) {
          p=ix/disx;
          axpml[      ix]=exp(-ux*v[IV(jz, kxr+ix)]*p*p);
          axpml[nxx-1-ix]=exp(-ux*v[IV(jz, kxl-ix)]*p*p);
      }

      for (iz=0; iz<nzb; iz++) {
          p=iz/disz;
          azpml[      iz]=exp(-uz*v[IV(kzb+iz, jx)]*p*p);
          azpml[nzz-1-iz]=exp(-uz*v[IV(kzt-iz, jx)]*p*p);
      }
   }
}

void fd2d::getfdcoef() {
  switch (nor) {
    case 1:
       coef[0] = -2.0;
       coef[1] = 1;
       break;
    case 2:
       coef[0] = -2.55567466;
       coef[1] =  1.37106192;
       coef[2] = -0.09322459;
       break;
    case 3:
       coef[0] = -2.81952122;
       coef[1] =  1.57500756;
       coef[2] = -0.18267338;
       coef[3] =  0.01742643;
       break;
    case 4:
       coef[0] = -2.97399944;
       coef[1] =  1.70507669;
       coef[2] = -0.25861812;
       coef[3] =  0.04577745;
       coef[4] = -0.00523630;
       break;
    case 5:
       coef[0] = -3.05450492;
       coef[1] =  1.77642739;
       coef[2] = -0.30779013;
       coef[3] =  0.07115999;
       coef[4] = -0.01422784;
       coef[5] =  0.00168305;
       break;
    default:
       cout << " Warning: no fd coefficients are given ...\n";
       exit(0);
  }
}

void fd2d::freevar() {
   delete [] p1;
   delete [] p2;
   delete [] p3;

   delete [] DPxx;
   delete [] DPzz;

   delete [] dpdx;
   delete [] dedx;
   delete [] dpdz;
   delete [] dedz;

   delete [] Expml;
   delete [] Hxpml;
   delete [] Ezpml;
   delete [] Hzpml;

   delete [] axpml;
   delete [] azpml;
}
