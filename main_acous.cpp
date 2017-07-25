#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include "omp.h"
#include "math.h"
#include "array.h"
#include "inout.h"
#include "numpy.h"
#include "para.h"
#include "getpara.h"
#include "acoustic.h"
#include "wavelet.h"
#include "shotdomain.h"
#include "Afd2d.h"
#define NSAFE 10
using namespace std;


int main (int argc, char **argv) {

   int nprocs, myid, namelen;
   char processor_name[MPI_MAX_PROCESSOR_NAME];

   //
   MPI_Init (&argc, &argv);
   MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
   MPI_Comm_rank (MPI_COMM_WORLD, &myid);
   MPI_Get_processor_name (processor_name, &namelen);
   //
   //int nprocs = 1;
   //int myid = 0;

   // *** read parameters from file
   time_t tstart, tend;
   Geometry geo;
   Fwdpar par;

   tstart = time(0);
   par.nproc = nprocs;

   GetPara PARA("fwd.par", &geo, &par);

   // *** openmp
   if (nprocs == 1 && par.nthread > 0) {
       if (par.nthread > omp_get_max_threads()) {
           par.nthread = omp_get_max_threads();
           cout << " --- Max threads = " << par.nthread << endl;
           cout << " --- Warn: number of threads is too large ...\n";
           cout << " --- Warn: only " << par.nthread << " threads are used \n";
       }
       omp_set_num_threads (par.nthread); 
   }
   else {
       #pragma omp parallel
       {
       int nthread = omp_get_num_threads();
       int ithread = omp_get_thread_num();
       if (ithread == 0) par.nthread = nthread;
       }
   }
   if (myid == 0) PARA.PrintPara (par);

   int iso = par.ishot1 - 1;
   int ns = geo.ns;
   int msize = geo.nx * geo.ny * geo.nz;
   float *vp = new float[msize];
   float *vs = new float[msize];
   float *mu = new float[msize];
   float *lam = new float[msize];
   float *rho = new float[msize];
   Shot *src = new Shot[ns];

   // *** read vel/rho and source-receiver configuration
   PARA.GetAcouVel (vp, vs, rho);
   PARA.GetHeader (src, myid);
   MinMaxVel (&par.vmin, &par.vmax, vp, msize);

   // *** modeling
   for (int is = iso + myid; is < ns; is += nprocs) {

        cout << " OOO myid=" << myid << " is=" << is+1 << endl;

        stringstream sidx;
        sidx << setfill('0') << setw(4) << is + 1;
        string fo("shot1/atr_");
        string fo1("shot1/atr_");
        string fo2("shot1/atr_");
        fo = fo + sidx.str() + ".bin";
        fo1 = fo1 + sidx.str() + "d.bin";
        fo2 = fo2 + sidx.str() + "b.bin";

        int iw = 0;
        if (ns == 1) iw = 1;

        int nxo = src[is].nx;
        int nto = src[is].nt;
        int m = nxo * nto;
        float *uo = new float[m];
        float *ud = new float[m];

        forward (uo, vp, rho, src[is], par, geo, iw, is, myid);
        writedata((char *)fo.c_str(), uo, m);

        if (par.iway == 1) {
            if (myid == 0) cout << " *** Only direct waves ... \n";
            MakeDirect (ud, vp, rho, src[is], par, geo, 0, is, myid);
            asub(uo, uo, ud, m);
            writedata((char *)fo1.c_str(), uo, m);
        }

        delete [] uo;
        delete [] ud;
   }

   MPI_Barrier (MPI_COMM_WORLD);

   tend = time(0);
   if (myid == 0) {
      float tcost = tend - tstart;
      if (tcost < 60)
          cout << " Elapsed time: " << tcost << " sec\n";
      else if (tcost < 3600)
          cout << " Elapsed time: " << tcost/60 << " min\n";
      else
          cout << " Elapsed time: " << tcost/3600 << " hour\n";
   }

   delete [] vp;
   delete [] vs;
   delete [] mu;
   delete [] lam;
   delete [] rho;
   delete [] src;

   MPI_Barrier (MPI_COMM_WORLD);
   MPI_Finalize ();
}



void forward (float *uo, float *v, float *rho, Shot src, Fwdpar par, Geometry geo, int iw, int is, int myid) {

        // output size
        int nxo = src.nx;
        int nzo = src.nz;
        int nto = src.nt;
        float dxo = src.dx;
        float dzo = src.dz;
        float dto = src.dt;

        float fp = par.fp;
        float fm = fp * 2.0;
        float vmin = par.vmin;
        float vmax = par.vmax;

        // computation grid
        Shot srcfd;
        GetShotDomain (&srcfd, src, fp, fm, vmin, vmax, par.nxor, par.ntor);

        // computation size 
        int nxb = par.nxpml + NSAFE;
        int nzb = par.nzpml + NSAFE;
        int nxg = srcfd.nx + 2 * nxb;
        int nzg = srcfd.nz + 2 * nzb;
        int ntc = srcfd.nt;
        float dx = srcfd.dx;
        float dz = srcfd.dz;
        float dt = srcfd.dt;

        if (is == 0) {
            cout << " nx = " << nxg << " (" << srcfd.nx << ", " << nxb << ")  " 
                 << " nz = " << nzg << " (" << srcfd.nz << ", " << nzb << ")  " 
                 << " nt = " << ntc << " (" << srcfd.nt << ")\n"; 
            cout << " dx = " << dx << " dz = " << dz << " dt = " << dt << endl;
            //for (int ix = 0; ix < srcfd.nx; ix ++)
            //     cout << " is = " << is << " ix = " << ix << " iz = " << srcfd.izr[ix] << endl;
        }

        // local coordinate (with 1st receiver as origin)
        float fxl = - nxb * srcfd.dx;
        float fzl = - nzb * srcfd.dz;
 
        // global coordinate (whole model)
        float fx = fxl + srcfd.xg;
        float fz = fzl;

        // free surface index
        int izf = -1;
        if (par.ifs > 0) izf = nzb;

        int ixs = srcfd.ixs + nxb;
        int izs = srcfd.izs + nzb;
        int *izg = new int[nxg];
        for (int ix = 0; ix < srcfd.nx; ix ++)
             izg[ix + nxb] = srcfd.izr[ix] + nzb;
        for (int ix = 0; ix < nxb; ix ++) {
             izg[ix] = izg[nxb];
             izg[nxg - 1 - ix] = izg[nxg - 1 - nxb];
        }
            
        // prepare wavelet
        Wavelet wvl(par.fp, par.td, dt, ntc, par.iwvl);
        wvl.shift();
        wvl.mkwvl(1.0);
        //wvl.time_derivative(1.0);
        if (myid == 0 && is == 0) writedata("wavelet.bin", wvl.wv, ntc);
        float ftl = -wvl.nshift * wvl.dt;

        // prepare source
        float *uw = new float[ntc * nxg];
        memset(uw, 0, nxg * ntc * sizeof(float));
        if (par.iplw == 0)
            memcpy(&uw[ixs*ntc], wvl.wv, ntc * sizeof(float));
        else
            StackSource(uw, wvl.wv, src.p, src.xo - fxl, dx, dt, nxg, ntc);

        // prepare media parameters
        float *rho1 = new float[nxg * nzg];
        float *v1 = new float[nxg * nzg];

        interp2d (rho1, fxl, fzl, dx, dz, nxg, nzg, rho, geo.bx, 0.0, geo.dx, geo.dz, geo.nx, geo.nz); 
        interp2d (  v1, fxl, fzl, dx, dz, nxg, nzg,   v, geo.bx, 0.0, geo.dx, geo.dz, geo.nx, geo.nz); 
        
        // FD simulation
        fd2d fdshot (nxg, nzg, ntc, dx, dz, dt, par.nxpml, par.nzpml, par.nxor);
        fdshot.snapset (dxo, 0, dzo, 10*dt, nxo, 1, nzo, nto);
        fdshot.initabc (izf, par.fp, v1);
        fdshot.forward (uw, v1, ixs, izs, izg, iw);

        interp2d ( uo, 0.0, 0.0, dxo, dto, nxo, nto, uw, fxl, ftl, dx, dt, nxg, ntc);

        delete [] izg;
        delete [] rho1;
        delete [] v1;
        delete [] uw;
        wvl.delwvl();
}




void MakeDirect (float *ud, float *vp, float *rho, Shot src, Fwdpar par, Geometry geo, int iw, int is, int myid) {
    int m, m0, nv = geo.nx * geo.nz;
    float *rho1 = new float[nv];
    float *vp1 = new float[nv];

    memcpy (rho1, rho, nv * sizeof(float));
    memcpy ( vp1,  vp, nv * sizeof(float));

    for (int ix = 0; ix < geo.nx; ix ++) {
         m = ix * geo.nz;
         m0 = m + geo.kzmin;
         for (int iz = geo.kzmin + 1; iz < geo.nz; iz ++) {
              rho1[m + iz] = rho[m0];
               vp1[m + iz] =  vp[m0];
    }
    }

    forward (ud, vp1, rho1, src, par, geo, iw, is, myid);

    delete [] rho1;
    delete [] vp1;
}

void MinMaxVel (float *vmin, float *vmax, float *v, int n) {
  float v1 = v[0];
  float v2 = v[0];
  for (int i = 1; i < n; i ++) {
       if (v1 < v[i]) v1 = v[i];
       if (v2 > v[i]) v2 = v[i];
  }

  *vmax = v1;
  *vmin = v2;
}

void StackSource (float *u, float *wv, float p, float xo, float dx, float dt, int nx, int nt) {
  int its;
  int nxh = 0.1 * nx;
  float *xtap = new float[nx];
  hamfunc(xtap, 0.5, nx, 0, nxh, nx - 1 - nxh, nx - 1, 0);

  for (int ix = 0; ix < nx; ix ++) {
       its = p * (ix * dx - xo) / dt + 0.5;
       int jt1 = its;
       int jt2 = its + nt - 1;
       if (jt1 < 0) jt1 = 0;
       if (jt2 >= nt) jt2 = nt - 1;
       int it1 = jt1 - its;
       for (int it = it1, jt = jt1; jt <= jt2; it ++, jt ++)
            u[ix*nt+jt] += xtap[ix] * wv[it];
  }
  delete [] xtap;
}
