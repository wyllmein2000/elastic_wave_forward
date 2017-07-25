#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "getpara.h"
#include "inout.h"
using namespace std;

GetPara::GetPara (char *filename, Geometry *geo, Fwdpar *par) {
   geo->ny = 1;
   geo->by = 0;
   geo->dy = 0;
   geo->kzmin = 0;
   geo->kzmax = 9999;

   par->ifs = 0;
   par->nxpml = 10;
   par->nzpml = 10;
   par->nxor = 6;
   par->ntor = 2;
   par->dim = 2;

   par->iwvl = 1;
   par->isrc = 0;

   readpar (filename, this->dataf, this->headf, this->vf1, this->vf2, this->vf3, geo, par);

   geo->ns = par->ishot2;

   if (par->iplw < 0 || par->iplw > 1) {
       cout << " Error: iplw = " << par->iplw << endl;
       exit(0);
   }
   if (par->ishot1 < 1 || par->ishot2 < par->ishot1) {
       cout << " Error: shot index = (" << par->ishot1 << " ," << par->ishot2 << endl;
       exit(0);
   }


   this->nx = geo->nx;
   this->ny = geo->ny;
   this->nz = geo->nz;
   this->nt = geo->nt;
   this->ns = geo->ns;
   this->dx = geo->dx;
   this->dy = geo->dy;
   this->dz = geo->dz;
   this->dt = geo->dt;
   this->bx = geo->bx;
   this->by = geo->by;
   this->bz = geo->bz;
}

void GetPara::GetAcouVel (float *vp, float *vs, float *rho) {
   int v_size = nx * ny * nz;
   readdata(vf1, vp, v_size);
   for (int i = 0; i < v_size; i ++) {
        vs[i] = vp[i] / sqrt(3);
        rho[i] = 3000.0;
   }
   //readdata(vf3, v3, v_size);
}

void GetPara::GetElasVel (float *v1, float *v2, float *v3) {
   int v_size = nx * ny * nz;
   readdata(vf1, v1, v_size);
   readdata(vf2, v2, v_size);
   readdata(vf3, v3, v_size);
}

/* vp, vs: m/s; rho: kg/m3; lam, mu: MPa  */
void GetPara::MapVpToLam (float *lam, float *mu, float *rho, float *vp, float *vs) {
   int v_size = nx * ny * nz;
   for (int i = 0; i < v_size; i ++) {
        rho[i] = rho[i] * 1e-6;
         mu[i] = rho[i] * vs[i] * vs[i];
        lam[i] = rho[i] * vp[i] * vp[i] - 2 * mu[i];
   }
   //writedata("junkrho.bin",rho,v_size);
   //writedata("junklam.bin",lam,v_size);
   //writedata("junkmu.bin",mu,v_size);
}

void GetPara::GetHeader (Shot *src, int myid) {
   int is, nshot;
   int *kelev;
   float xm = bx + (nx - 1) * dx;

   ifstream infile(headf);
   if (!infile) {
       cerr << " Error: uable to open input file: " << headf << endl;
       exit(0);
   }
   infile >> nshot;
   if (nshot < ns) {
       cerr << " Error: incomplete header infor ...\n";
       exit(0);
   }
   if (myid == 0) {
       cout << endl;
       cout << " *** Shots infor ...\n";
   }
   for (is = 0; is < ns; is++) {
       infile >> src[is].xs >> src[is].zs >> src[is].xg >> src[is].zg >> src[is].nx;
       src[is].ishot = is;
       src[is].topo = 0;
       src[is].ixs = (src[is].xs - src[is].xg) / dx + 0.5;
       src[is].ixg = 0;
       src[is].izs = (src[is].zs - bz)/ dz + 0.5;
       src[is].izg = (src[is].zg - bz) / dz + 0.5;
       src[is].nz = nz;
       src[is].nt = nt;
       src[is].dx = dx;
       src[is].dy = dy;
       src[is].dz = dz;
       src[is].dt = dt;

       src[is].p = src[is].xs;
       src[is].xo = bx;
       if (src[is].p <= 0) src[is].xo = xm;

       src[is].zr = new float[src[is].nx];
       src[is].izr = new int[src[is].nx];
       for (int ix = 0; ix < src[is].nx; ix++) {
            src[is].zr[ix] = src[is].zg;
            src[is].izr[ix] = src[is].izg;
       }

       if (myid == 0 && is % 20 == 0)
           cout << "     ishot=" << is << ", ixs=" << src[is].ixs << ", izs=" << src[is].izs << ", nxg=" << src[is].nx << endl;
   }
   infile.close();

   if (src[0].zs < 0) {
       kelev = new int[nx];
       readdata("gelevation.dat", kelev, nx); 
       for (int is = 0; is < ns; is ++) {
            int jx = (src[is].xs - bx) / dx;
            src[is].izs = kelev[jx];
            src[is].zs  = kelev[jx] * dz;
            src[is].izg = src[is].izs;
            src[is].zg  = src[is].zs;
            for (int ix = 0; ix < src[is].nx; ix ++) {
                 int jx = (src[is].xg + ix * dx - bx) / dx;
                 src[is].izr[ix] = kelev[jx];
                 src[is].zr[ix]  = kelev[jx] * dz;
                 if (src[is].izg > src[is].izr[ix]) {
                     src[is].izg = src[is].izr[ix];
                     src[is].zg  = src[is].zr[ix];
                 }
            }
       }
       delete [] kelev;
   }

}

void GetPara::PrintPara (Fwdpar par) {
   cout << endl;
   cout << " *** Input Parameters  ***" << endl;
   cout << "     model size:" << endl;
   cout << "     model dim =" << par.dim << endl;
   cout << "     nx = " << nx << ", bx = " << bx << ", dx = " << dx << endl;
   cout << "     ny = " << ny << ", by = " << by << ", dy = " << dy << endl;
   cout << "     nz = " << nz << ", bz = " << bz << ", dz = " << dz << endl;
   cout << "     \n";
   cout << "     Data infor:" << endl;
   cout << "      source type: isrc = " << par.isrc << " (=0, explosive; =1, z-force)\n";
   cout << "     wavelet type: iwvl = " << par.iwvl << " (=-1, read; =0, impulse; =1, ricker)\n";
   cout << "      dominant freq: fp = " << par.fp << " Hz\n";
   cout << "      gather Type: iplw = " << par.iplw << " (=0, shot; =1, plane-wave)\n";
   cout << "      Shot Index: ishot = (" << par.ishot1 << ", " << par.ishot2 << ")\n";
   cout << "     \n";
   cout << "     computing ... " << endl;
   cout << "           PML SIZE: nxpml = " << par.nxpml << ", nzpml = " << par.nzpml << endl;
   cout << "         Free Surface: ifs = " << par.ifs << " (=0, no; =1, yes)\n";
   cout << "      FD Order: nxor, ntor = (" << par.nxor << ", " << par.ntor << ")\n";
   cout << "        Number of processors: " << par.nproc << endl;
   cout << "     Num of threads per proc: " << par.nthread << endl;
   cout << endl;
   cout << " *** Output Parameters  ***" << endl;
   cout << "     data size:" << endl;
   cout << "     nt = " << nt << ", dt = " << dt << " sec\n";
}

void GetPara::readpar (char *filename, char *dfile, char *hfile, char *vf1, char *vf2, char *vf3, Geometry *geomd, Fwdpar *par) {
     FILE *fp;
     char *sp, str[132], atom[132];

     fp = fopen(filename, "r");
     if (fp == NULL)
        cout << "Input parameter file: " << filename << ", open error!\n";

     while (fgets(str, 132, fp)) {
        if (str[0] == '#' || str[0] == '!') continue;
        sp = str;

        if (!next_atom( &sp, atom)) continue;
        if (!strcmp( atom, "TRACE")) {
            next_atom( &sp, atom);
            sprintf(dfile, "%s", atom);
        } else if (!strcmp( atom, "HEADFILE")) {
            next_atom( &sp, atom);
            sprintf(hfile, "%s", atom);
        } else if (!strcmp( atom, "VPFILE")) {
            next_atom( &sp, atom);
            sprintf(vf1, "%s", atom);
        } else if (!strcmp( atom, "VSFILE")) {
            next_atom( &sp, atom);
            sprintf(vf2, "%s", atom);
        } else if (!strcmp( atom, "RHOFILE")) {
            next_atom( &sp, atom);
            sprintf(vf3, "%s", atom);
        } else if (!strcmp( atom, "FORWARD")) {
            next_atom( &sp, atom);
            par->iway = atoi(atom);
        } else if (!strcmp( atom, "RANGE")) {
            next_atom( &sp, atom);
            geomd->kzmin = atoi(atom);
            next_atom( &sp, atom);
            geomd->kzmax = atoi(atom);
        } else if (!strcmp( atom, "VELZ")) {
            next_atom( &sp, atom);
            geomd->bz = atof(atom);
            next_atom( &sp, atom);
            geomd->dz = atof(atom);
            next_atom( &sp, atom);
            geomd->nz = atoi(atom);
        } else if (!strcmp( atom, "VELY")) {
            next_atom( &sp, atom);
            geomd->by = atof(atom);
            next_atom( &sp, atom);
            geomd->dy = atof(atom);
            next_atom( &sp, atom);
            geomd->ny = atoi(atom);
        } else if (!strcmp( atom, "VELX")) {
            next_atom( &sp, atom);
            geomd->bx = atof(atom);
            next_atom( &sp, atom);
            geomd->dx = atof(atom);
            next_atom( &sp, atom);
            geomd->nx = atoi(atom);
        } else if (!strcmp( atom, "TIME")) {
            next_atom( &sp, atom);
            geomd->dt = atof(atom);
            next_atom( &sp, atom);
            geomd->nt = atoi(atom);
        } else if (!strcmp( atom, "WAVELET")) {
            next_atom( &sp, atom);
            par->fp = atof(atom);
            next_atom( &sp, atom);
            par->td = atof(atom);
            next_atom( &sp, atom);
            par->iwvl = atoi(atom);
            next_atom( &sp, atom);
            par->isrc = atoi(atom);
        } else if (!strcmp( atom, "FREQUENCY")) {
            next_atom( &sp, atom);
            par->f1 = atof(atom);
            next_atom( &sp, atom);
            par->f2 = atof(atom);
            next_atom( &sp, atom);
            par->f3 = atof(atom);
            next_atom( &sp, atom);
            par->f4 = atof(atom);
        } else if (!strcmp( atom, "SHOTS")) {
            next_atom( &sp, atom);
            par->iplw = atoi(atom);
            next_atom( &sp, atom);
            par->ishot1 = atoi(atom);
            next_atom( &sp, atom);
            par->ishot2 = atoi(atom);
        } else if (!strcmp( atom, "PML")) {
            next_atom( &sp, atom);
            par->nzpml = atoi(atom);
            next_atom( &sp, atom);
            par->nxpml = atoi(atom);
            next_atom( &sp, atom);
            par->ifs = atoi(atom);
        } else if (!strcmp( atom, "ORDER")) {
            next_atom( &sp, atom);
            par->nxor = atoi(atom);
            next_atom( &sp, atom);
            par->ntor = atoi(atom);
        } else if (!strcmp( atom, "nthread")) {
            next_atom( &sp, atom);
            par->nthread = atoi(atom);
        }
     }
     

}

//===================================================================
/*
 *  * next_atom -- gets next field in a GEOMON free-format card. Returns 0
 *  * if *startp points to end of card (NULL). The field is copied into String
 *  * atom, and *startp is updated to point to the start of the next field.
 *  * Fields are delimited by commas alone; leading and trailing blanks are
 *  * ignored.
 *  */
int GetPara::next_atom( char **startp, char* atom)
        /*String *startp;*/                 /* Beginning of next field in card */
        /*String atom;*/                    /* Returned null-terminated field */
{
        char *sp, *ap, *bp;
        int got_it = 0,count;
        if (!*(sp = *startp)) return (0);

        count = 0;
        for (ap = bp = atom; *sp && !got_it; sp++) {
                switch (*sp) {
                case ' ':
                        if(count > 0){
                           got_it++;
                           break;
                        }
                case '\t':
                case '\n':
                        //if(count == 0)return(0);
                        if (ap == atom) continue;
                        *ap++ = *sp;
                        break;
                case ',':
                        got_it++;
                        break;
                default:
                        //*ap++ = (char)toupper( (int)(*sp));
                        *ap++ = (char)(*sp);
                        count += 1;
                        bp = ap;
                }
        }
        *bp = 0;
        *startp = sp;
        return (1);
}
//===================================================================
