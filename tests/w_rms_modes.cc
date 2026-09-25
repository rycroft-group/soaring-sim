#include <cstdio>
#include <cmath>

#include "fileinfo.hh"

int nmodes[]={10,12,16,20,24,32,40,48,64,80,96,128,160};

int main() {
    fileinfo fi("sims/lfd64small.cfg");
    for(int k=0;k<13;k++) {
        fi.nx=fi.ny=fi.nz=nmodes[k];
        fi.calculate_wind_param(true);
        printf("%d %.7g\n",nmodes[k],fi.w_rms*fi.v_phys/sqrt(3.));
    }
}
