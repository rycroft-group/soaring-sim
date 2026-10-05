#include <cstdio>
#include <cmath>

#include "fileinfo.hh"

// The number of modes to consider
int nmodes[]={10,12,16,20,24,32,40,48,64,80,96,128,160};

int main() {

    // Read in the parameter set from the tau_c=64 s example simulation
    fileinfo fi("../sims/lfd64small.cfg");

    // Loop over the different numbers of modes
    for(int k=0;k<13;k++) {

        // Override the mode number values in the fileinfo class, and use it to
        // calculate the RMS wind speed.
        fi.nx=fi.ny=fi.nz=nmodes[k];
        fi.calculate_wind_param(true);

        // Print the number of modes and the RMS wind component velocity in
        // m/s. The factor of sqrt(3) is due to the conversion from the wind
        // speed to a single velocity component.
        printf("%d %.7g\n",nmodes[k],fi.w_rms*fi.v_phys/sqrt(3.));
    }
}
