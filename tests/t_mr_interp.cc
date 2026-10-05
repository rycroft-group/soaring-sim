#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.hh"
#include "tf_grid_mr.hh"

// The number of Fourier modes to use in the turbulent wind field
const int nmode=64;

// The side length of the cube for the turbulent wind field
const double L=50;

// Parameters controlling the wind field mode amplitude and temporal
// fluctuations. (Based on tau_c=64 s, and 1 m/s RMS per component.)
const double alpha=0.00913104,C=63.3783,Cinv=1/C;

// The number of Hermite interpolation segments to use in time
const int hseg=10;

// The duration over which to set up the mean-reversion interpolation field.
// Corresponds to 50 s.
const double Tmr_phys=50,
             Tmr=Tmr_phys*0.9902853;

// The number of sample intervals to divide the duration into, and the corresponding
// interval spacing
const int niv=256;
const double h=Tmr/niv;

int main() {

    // Set up multi-threaded FFTW computations, if available
#ifdef FFTW_OMP
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    // Create the turbulence simulation, and initialize the velocity modes in
    // steady state
    turb_fluid_grid tf(nmode,nmode,nmode,0,L,0,L,0,L,Cinv,alpha);
    turb_fluid_grid_mr tf2(nmode,nmode,nmode,0,L,0,L,0,L,Cinv,alpha,hseg,Tmr);
    turb_fluid_grid tf3(nmode,nmode,nmode,0,L,0,L,0,L,Cinv,alpha);
    tf.init_steady_state();

    // Calculate the mean-reversion interpolation table
    tf2.copy_modes(tf);
    tf2.transform();
    tf2.calc_mr_fields();

    // For a range of values in time, check that the interpolated value of
    // velocity at one point matches a direct computation.
    FILE *fp=safe_fopen("t_mr_interp.dat","w");
    double T,ux,uy,uz,ux2,uy2,uz2;
    for(int i=0;i<=niv;i++) {
        T=i*h;

        // Explicitly mean-revert the modes to perform
        tf3.copy_modes(tf);
        tf3.mean_revert(T);
        tf3.transform();
        tf3.lin_interp(0,0,0,ux,uy,uz);

        // Compute the interpolated velocity and output both
        tf2.lin_interp_mr(T,0,0,0,ux2,uy2,uz2);
        fprintf(fp,"%g %g %g %g %g %g %g\n",T,ux,uy,uz,ux2,uy2,uz2);
    }
    fclose(fp);
}
