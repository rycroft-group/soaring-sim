#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.hh"
#include "tf_grid.hh"

// The total number of frames
const int nframes=2000;

// The number of Fourier modes to use in the turbulent wind field
const int nmode=64;

// The side length of the cube for the turbulent wind field
const double L=50;

// Parameters controlling the wind field mode amplitude and temporal
// fluctuations. (Based on tau_c=64 s, and 1 m/s RMS per component.)
const double alpha=0.00913104,C=63.3783,Cinv=1/C;

// The duration over which to simulate. Corresponds to 512 s.
const double dur_phys=8000,
             duration=dur_phys*0.9902853;

// The padding factor to apply to the timestep
const double dt_pad=0.2;

int main() {

    // Set up multi-threaded FFTW computations, if available
#ifdef _OPENMP
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    // Create the turbulence simulation, and initialize the velocity modes in
    // steady state
    turb_fluid_grid tf(nmode,nmode,nmode,0,L,0,L,0,L,Cinv,alpha);
    double ubar,vbar,wbar,urms,vrms,wrms,srms=0,trms;
    tf.init_steady_state();

    // Output the initial cross-section
    FILE *fp=safe_fopen("vstats.dat","w");
    tf.transform();
    puts("# Output frame 0");
    tf.grid_vel_stats(ubar,vbar,wbar,urms,vrms,wrms);
    trms=urms*urms+vrms*vrms+wrms*wrms;
    srms+=trms;
    fprintf(fp,"0 0 %g %g %g %g %g %g %g\n",ubar,vbar,wbar,urms,vrms,wrms,sqrt(trms));

    // Compute a timestep
    double dt=dt_pad*tf.est_max_timestep(),sint=duration/nframes;
    int l=static_cast<int>(sint/dt)+1;
    dt=sint/l;

    // Integrate the modes and output additional snapshots
    double t0=wtime(),t1,t2;
    for(int k=1;k<=nframes;k++) {

        // Perform stochastic integration timesteps
        for(int i=0;i<l;i++) tf.step_forward(dt);
        t1=wtime();

        // Output a cross-section
        tf.transform();
        tf.grid_vel_stats(ubar,vbar,wbar,urms,vrms,wrms);
        trms=sqrt(urms*urms+vrms*vrms+wrms*wrms);
        srms+=trms;
        fprintf(fp,"%d %g %g %g %g %g %g %g %g\n",k,sint*k,ubar,vbar,wbar,urms,vrms,wrms,
                trms);

        // Print diagnostic information
        t2=wtime();
        if(k%100==0) printf("# Output frame %d [%d, %.4g s, %.4g s]\n",k,l,t1-t0,t2-t1);
        t0=t2;
    }
    fclose(fp);
    printf("# Long term average: %g\n",srms/(nframes+1));
}
