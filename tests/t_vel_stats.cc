#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.hh"
#include "tf_grid.hh"

// The total number of frames
const int nframes=50000;

// The duration to integrate over
const double duration=1000.;

// The padding factor to apply to the timestep
const double dt_pad=0.2;

int main() {

#ifdef FFTW_OMP
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    // Create the turbulence simulation, and initialize the velocity modes in
    // steady state
    int nmode=32;
    double time_scale=32.;
    double alpha=M_PI*M_PI/(9.*(pow(sqrt(3.)*M_PI,-2/3.)-pow(sqrt(3.)*M_PI*nmode,-2/3.))),
           C=pow(sqrt(3)*M_PI,2/3.)*time_scale;
    turb_fluid_grid tf(nmode,nmode,nmode,0,50,0,50,0,50,C,alpha);
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
    printf("%g\n",tf.mean_rms());
  /*  tf.output_x("diagx",0,2);
    tf.output_y("diagy",0,2);
    tf.output_z("diagz",0,2);*/
    for(int k=1;k<=nframes;k++) {

        // Perform stochastic integration timesteps
        for(int i=0;i<1;i++) tf.init_steady_state();
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
        if(k%1000==0) printf("# Output frame %d [%d, %.4g s, %.4g s]\n",k,l,t1-t0,t2-t1);
        t0=t2;
    }
    fclose(fp);
    printf("Long term average: %g\n",srms/(nframes+1));
}
