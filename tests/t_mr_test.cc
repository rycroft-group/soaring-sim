#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.hh"
#include "tf_grid_mr.hh"

// The output directory to write to
const char fn[]="t_mr_test.odr";

// The total number of instances to simulate
const int n_inst=1024;

// The duration to integrate over. Corresponds to 2 s.
const double dur_phys=8.,
             duration=dur_phys*0.9902853;

// The padding factor to apply to the timestep
const double dt_pad=0.2;

// The number of Fourier modes to use in the turbulent wind field
const int nmode=64;

// The side length of the cube for the turbulent wind field
const double L=50;

// Parameters controlling the wind field mode amplitude and temporal
// fluctuations. (Based on tau_c=64 s, and 1 m/s RMS per component.)
const double alpha=0.00913104,C=63.3783,Cinv=1/C;

int main() {

    // Set up multi-threaded FFTW computations, if available
#ifdef _OPENMP
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    // Create the output directory
    mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Create the turbulent wind field, and initialize the velocity modes in
    // steady state
    turb_fluid_grid tf(nmode,nmode,nmode,0,1,0,1,0,1,Cinv,alpha);
    tf.init_steady_state();

    // Output the initial cross-section
    char buf[256];
    tf.transform();
    sprintf(buf,"%s/uz.init",fn);
    tf.output_z(buf,0,0);

    // Create a copy of the turbulent wind field, and mean-revert the modes
    turb_fluid_grid tf2(nmode,nmode,nmode,0,1,0,1,0,1,Cinv,alpha);
    tf2.copy_modes(tf);
    tf2.mean_revert(duration);
    tf2.transform();
    sprintf(buf,"%s/uz.mr",fn);
    tf2.output_z(buf,0,0);
    puts("# Output base frames");

    // Compute a timestep
    double dt=dt_pad*tf.est_max_timestep();
    int l=static_cast<int>(duration/dt)+1;
    dt=duration/l;

    // Simulate many instances of the modes evolving stochastically
    for(int i=0;i<3*tf.mno;i++) tf.uu[i]=0;
    double t0=wtime(),t1,t2;
    for(int k=0;k<n_inst;k++) {

        // Initialize the turbulent wind field with the reference set of modes
        // and time-integrate them forward
        tf2.copy_modes(tf);
        for(int i=0;i<l;i++) tf2.step_forward(dt);
        t1=wtime();

        // Store the wind field contribution to the average
        tf2.transform();
        for(int i=0;i<3*tf.mno;i++) tf.uu[i]+=tf2.uu[i];

        // Output a cross-section
        sprintf(buf,"%s/uz.%d",fn,k);
        tf2.output_z(buf,0,0);

        // Print diagnostic information
        t2=wtime();
        printf("# Output instance %d [%d, %.4g s, %.4g s]\n",k,l,t1-t0,t2-t1);
        t0=t2;
    }

    // Output the wind field averaged across all of the instances, which should
    // be a close match to the deterministically mean-reverted field
    double fac=1./n_inst;
    for(int i=0;i<3*tf.mno;i++) tf.uu[i]*=fac;
    sprintf(buf,"%s/uz.avg",fn);
    tf.output_z(buf,0,0);
}
