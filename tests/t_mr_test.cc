#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.hh"
#include "tf_grid_mr.hh"

// The output directory to write to
const char fn[]="ts2.odr";

// The total number of frames
const int nframes=1024;

// The duration to integrate over
const double duration=12.;

// The padding factor to apply to the timestep
const double dt_pad=0.2;

int main() {

#ifdef FFTW_OMP
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    // Create the output directory
    mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Create the turbulence simulation, and initialize the velocity modes in
    // steady state
    int nmode=64;
    double time_scale=32.;
    double alpha=M_PI*M_PI/(9.*(pow(sqrt(3.)*M_PI,-2/3.)-pow(sqrt(3.)*M_PI*nmode,-2/3.))),
           C=pow(sqrt(3)*M_PI,2/3.)*time_scale;
    turb_fluid_grid tf(nmode,nmode,nmode,0,1,0,1,0,1,C,alpha);
    turb_fluid_grid_mr tf2(nmode,nmode,nmode,0,1,0,1,0,1,C,alpha,10,20);
    turb_fluid_grid tf3(nmode,nmode,nmode,0,1,0,1,0,1,C,alpha);
    tf.init_steady_state();

    // Output the initial cross-section
    char buf[256];
    tf.transform();
    sprintf(buf,"%s/uz.init",fn);
    tf.output_z(buf,0,0);

    double ux,uy,uz,ux2,uy2,uz2;
    tf2.copy_modes(tf);
    tf2.calc_mr_fields();
    for(double T=0;T<20;T+=0.05) {
        tf3.copy_modes(tf);
        tf3.mean_revert(T);
        tf3.transform();
        tf3.lin_interp(0,0,0,ux,uy,uz);
        tf2.lin_interp_mr(T,0,0,0,ux2,uy2,uz2);
        printf("%g %g %g %g %g %g %g\n",T,ux,uy,uz,ux2,uy2,uz2);
    }
}
/*
    tf2.mean_revert(duration);
    tf2.transform();
    sprintf(buf,"%s/uz.mr",fn);
    tf2.output_z(buf,0,0);

    puts("# Output base frames");

    // Compute a timestep
    double dt=dt_pad*tf.est_max_timestep();
    int l=static_cast<int>(duration/dt)+1;
    dt=duration/l;

    // Integrate the modes and output additional snapshots
    for(int i=0;i<3*tf.mno;i++) tf.uu[i]=0;
    double t0=wtime(),t1,t2;
    for(int k=0;k<nframes;k++) {

        tf2.copy_modes(tf);
        for(int i=0;i<l;i++) tf2.step_forward(dt);
        t1=wtime();

        // Output a cross-section
        tf2.transform();
        for(int i=0;i<3*tf.mno;i++) tf.uu[i]+=tf2.uu[i];
        //sprintf(buf,"%s/uz.%d",fn,k);
        //tf2.output_z(buf,0,0);

        // Print diagnostic information
        t2=wtime();
        printf("# Output frame %d [%d, %.4g s, %.4g s]\n",k,l,t1-t0,t2-t1);
        t0=t2;
    }

    double fac=1./nframes;
    for(int i=0;i<3*tf.mno;i++) tf.uu[i]*=fac;
    sprintf(buf,"%s/uz.avg",fn);
    tf.output_z(buf,0,0);
}*/
