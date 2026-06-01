#include <cstdio>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include "gsl/gsl_rng.h"

#include "common.hh"
#include "turb_fluid.hh"

// The total number of output points
const int nframes=256;

// The duration to integrate over
const double duration=1;

// The padding factor to apply to the timestep
const double dt_pad=0.2;

int main(int argc,char **argv) {

    // Check for a single command-line argument
    if(argc!=2) {
        fputs("Syntax: ./tracer <num_tracers>\n",stderr);
        return 1;
    }

    // Read the number of tracers from the command line, and check it is in
    // range
    int tra=atoi(argv[1]);
    if(tra<1||tra>256) {
        fputs("Tracer number out of range\n",stderr);
        return 1;
    }

    // Set whether to use the multi-velocity evaluation routine. If the number
    // of tracers is larger than one or two, this routine becomes advantageous.
    bool multi=tra>1;

    // Create the turbulence simulation, and initialize the velocity modes in
    // steady state
    turb_fluid tf(256,256,256,0,160,0,160,0,160,1.,10.);
    tf.init_steady_state();

    // Allocate and initialize tracers
    tf.allocate_vel_table(tra);

    double *x=new double[6*tra],*u=x+3*tra;
    gsl_rng* r=gsl_rng_alloc(gsl_rng_taus2);
    for(int i=0;i<3*tra;i++) x[i]=0.499+0.002*gsl_rng_uniform(r);
    gsl_rng_free(r);

    // Open output file and store initial positions
    FILE *fp=safe_fopen("track.dat","w");
    fputc('0',fp);
    for(double *xp=x;xp<x+3*tra;xp+=3) fprintf(fp," %g %g %g",*xp,xp[1],xp[2]);
    fputc('\n',fp);
    puts("# Output frame 0");

    // Compute a timestep
    double dt=dt_pad*tf.est_max_timestep(),sint=duration/nframes;
    int l=static_cast<int>(sint/dt)+1;
    dt=sint/l;

    // Integrate the modes and output additional snapshots
    double t0=wtime(),t1;
    for(int k=1;k<=nframes;k++) {

        // Perform stochastic integration timesteps
        for(int i=0;i<l;i++) {

            // Integrate the tracers
            if(multi) tf.vel_multi(tra,x,u);
            for(double *xp=x,*up=u;xp<x+3*tra;xp+=3,up+=3) {
                if(!multi) tf.vel(*x,x[1],x[2],*u,u[1],u[2]);
                *xp+=dt*(*up);
                xp[1]+=dt*up[1];
                xp[2]+=dt*up[2];
            }

            // Integrate the velocity field
            tf.step_forward(dt);
        }

        // Output tracer position and print diagnostic information
        t1=wtime();
        fprintf(fp,"%g",k*sint);
        for(double *xp=x;xp<x+3*tra;xp+=3) fprintf(fp," %g %g %g",*xp,xp[1],xp[2]);
        fputc('\n',fp);
        printf("# Output frame %d [%d, %.4g s]\n",k,l,t1-t0);
        t0=t1;
    }
    fclose(fp);

    // Free the dynamically allocated memory
    delete [] x;
}
