#include <cstdio>

#include "common.hh"
#include "turb_fluid.hh"

// Number of histogram bins
const int nbin=200;

// The total number of histograms
const int nhist=100;

// The duration to integrate over
const double duration=100;

// The padding factor to apply to the timestep
const double dt_pad=0.2;

int main(int argc,char **argv) {

    if(argc!=2) {
        fputs("Syntax: ./m_histogram <modes>\n",stderr);
        return 1;
    }

    int nmode=atoi(argv[1]);
    if(nmode<=0||nmode>4096) {
        fputs("Mode number out of range\n",stderr);
        return 1;
    }

    // Create the turbulence simulation, and initialize the velocity modes in
    // steady state
    double alpha=M_PI*M_PI/(9.*(pow(sqrt(3.)*M_PI,-2/3.)-pow(sqrt(3.)*M_PI*nmode,-2/3.)));
    turb_fluid tf(nmode,nmode,nmode,0,10,0,10,0,10,1.,alpha);
    tf.init_steady_state();

    // Compute the histogram for the initial state
    double *hi=new double[nbin*(nhist+1)],hmax;
    tf.histogram(hi,nbin,hmax);
    puts("# Compute frame 0");

    // Compute the timestep
    double dt=dt_pad*tf.est_max_timestep(),sint=duration/nhist,
           ubar,vbar,wbar,uvar,vvar,wvar;
    int l=static_cast<int>(sint/dt)+1;
    dt=sint/l;
    printf("dt: %g\n",dt);

    // Integrate the Fourier modes
    double t0=wtime(),t1,t2;
    for(int i=1;i<=nhist;i++) {

        // Perform stochastic integration timesteps
        for(int k=0;k<l;k++) tf.step_forward(dt);
        t1=wtime();

        // Compute the histogram
        tf.histogram(hi+nbin*i,nbin,hmax);
//        tf.velocity_stats(1000,ubar,vbar,wbar,uvar,vvar,wvar);

        // Print diagnostic information
        t2=wtime();
//        printf("%g %g %g %g %g %g %g\n",ubar,vbar,wbar,uvar,vvar,wvar,uvar+vvar+wvar);
        printf("# Compute frame %d [%d, %.4g s, %.4g s]\n",i,l,t1-t0,t2-t1);
        t0=t2;
    }

    // Output the histograms
    char buf[256];
    sprintf(buf,"mhist%d.dat",nmode);
    FILE *fp=safe_fopen(buf,"w");
    for(int j=0;j<nbin;j++) {
        fprintf(fp,"%g",(j+0.5)/nbin*hmax);
        for(int i=0;i<=nhist;i++) fprintf(fp," %g",hi[nbin*i+j]);
        fputc('\n',fp);
    }

    // Close the file and remove the dynamically allocated memory
    fclose(fp);
    delete [] hi;
}
