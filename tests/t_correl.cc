#include <cstdio>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.hh"
#include "tf_grid.hh"

// The duration to integrate over
const double duration=0.1;

// The padding factor to apply to the timestep
const double dt_pad=0.2;

// The number of instances
const int ins=81;

int main(int argc,char **argv) {

    if(argc!=3) {
        fputs("Syntax: ./t_correl <modes> <max_radius>\n",stderr);
        return 1;
    }

    int nmode=atoi(argv[1]);
    if(nmode<=0||nmode>4096) {
        fputs("Mode number out of range\n",stderr);
        return 1;
    }

    const int nbin=201,nsamp=4096;
    double mrad=atof(argv[2]),h=mrad/(nbin-1);

#ifdef FFTW_OMP
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    // Create the turbulence simulation, and initialize the velocity modes in
    // steady state
    double alpha=1,C=1;
    turb_fluid_grid tf(nmode,nmode,nmode,0,50,0,50,0,50,C,alpha);
    tf.init_steady_state();

    double w[ins][nbin],*samp=new double[6*nsamp];

    // Compute a timestep
    double dt=dt_pad*tf.est_max_timestep();
    int l=static_cast<int>(duration/dt)+1;
    dt=duration/l;

    // Perform stochastic integration timesteps
    //printf("dt=%g, steps=%d\n",dt,l);
    tf.correl_init(samp,nsamp);

    for(int o=0;o<ins;o++) {
        for(int i=0;i<l;i++) tf.step_forward(dt);

        // Compute correlations at a later time
        tf.correl_function(w[o],nbin,mrad,samp,nsamp);

        printf("%d %g",o,w[o][0]);
        for(int k=1;k<6;k++) printf(" %g",w[o][k]);
        putchar('\n');
    }
    delete [] samp;

    char buf[256];
    sprintf(buf,"tcor_%d_tch2_r%s.dat",nmode,argv[2]);
    FILE *fp=safe_fopen(buf,"w");
    for(int q=0;q<nbin;q++) {
        double s=w[0][q];
        for(int o=1;o<ins;o++) s+=w[o][q];

        fprintf(fp,"%g %g",(q+0.5)*h,s/ins);
        for(int o=0;o<ins;o++) fprintf(fp," %g",w[o][q]);

        fputc('\n',fp);
    }
    fclose(fp);

}
