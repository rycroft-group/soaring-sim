#include <cstdio>
#include <cmath>

#include "common.hh"
#include "en_spec.hh"
#include "turb_fluid.hh"

// The total number of histograms during time integration
const int nhist=8;

// Physical scales in SI units for scaling
const double g_phys=9.80665,l_phys=10,
             t_phys=sqrt(l_phys/g_phys),
             v2_phys=l_phys*g_phys,v_phys=sqrt(v2_phys);

// The duration to integrate over. Corresponds to 800 s.
const double dur_phys=800,duration=dur_phys/t_phys;

// The padding factor to apply to the timestep
const double dt_pad=0.075;

// The side length of the cube for the turbulent wind field
const double L=50;

// Parameters controlling the wind field mode amplitude and temporal
// fluctuations. (Based on tau_c=32 s, and 1 m/s RMS per component at N=64.)
const double alpha=0.00913104,C=32/t_phys,Cinv=1./C;

int main(int argc,char **argv) {

    // Print a syntax message if the number of modes hasn't been supplied
    if(argc!=2) {
        fputs("Syntax: ./en_spectrum <modes>\n",stderr);
        return 1;
    }

    // Check that the number of modes is in a reasonable range
    int nmode=atoi(argv[1]);
    if(nmode<=0||nmode>4096) {
        fputs("Mode number out of range\n",stderr);
        return 1;
    }

    // Create the turbulence simulation, and initialize the velocity modes in
    // steady state
    turb_fluid tf(nmode,nmode,nmode,0,L,0,L,0,L,Cinv,alpha);
    tf.init_steady_state();

    // Set up the nonlinear mapping in the energy spectrum computation
    en_spec_param es(72,15);
    tf.setup_es_param(es);

    // Compute the histogram for the initial state
    en_spec_data *ed=new en_spec_data[(nhist+1)*es.nbin];
    double E_integ=tf.energy_spectrum(es,ed);
    printf("# Compute frame 0 {%g m^2/s^2}\n",v2_phys*E_integ);

    // Compute the timestep
    double dt=dt_pad*tf.est_max_timestep(),sint=duration/nhist;
    int l=static_cast<int>(sint/dt)+1;
    dt=sint/l;

    // Integrate the Fourier modes
    double t0=wtime(),t1,t2;
    for(int i=1;i<=nhist;i++) {

        // Perform stochastic integration timesteps
        for(int k=0;k<l;k++) tf.step_forward(dt);
        t1=wtime();

        // Compute the histogram
        E_integ=tf.energy_spectrum(es,ed+i*es.nbin);

        // Print diagnostic information
        t2=wtime();
        printf("# Compute frame %d [%d, %.4g s, %.4g s] {%g m^2/s^2}\n",i,l,t1-t0,t2-t1,v2_phys*E_integ);
        t0=t2;
    }

    // Output the histograms
    char buf[128];
    sprintf(buf,"en_spec_%d.dat",nmode);
    FILE *fp=safe_fopen(buf,"w");
    for(int j=0;j<es.nbin;j++) {
        fprintf(fp,"%d",j);
        for(int i=0;i<=nhist;i++) {
            en_spec_data &ee=ed[i*es.nbin+j];
            if(ee.n!=0) fprintf(fp," %.10g %.10g",ee.k,ee.E);
            else fputs(" NaN NaN",fp);
        }
        fputc('\n',fp);
    }

    // Close the file and remove the dynamically allocated memory
    fclose(fp);
    delete [] ed;
}
