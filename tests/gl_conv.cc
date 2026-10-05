#include <cstdio>
#include <cstdlib>

#include "glider_test.hh"

int main() {

    // The total number of output points
    const int nframes=1;

    // The duration to integrate over
    const double duration=1;

    // The glider model
    double c_L=1,c_D=1/15.;
    glider_model gm(c_L,c_D,-4,4,-1,1,10);

    // Optional: turbulent wind field can be added
    //turb_fluid tf(64,64,64,0,50,0,50,0,50,1,0.02);
    //tf.init_steady_state();

    // Set up the glider test class. To add the turbulent wind field, change
    // the second argument to &tf.
    glider_test gt(gm,NULL,3,duration,nframes);

    // Simulate a glider with a very small timestep to serve as a reference for
    // a convergence test
    const integration_type itype=it_improv_e;
    glider g_ref,g;
    double dt,dt_pad=0.1;
    gt.integrate(itype,1e-5,dt,g_ref);

    // Simulate the glider with a range of timesteps and calculate the two-norm
    // of the end glider position to the reference position
    while(dt_pad>1e-4) {
        gt.integrate(itype,dt_pad,dt,g);
        printf("%g %g %g\n",dt_pad,dt,two_norm(g,g_ref));
        dt_pad*=0.8;
    }
}
