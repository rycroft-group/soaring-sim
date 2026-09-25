#include <cstdio>
#include <cstdlib>

#include "glider_test.hh"

int main() {

    // The total number of output points
    const int nframes=2000;

    // The duration to integrate over
    const double duration=50;

    // The glider model
    double c_L=1,c_D=1/15.;
    glider_model gm(c_L,c_D,-4,4,-1,1,10);

    turb_fluid tf(64,64,64,0,50,0,50,0,50,1,0.02);
    tf.init_steady_state();

    // Set up the glider test class
    glider_test gt(gm,&tf,3,duration,nframes);

    const integration_type itype=it_euler;
    glider g_ref,g;
    double dt,dt_pad=0.01;
    //gt.integrate(itype,1e-3,dt,g_ref,"gd.dat");
    gt.integrate(itype,2e-3,dt,g_ref,"gd2.dat");
    //gt.integrate(itype,4e-3,dt,g_ref,"gd4.dat");

    return 0;
    while(dt_pad>1e-3) {

        gt.integrate(itype,dt_pad,dt,g);
        printf("%g %g\n",dt,two_norm(g,g_ref));

        dt_pad*=0.8;
    }
}
