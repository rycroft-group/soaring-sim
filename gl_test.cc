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

    // Set up the glider test class
    glider_test gt(gm,NULL,4,duration,nframes);

    glider g;
    double dt;
    gt.integrate(it_improv_e,0.2,dt,g,"gt.dat");
}
