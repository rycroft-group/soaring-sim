#include "gpr.hh"
#include "tf_grid.hh"
#include "file_output.hh"
#include "common.hh"

#include <cstdio>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <sys/stat.h>

/* The number of samples of the velocity field to make in the GPR class, in one
 * direction. */
const int m=24;

/** The number of samples in one direction to make in the output grid. */
const int n=26;

int main(){
     // Create the Gaussian process regression class and initialize the lookup
    // table for kernel
    gpr g(m*m);
    g.init(1000,50,sqrt(3.),5);

    double box_size=1.;
    int nmode = 64;

    // Create the turbulent fluid field and initialize the modes in steady state
    turb_fluid_grid tf(nmode,nmode,nmode,0,box_size,0,box_size,0,box_size,0.67,3.47);
    tf.init_steady_state();

    // Compute the velocity on a grid for initializing the Gaussian process
    double h=box_size/m;
    double x,y,z,wx,wy,wz;
    z=0.5;
    for(int j=0;j<m;j++){
        y=j*(h);
        for(int i=0;i<m;i++){
            x=i*(h);
            tf.vel(x,y,z,wx,wy,wz);
            printf("w %g %g %g\n",wx,wy,wz);
            g.add_measurement(x/box_size,y/box_size,z/box_size,0.,wx,wy,wz);
        }

    }

    h=box_size/n;
    for(int j=0;j<n;j++){
        y=j*(h);
        for(int i=0;i<n;i++){
            x=i*(h);
            g.orig_predict(x/box_size,y/box_size,z/box_size,0.,wx,wy,wz);
            printf("predict %g %g %g\n",wx,wy,wz);
            tf.vel(x,y,z,wx,wy,wz);
            printf("actual %g %g %g\n",wx,wy,wz);
        }

    }

}
