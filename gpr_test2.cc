#include "k_func.hh"
#include "gpr.hh"
#include "tf_grid.hh"
#include "file_output.hh"
#include "common.hh"

#include <cstdio>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <sys/stat.h>

/** This test verify the gpr prediction accuracy after 0.5s given a coarse sampling/ lines of samples. */

/* The number of samples of the velocity field to make in the GPR class, in one
 * direction. */
const int m=24,mm=m*m;

/** The number of samples in one direction to make in the output grid. */
const int n=128;

/** The memory allocation for assembling the array to pass to the vel_multi routine. */
const int mem=3*(m>n?m*m:n*n);
// const int mem=50;

const double box_size=50.;
const double dt=0.01;

/** The number of modes in a grid is nmode^3. */
const int nmode=64;
/** The time scale of dynamic wind field. */
const double C=32.;

/** Computes the velocity field on a grid.
 * \param[in] pos a pointer to the memory to use for the positions.
 * \param[in] vel a pointer to the memory to use for the velocities.
 * \param[in] tf a reference to the turbulent fluid class.
 * \param[in] g the number of grid points in one direction. */
void compute_grid_velocity(double *pos,double *vel,turb_fluid_grid &tf,int d) {
    double h=box_size/d,*pp=pos,z;

    // Create a grid of positions in a 2D slice
    for(int j=0;j<d;j++) {
        z=h*(j);
	    for(int i=0;i<d;i++) {
            *(pp++)=h*(i);
            *(pp++)=0.5*box_size;
            *(pp++)=z;
        }
    }

    // Compute the velocities at the grid of positions
    tf.allocate_vel_table(d*d);
    tf.vel_multi(d*d,pos,vel);
}

int main() {

    // Create the Gaussian process regression class and initialize the lookup
    // table for kernel
    // gpr g(m*m);
    kernel_rt KF(nmode,box_size,C,100,100,sqrt(3.)*box_size,50);

    gpr g(mm,KF),g2(mm,KF),g3(mm,KF);

    // Create the turbulent fluid field and initialize the modes in steady state
    double alpha=1;
    turb_fluid_grid tf(nmode,nmode,nmode,0,box_size,0,box_size,0,box_size,C,alpha);
    tf.init_steady_state();

    // Compute the velocity on a grid for initializing the Gaussian process
    // regression class
    double pos[mem],vel[mem];
    double wx,wy,wz,t0=0.;
    compute_grid_velocity(pos,vel,tf,m);
    double fsamp[m*m],*fsp=fsamp,fsamp2[m*m],*fsp2=fsamp2;
    for(double *pp=pos,*vp=vel;pp<pos+3*m*m;pp+=3,vp+=3){
        //for(int i=0;i<1;i++) tf.step_forward(0.01);
        //t0+=0.01;
        //tf.vel(*pp,pp[1],pp[2],wx,wy,wz);
        //g.add_measurement(*pp,pp[1],pp[2],t0,wx,wy,wz);
        *(fsp++)=vp[2];
       // printf("wx %g\n",wx);
        g.add_measurement(*pp,pp[1],pp[2],t0,*vp,vp[1],vp[2]);

    }
    for(double *pp=pos,*vp=vel,i=0;pp<pos+3*m*m;pp+=3,vp+=3,i++){
        // tf.vel(*pp,pp[1],pp[2],wx,wy,wz);
        // g.add_measurement(*pp,pp[1],pp[2],0.,wx,wy,wz);
        if(int(i/m)==4 or int(i/m)==8 or int(i/m)==12 or int(i/m)==16){
            *(fsp2++)=vp[2];
            g2.add_measurement(*pp,pp[1],pp[2],0,*vp,vp[1],vp[2]);

        }
        else{
            *(fsp2++)=0.;
        }
    }

    double l=0.5*box_size/m,u=box_size-l;
    gnuplot_output("sample_z.dat",fsamp,m,m,l,u,l,u);
    gnuplot_output("sample_z_line.dat",fsamp2,m,m,l,u,l,u);
    // Compute the velocity on the output grid
    double fwind[n*n],fpred[n*n],fpred2[n*n],*fwp=fwind,*fpp=fpred,*fpp2=fpred2;
    for(int i=0;i<50;i++){
	    t0+=dt;
	    tf.step_forward(dt);
    }
    compute_grid_velocity(pos,vel,tf,n);

    g.calculate();
    g2.calculate();

    // Assemble the output grids with the actual x velocity and the prediction
    // from Gaussian process regression
    for(double *pp=pos,*vp=vel;pp<pos+3*n*n;pp+=3,vp+=3) {
        *(fwp++)=vp[2];
        // printf("position %g %g %g\n",*pp,pp[1],pp[2]);
	   // tf.vel(*pp+0.01,pp[1],pp[2]+0.01,wx,wy,wz);
	   // g.add_measurement(*pp,pp[1],pp[2],0,wx,wy,wz);
        g.predict(*pp,pp[1],pp[2],t0,wx,wy,wz);
        *(fpp++)=wz;
        g2.predict(*pp,pp[1],pp[2],t0,wx,wy,wz);
        *(fpp2++)=wz;
        //printf("acutal %g predict %g\n",*vp,wx);
    }

    // Output the fields in the Gnuplot matrix binary format
    l=0.5*box_size/n,u=box_size-l;
    gnuplot_output("wind_z.dat",fwind,n,n,l,u,l,u);
    gnuplot_output("pred_z.dat",fpred,n,n,l,u,l,u);
    gnuplot_output("pred_z_line.dat",fpred2,n,n,l,u,l,u);

}
