#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.hh"
#include "tf_grid.hh"

// The side length of the cube for the turbulent wind field
const double L=50;

// Parameter controlling the wind field mode amplitude, based on 1 m/s per
// component for N=64.
const double alpha=0.00913104;

int main(int argc,char **argv) {

    // Print syntax message
    if(argc!=2&&argc!=9) {
        fputs("Syntax: ./t_wind_interp <modes>\n"
              "        ./t_wind_interp <modes> <ax> <ay> <az> <bx> <by> <bz> <steps>\n"
              "\nCompares the exact expansion of the wind field with <modes> Fourier\n"
              "modes to trilinear, Lanczos-2, and tricubic interpolants, by evaluating\n"
              "the different expressions at <steps> samples along a line from\n"
              "<ax>,<ay>,<az> to <bx>,<by>,<bz>. Output is stored to twi.dat.\n\n"
              "If only the <modes> argument is supplied, then a line of the form\n"
              "(lambda,0,0) is used for varying lambda, and the wind at the grid points\n"
              "is stored to twi_grid.dat.\n",stderr);
        return 1;
    }

    // Set up multi-threaded FFTW computations, if available
#ifdef FFTW_OMP
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    // Read the number of Fourier modes to use and check that the value is
    // within a reasonable range
    int st,nmode=atoi(argv[1]);
    if(nmode<=0||nmode>4096) {
        fputs("Mode number out of range\n",stderr);
        return 1;
    }
 
    // Read in additional command-line arguments if available; otherwise use
    // the default set
    double ax,ay,az,bx,by,bz,dx,dy,dz,fac;
    if(argc==2) {
        ax=0;ay=0;az=0;bx=L;by=0;bz=0;st=nmode*16;
    } else {

        // Read in the start and end points of the sample line
        ax=atof(argv[2]);
        ay=atof(argv[3]);
        az=atof(argv[4]);
        bx=atof(argv[5]);
        by=atof(argv[6]);
        bz=atof(argv[7]);
        st=atoi(argv[8]);

        // Read in the number of steps and check that it is a reasonable value
        if(st<=0||st>65536) {
            fputs("Number of steps out of range\n",stderr);
            return 1;
        }
    }
    printf("%d %g %g %g %g %g %g\n",st,ax,ay,az,bx,by,bz);
    
    // Calculate the delta vector between each sample point
    fac=1./st;
    dx=(bx-ax)*fac;dy=(by-ay)*fac;dz=(bz-az)*fac;

    // Create the turbulence simulation, and initialize the velocity modes in
    // steady state. The Cinv parameter is set to zero since this test only
    // requires a frozen wind field.
    turb_fluid_grid tf(nmode,nmode,nmode,0,L,0,L,0,L,0,alpha);
    tf.init_steady_state();

    // Compute the exact expansion of the Fourier modes at the sample points
    double *pos=new double[3*(st+1)],*vel=new double[3*(st+1)],*pp=pos,*vp=vel;
    for(int i=0;i<=st;i++) {
        *(pp++)=ax+i*dx;
        *(pp++)=ay+i*dy;
        *(pp++)=az+i*dz;
    }
    tf.allocate_vel_table(st+1);
    tf.vel_multi(st+1,pos,vel);

    // Output the exact Fourier expansion of the wind along the sample line,
    // along with the trilinear, Lanczos-2, and tricubic interpolants
    FILE *fp=safe_fopen("twi.dat","w");
    tf.transform();
    for(pp=pos;pp<pos+3*(st+1);pp+=3,vp+=3) {
        double u1,v1,w1,u2,v2,w2,u3,v3,w3;
        tf.lin_interp(*pp,pp[1],pp[2],u1,v1,w1);
        tf.la2_interp(*pp,pp[1],pp[2],u2,v2,w2);
        tf.cub_interp(*pp,pp[1],pp[2],u3,v3,w3);
        fprintf(fp,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                *pp,pp[1],pp[2],*vp,vp[1],vp[2],u1,v1,w1,u2,v2,w2,u3,v3,w3);
    }
    fclose(fp);

    // If the default sample line is used, then store the wind velocities at
    // the aligned grid points to another file
    if(argc==2) {
        double h=L/nmode;
        fp=safe_fopen("twi_grid.dat","w");
        vp=tf.uu;
        for(int i=0;i<nmode;i++,vp+=3)
            fprintf(fp,"%g 0 0 %g %g %g\n",i*h,*vp,vp[1],vp[2]);
        vp=tf.uu;
        fprintf(fp,"%g 0 0 %g %g %g",L,*vp,vp[1],vp[2]);
        fclose(fp);
    }

    // Free the dynamically allocated memory
    delete [] vel;
    delete [] pos;
}
