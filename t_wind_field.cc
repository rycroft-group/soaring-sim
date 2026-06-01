#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.hh"
#include "tf_grid.hh"

// The output directory to write to
const char fn[]="ts2.odr";

// The number of modes in each direction
const int nx=128,ny=128,nz=128;

// A small function for computing the square of a number
inline double sqr(double x) {
    return x*x;
}

/** Computes the scaling factor to apply to a mode, to take
 * into account symmetries.
 * \param[in] (i,j,k) the mode index to consider.
 * \return The scaling factor, which can be zero, one, or two. */
inline int f_mode(int i,int j,int k) {
    return i!=0?1:(2*k>=nz?(2*k==nz?1:(2*j==ny?1:0))
                          :(k>0?(2*j>=ny?(2*j==ny?1:2):2)
                               :(2*j>=ny?(2*j==ny?1:0):(j==0?0:2))));
}

int main() {

#ifdef FFTW_OMP
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    // Create the output directory
    mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Compute the alpha parameter required to achieve a given RMS wind speed
    double lx=1,ly=1,lz=1,
           s=0,facx=2*M_PI/lx,facy=2*M_PI/ly,facz=2*M_PI/lz,
           w_rms=sqrt(3.0);
#pragma omp parallel for reduction(+:s)
    for(int k=0;k<nz;k++) {
        double w=sqr(facz*(k>nz/2?nz-k:k));
        for(int j=0;j<ny;j++) {
            double ww=w+sqr(facy*(j>ny/2?ny-j:j));
            for(int i=0;i<((nx>>1)+1);i++) if(f_mode(i,j,k)!=0)
                s+=pow(ww+sqr(facx*i),-11/6.);
        }
    }
    double alpha=4*M_PI*w_rms*w_rms/(facx*facy*facz*s);

    // Create the turbulence simulation, and initialize the velocity modes in
    // steady state
    turb_fluid_grid tf(nx,ny,nz,0,lx,0,ly,0,lz,1.0,alpha);
    tf.init_steady_state();

    // Output cross-sections of the velocity
    char buf[256];
    tf.transform();
    for(int i=0;i<nz;i++) {
        sprintf(buf,"%s/uz.%d",fn,i);
        tf.output_z(buf,2,i);
        printf("# Output frame %d\n",i);
    }

    // Optional code to compute the mean squared value of each velocity
    // component. Note that these instantaneous values may not match the target
    // value precisely, because of random temporal fluctuations.
    //tf.mean_squared(mux,muy,muz);
    //printf("%g %g %g\n",mux,muy,muz);
}
