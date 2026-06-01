#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.hh"
#include "tf_grid.hh"

// The output directory to write to
const char fn[]="ts2.odr";

// The total number of frames
const int nframes=8;

// The duration to integrate over
const double duration=10.;

// The padding factor to apply to the timestep
const double dt_pad=0.2;

int main() {

#ifdef FFTW_OMP
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    // Create the output directory
    mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Create the turbulence simulation, and initialize the velocity modes in
    // steady state
    int nmode=32;
    double time_scale=32.;
    double alpha=M_PI*M_PI/(9.*(pow(sqrt(3.)*M_PI,-2/3.)-pow(sqrt(3.)*M_PI*nmode,-2/3.))),
           C=pow(sqrt(3)*M_PI,2/3.)*time_scale;
    turb_fluid_grid tf(nmode,nmode,nmode,0,1,0,1,0,1,C,alpha);
    tf.init_steady_state();

    // Output the initial cross-section
    char buf[256];
    tf.transform();
    sprintf(buf,"%s/uz.0",fn);
    tf.output_z(buf,0,0);
    puts("# Output frame 0");

    // Compute a timestep
    double dt=dt_pad*tf.est_max_timestep(),sint=duration/nframes;
    int l=static_cast<int>(sint/dt)+1;
    dt=sint/l;

    // Integrate the modes and output additional snapshots
    /*for(int k=1;k<=nframes;k++) {

        // Perform stochastic integration timesteps
        for(int i=0;i<l*100;i++) tf.step_forward(dt);
        t1=wtime();

        // Output a cross-section
        tf.transform();
        sprintf(buf,"%s/uz.%d",fn,k);
        tf.output_z(buf,2,1);

        // Print diagnostic information
        t2=wtime();
        printf("# Output frame %d [%d, %.4g s, %.4g s]\n",k,l,t1-t0,t2-t1);
        t0=t2;
    }*/

    tf.transform();
    for(int i=0;i<1001;i++) {
        double rx=0.001*i,ry=0.0002*i,rz=0.00067*i;
        double wx,wy,wz,vx,vy,vz;
        tf.vel(rx,ry,rz,vx,vy,vz);
        tf.cub_interp(rx,ry,rz,wx,wy,wz);
        printf("%g %g %g %g %g %g %g %g %g\n",rx,ry,rz,vx,vy,vz,wx,wy,wz);
    }
    return 0;

    const int N=131072;
    double *rr=new double[3*N],*vv=new double[3*N];
    for(int i=0;i<3*N;i++) rr[i]=static_cast<double>(rand())*(1./RAND_MAX);

    tf.allocate_vel_table(N);
    tf.vel_multi(N,rr,vv);

    double l2=0;
    for(int i=0;i<3*N;i++) l2+=vv[i]*vv[i];
    l2=1./sqrt(l2);

    int zz=0;
    double t0=wtime(),t1,t2;
    for(double aa=-1.2;aa<-0.3;aa+=0.01) {
        tf.aa=aa;
        double ss=0,wx,wy,wz;
        for(int i=0;i<N;i++) {
            tf.lin_interp(rr[3*i],rr[3*i+1],rr[3*i+2],wx,wy,wz);
            wx-=vv[3*i];
            wy-=vv[3*i+1];
            wz-=vv[3*i+2];
            ss+=wx*wx+wy*wy+wz*wz;
        }
        zz++;
        //printf("%g %g\n",aa,sqrt(ss)*l2);
    }
    printf("# Time: %g ns\n",1e9*(wtime()-t0)/(double(zz)*double(N)));

    delete [] rr;
    delete [] vv;
}
