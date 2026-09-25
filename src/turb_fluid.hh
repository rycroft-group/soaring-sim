#ifndef TURB_FLUID_HH
#define TURB_FLUID_HH

#include <cstdio>
#include <cstring>
#include <cmath>

#include <gsl/gsl_rng.h>
#include <fftw3.h>

#include "common.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

/** \brief A class for simulating a 3D turbulence model in terms of Fourier
 * modes. */
class turb_fluid {
    public:
        /** The number of gridpoints in the x direction. */
        const int m;
        /** The number of gridpoints in the y direction. */
        const int n;
        /** The number of gridpoints in the z direction. */
        const int o;
        /** The number of points in the x direction of the Fourier
         * transformed domain. */
        const int fftm;
        /** The total number of gridpoints in the xy plane. */
        const int mn;
        /** The total number of gridpoints. */
        const int mno;
        /** The lower bound of the domain the x direction. */
        const double ax;
        /** The upper bound of the domain the x direction. */
        const double bx;
        /** The lower bound of the domain the y direction. */
        const double ay;
        /** The upper bound of the domain the y direction. */
        const double by;
        /** The lower bound of the domain the z direction. */
        const double az;
        /** The upper bound of the domain the z direction. */
        const double bz;
        /** The grid spacing in the x direction. */
        const double dx;
        /** The grid spacing in the y direction. */
        const double dy;
        /** The grid spacing in the z direction. */
        const double dz;
        /** The trig factor in the x direction. */
        const double facx;
        /** The trig factor in the y direction. */
        const double facy;
        /** The trig factor in the z direction. */
        const double facz;
        /** The reciprocal of the constant so that timescales of modes go like C*k^{-2/3}. */
        const double Cinv;
        /** The constant so that energy scales like E(k)=alpha*k^{-5/3}. */
        const double alpha;
        /** The normalizing constant setting the overall scale of the Fourier
         * modes for velocity. */
        const double fnor;
        /** The Fourier transform of the velocity. */
        fftw_complex* kk;
        turb_fluid(int m_,int n_,int o_,double ax_,double bx_,double ay_,double by_,double az_,double bz_,double Cinv_,double alpha_,unsigned long seed=1);
        ~turb_fluid();
        void init_zero();
        template<int mode>
        void update_random(double dt);
        void mean_revert(double T);
        inline void step_forward(double dt) {update_random<0>(dt);}
        inline void init_steady_state() {update_random<1>(0.);}
        void histogram(double *hi,int nbin,double &hmax);
        double est_max_timestep();
        void vel(double x,double y,double z,double &ux,double &uy,double &uz);
        void vel_multi(int q,double *pos,double *vel);
        void vel_dot_multi(int q,double *pos,double *dir,double *vel,double *dvel);
        void sample_vel_stats(int nsamp,double &ubar,double &vbar,double &wbar,double &urms,double &vrms,double &wrms);
        double mean_rms();
        void correl_init(double *samp,int nsamp);
        void correl_function(double *w,int nbin,double mrad,double *samp,int nsamp);
        void allocate_vel_table(int slots,bool extended=false);
        /** Copies the Fourier mode coefficients from another turbulent fluid
         * class (with an identical size) to this one.
         * \param[in] tf a reference to the other turbulent fluid class. */
        inline void copy_modes(turb_fluid &tf) {
            memcpy(kk,tf.kk,sizeof(fftw_complex)*3*fftm*n*o);
        }
        void save(FILE *fp);
        /** Saves the complete state of the modes and the parameters
         * to a file.
         * \param[in] filename the name of the file to write to. */
        inline void save(const char* filename) {
            FILE *fp=safe_fopen(filename,"wb");
            save(fp);
            fclose(fp);
        }
    protected:
        /** Computes the square of a number.
         * \param[in] x the number.
         * \return The square. */
        inline double sqr(double x) {return x*x;}
        /** Computes the scaling factor to apply to a mode, to take
         * into account symmetries.
         * \param[in] (i,j,k) the mode index to consider.
         * \return The scaling factor, which can be zero, one, or two. */
        inline int f_mode(int i,int j,int k) {
            return i!=0?1:(2*k>=o?(2*k==o?1:(2*j==n?1:0))
                                 :(k>0?(2*j>=n?(2*j==n?1:2):2)
                                      :(2*j>=n?(2*j==n?1:0):(j==0?0:2))));
        }
        void vel_stats_internal(double *up,int nsamp,double &ubar,double &vbar,double &wbar,double &urms,double &vrms,double &wrms);
    private:
        /** Temporary space for computing the rescaled coordinates. */
        double* htab;
        /** Temporary space for assembling the Fourier coefficient table. */
        double* ftab;
        /** The total number of slots available in the Fourier coefficient
         * table. */
        int fslots;
        /** The number of threads. */
        const int nt;
        /** Temporary space for computing the row and slice sums in the
         * multi-velocity routine. */
        double** const rtab;
        /** The array of GSL random number generators. */
        gsl_rng **rng;
        void fourier_table(int q,double *pos);
        /** Computes the modulus squared of a complex (x,y,z) vector.
         * \param[in] kp a pointer to the start of the vector.
         * \return The modulus squared. */
        inline double complex_msq(fftw_complex *kp) {
            return sqr(kp[0][0])+sqr(kp[0][1])
                  +sqr(kp[1][0])+sqr(kp[1][1])
                  +sqr(kp[2][0])+sqr(kp[2][1]);
        }
        inline void xfill(double *rp,fftw_complex *kp,double *xp) {
            rp[0]+=kp[0][0]*xp[0]-kp[0][1]*xp[1];
            rp[1]+=kp[0][0]*xp[1]+kp[0][1]*xp[0];
            rp[2]+=kp[1][0]*xp[0]-kp[1][1]*xp[1];
            rp[3]+=kp[1][0]*xp[1]+kp[1][1]*xp[0];
            rp[4]+=kp[2][0]*xp[0]-kp[2][1]*xp[1];
            rp[5]+=kp[2][0]*xp[1]+kp[2][1]*xp[0];
        }
        inline void yfill(double *sp,double *rp,double *yp) {
            sp[0]+=rp[0]*yp[0]-rp[1]*yp[1];
            sp[1]+=rp[0]*yp[1]+rp[1]*yp[0];
            sp[2]+=rp[2]*yp[0]-rp[3]*yp[1];
            sp[3]+=rp[2]*yp[1]+rp[3]*yp[0];
            sp[4]+=rp[4]*yp[0]-rp[5]*yp[1];
            sp[5]+=rp[4]*yp[1]+rp[5]*yp[0];
        }
        inline void zfill(double *vp,double *sp,double zcos,double zsin) {
            double vx=sp[0]*zcos-sp[1]*zsin,
                   vy=sp[2]*zcos-sp[3]*zsin,
                   vz=sp[4]*zcos-sp[5]*zsin;
#pragma omp atomic
            *vp+=vx;
#pragma omp atomic
            vp[1]+=vy;
#pragma omp atomic
            vp[2]+=vz;
        }
        /** Returns the thread number, if the code was compiled with OpenMP
         * support. Otherwise it returns zero.
         * \return The thread number. */
#ifdef _OPENMP
        inline int thread_num() {return omp_get_thread_num();}
#else
        inline int thread_num() {return 0;}
#endif
};

#endif
