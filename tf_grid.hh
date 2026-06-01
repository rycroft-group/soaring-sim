#ifndef TURB_FLUID_GRID_HH
#define TURB_FLUID_GRID_HH

#include <cstring>
#include <cmath>

#include "turb_fluid.hh"
#include <fftw3.h>

/** \brief An extension of the 3D turbulence class that can also perform
 * Fourier transforms and evaluate the velocity on a spatial grid. */
class turb_fluid_grid : public turb_fluid {
    public:
        double aa;
        /** The inverse grid spacing in the x direction. */
        const double xsp;
        /** The inverse grid spacing in the y direction. */
        const double ysp;
        /** The inverse grid spacing in the z direction. */
        const double zsp;
        /** The size of the Fourier mode array. */
        const size_t ksize;
        /** A copy of the Fourier modes. */
        fftw_complex* const kcopy;
        /** The velocity. */
        double* const uu;
        turb_fluid_grid(int m_,int n_,int o_,double ax_,double bx_,double ay_,double by_,double az_,double bz_,double C_,double alpha_,unsigned long seed=1);
        ~turb_fluid_grid();
        /** Outputs a cross-section fo the velocity for a slice where
         * x=constant.
         * \param[in] filename the name of the file to write to.
         * \param[in] fld the velocity component to write (x=0, y=1, z=2).
         * \param[in] i the x grid index to write. */
        inline void output_x(const char* filename,int fld,int i) {
            output_cross_section(filename,uu+3*i+fld,ay,dy,az,dz,n,o,3*m,3*mn);
        }
        /** Outputs a cross-section fo the velocity for a slice where
         * y=constant.
         * \param[in] filename the name of the file to write to.
         * \param[in] fld the velocity component to write (x=0, y=1, z=2).
         * \param[in] j the y grid index to write. */
        inline void output_y(const char* filename,int fld,int j) {
            output_cross_section(filename,uu+3*j*m+fld,ax,dx,az,dz,m,o,3,3*mn);
        }
        /** Outputs a cross-section fo the velocity for a slice where
         * z=constant.
         * \param[in] filename the name of the file to write to.
         * \param[in] fld the velocity component to write (x=0, y=1, z=2).
         * \param[in] k the z grid index to write. */
        inline void output_z(const char* filename,int fld,int k) {
            output_cross_section(filename,uu+3*k*mn+fld,ax,dx,ay,dy,m,n,3,3*m);
        }
        /** Performs the fast Fourier transform to evaluate the physical
         * velocity on a grid. */
        inline void transform() {
            memcpy(kcopy,kk,ksize);
            fftw_execute(plan_bck);
        }
        void lin_interp(double x,double y,double z,double &ux,double &uy,double &uz);
        void cub_interp(double x,double y,double z,double &ux,double &uy,double &uz);
        void full_dft();
        inline void grid_vel_stats(double &ubar,double &vbar,double &wbar,double &urms,double &vrms,double &wrms) {
            vel_stats_internal(uu,mno,ubar,vbar,wbar,urms,vrms,wrms);
        }
        void mean_squared(double &mux,double &muy,double &muz);
    protected:
        void interp_slice(double &ux,double &uy,double &uz,double *s,int *q,double *up,double sf);
        void bic_basis(double x,int i,int d,double *s,int *q,int mul);
        inline int step_mod(int a,int b) {return a>=0?a%b:b-1-(b-1-a)%b;}
        inline int step_int(double a) {return a<0?int(a)-1:int(a);}
        inline void grid_remap(double &x,double &y,double &z,int &i,int &j,int &k) {

            // Determine which box of the fluid grid we are in
            i=step_int(x=(x-ax)*xsp);
            j=step_int(y=(y-ay)*ysp);
            k=step_int(z=(z-az)*zsp);

            // Calculate the fractional position within the box
            x-=i;y-=j;z-=k;

            // Ensure that the box indices are in the primary grid,
            // i.e. 0 <= i < m, 0 <= j < n, 0 <= k < o
            i=step_mod(i,m);j=step_mod(j,n);k=step_mod(k,o);
        }
        /** The FFTW plan for converting the velocity from the frequency
         * space to the physical space. */
        fftw_plan plan_bck;
    private:
        inline void interp_line(double &vx,double &vy,double &vz,double *s,int *q,double *up,double sf);
        void output_cross_section(const char* filename,double *uu,double a1,double d1,double a2,double d2,int l1,int l2,int s1,int s2);
};

#endif
