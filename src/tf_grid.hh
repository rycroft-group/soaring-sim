#ifndef TURB_FLUID_GRID_HH
#define TURB_FLUID_GRID_HH

#include <cstring>
#include <cmath>

#include "turb_fluid.hh"
#include <fftw3.h>

/** A default value of the free parameter used in tricubic interpolation. */
const double turb_fluid_grid_a_c_default=-0.84;

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
         /** Performs tricubic interpolation of the velocity.
         * \param[in] (x,y,z) the position at which to interpolate.
         * \param[out] (ux,uy,uz) the velocity vector. */
        inline void cub_interp(double x,double y,double z,double &ux,double &uy,double &uz) {
            four_pt_interp(true,x,y,z,ux,uy,uz);
        }
        /** Performs interpolation of the velocity using the Lanczos2 kernel.
         * \param[in] (x,y,z) the position at which to interpolate.
         * \param[out] (ux,uy,uz) the velocity vector. */
        inline void la2_interp(double x,double y,double z,double &ux,double &uy,double &uz) {
            four_pt_interp(false,x,y,z,ux,uy,uz);
        }
        void full_dft();
        inline void grid_vel_stats(double &ubar,double &vbar,double &wbar,double &urms,double &vrms,double &wrms) {
            vel_stats_internal(uu,mno,ubar,vbar,wbar,urms,vrms,wrms);
        }
        void mean_squared(double &mux,double &muy,double &muz);
        inline void set_cubic_a_param(double a_c_) {a_c=a_c_;}
    protected:
        void interp_slice(double &ux,double &uy,double &uz,double *s,int *q,double *up,double sf);
        void setup_cubic_basis(double x,double y,double z,double *s);
        void setup_lanczos2_basis(double x,double y,double z,double *s);
        inline int step_mod(int a,int b) {return a>=0?a%b:b-1-(b-1-a)%b;}
        inline int step_int(double a) {return a<0?int(a)-1:int(a);}
        /** Takes a given position at which to perform interpolation, and
         * determines which grid box it is within. Updates the position
         * coordinates to cover the ranges [0,1]^3 within the grid box.
         * \param[in,out] (x,y,z) the position to consider, which is remapped
         *                        into the box coordinate system upon
         *                        completion.
         * \param[out] (i,j,k) the grid index of the box to consider, remapped
         *                     into the primary domain if necessary using
         *                     periodicity. */
        inline void grid_remap(double &x,double &y,double &z,int &i,int &j,int &k) {

            // Determine which box of the fluid grid we are in
            i=step_int(x=(x-ax)*xsp);
            j=step_int(y=(y-ay)*ysp);
            k=step_int(z=(z-az)*zsp);

            // Calculate the fractional position within the box
            x-=i;y-=j;z-=k;

            // Ensure that the box indices are in the primary domain,
            // i.e. 0 <= i < m, 0 <= j < n, 0 <= k < o
            i=step_mod(i,m);j=step_mod(j,n);k=step_mod(k,o);
        }
        /** Calculates the basis function coefficients in one coordinate
         * direction for tricubic interpolation.
         * \param[in] x the fractional position of the point.
         * \param[out] s a pointer to store the interpolant coefficients to. */
        inline void cubic_coeffs(double x,double *s) {
            double xx=x*x;
            *s=a_c*x*(1-2*x+xx);
            s[1]=((a_c+2)*x-(a_c+3))*xx+1;
            s[2]=(-(a_c+2)*xx+(2*a_c+3)*x-a_c)*x;
            s[3]=a_c*xx*(1-x);
        }
        /** Evaluates the Lanczos2 kernel near zero via a fourth-order Taylor
         * expansion to avoid floating point inaccuracies. This formula should
         * be precise to within machine precision for |x|<1e-3.
         * \param[in] x the function argument.
         * \return The Taylor expansion. */
        inline double lanczos2_tay_zero(double x) {
            double xx=x*x;
            return 1+xx*(-2.0561675835602821+xx*1.5389283479328381);
        }
        /** Calculates the basis function coefficients in one coordinate
         * direction for Lanczos2 interpolation.
         * \param[in] x the fractional position of the point.
         * \param[out] s a pointer to store the interpolant coefficients to. */
        inline void lanczos2_coeffs(double x,double *s) {
            const double k=4/(M_PI*M_PI),tol=1e-3;
            double t=0.5*M_PI*x,sn=sin(t),cs=cos(t),
                   sc=k*sn*cs,vs2=sc*sn,vc2=sc*cs,
                   xp=x+1,xm=x-1,xmm=x-2;
            *s=-vc2/(xp*xp);
            s[1]=fabs(x)>tol?vs2/(x*x):lanczos2_tay_zero(x);
            s[2]=fabs(xm)>tol?vc2/(xm*xm):lanczos2_tay_zero(xm);
            s[3]=-vs2/(xmm*xmm);
        }
        /** Sets up the memory strides to be used in a four-point interpolation
         * scheme (cubic or Lanczos2) in one direction.
         * \param[in] (i,j,k) the grid box where the interpolation is taking place.
         * \param[out] q an array of length 12 for storing the memory slides (4
         *               for each coordinate direction). */
        inline void setup_memory_strides(int i,int j,int k,int *q) {
            memory_strides(i,m,q,3);
            memory_strides(j,n,q+4,3*m);
            memory_strides(k,o,q+8,3*m*n);
        }
        /** The FFTW plan for converting the velocity from the frequency
         * space to the physical space. */
        fftw_plan plan_bck;
    private:
        /** The free parameter in the tricubic interpolation. */
        double a_c;
        /** Computes the memory strides for four-point interpolation (cubic or
         * Lanczos2) in one direction.
         * \param[in] i the grid index of the block being interpolated.
         * \param[in] d the number of gridpoints in the dimension being
         *              considered.
         * \param[out] q a pointer to the memory stride lengths.
         * \param[in] mul a base memory stride length. */
        inline void memory_strides(int i,int d,int *q,int mul) {
            *q=i>0?-mul:mul*(d-1);
            q[1]=0;
            q[2]=i<d-1?mul:mul*(1-d);
            q[3]=i<d-2?2*mul:mul*(2-d);
        }
        void four_pt_interp(bool cub,double x,double y,double z,double &ux,double &uy,double &uz);
        inline void interp_line(double &vx,double &vy,double &vz,double *s,int *q,double *up,double sf);
        void output_cross_section(const char* filename,double *uu,double a1,double d1,double a2,double d2,int l1,int l2,int s1,int s2);
};

#endif
