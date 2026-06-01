#ifndef K_FUNC_HH
#define K_FUNC_HH

#define _USE_MATH_DEFINES
#include <cmath>
#include <gsl/gsl_integration.h>

/** The size of the adaptive integration workspace. */
const int kernel_func_ws_size=1024;

/** The relative integration tolerance. */
const double kernel_func_integ_tol=1e-11;

/** Whether to give a warning about bilinear interpolation values
 * being out of range. */
const bool kernel_func_warn_too_large=true;

class kernel_func {
    public:
        /** The number of modes in one direction in the fluid field. */
        const int N;
        /** The number of R bins in the covariance table. */
        const int rs;
        /** The number of out-of-range table calls. */
        int miss;
        /** The maximum R value in the covariance table. */
        const double rmax;
        /** The turbulent fluid box size. */
        const double B;
         /** The minimum resolved wave number in the turbulent fluid. */
        const double kmin;
        /** The maximum resolved wave number in the turbulent fluid. */
        const double kmax;
        /** Normalization constant used in the covariance function,
         * which incorporates the definition of alpha from the paper. */
        const double anor;
        kernel_func(int N_,double B_,int rs_,double rmax_,int blsize);
        virtual ~kernel_func();
        virtual double eval(double r,double t) = 0;
    protected:
        /** A factor used in transforming the r coordinate. */
        const double rfac;
        /** A factor used in inverse transforming the r coordinate. */
        const double irfac;
        /** The covariance function evaluated on a grid for rapid bilinear
         * interpolation. */
        double* const bl;
        inline double rstretch(double r) {return rfac*r*r;}
        inline double irstretch(double r) {return irfac*sqrt(r);}
};

class kernel_rt : public kernel_func {
    public:
        /** The number of T bins in the covariance table. */
        const int ts;
        /** The maximum T value in the covariance table. */
        const double tmax;
        /** The reciprocal of the time scale of the wind field. */
        const double Cinv;
        kernel_rt(int N_,double B_,double Cinv_,int rs_,int ts_,double rmax_,double tmax_);
        virtual double eval(double r,double t);
        double integrate(double r,double t);
        void output_table(const char* filename);
    private:
        /** A factor used in transforming the t coordinate. */
        const double tfac;
         /** A factor used in inverse transforming the t coordinate. */
        const double itfac;
        inline double tstretch(double t) {return tfac*t*t;}
        inline double itstretch(double t) {return itfac*sqrt(t);}
};

class kernel_r : public kernel_func {
    public:
        kernel_r(int N_,double B_,int rs_,double rmax_);
        virtual double eval(double r,double t);
        double integrate(double r);
        void output_table(const char* filename);
};

#endif
