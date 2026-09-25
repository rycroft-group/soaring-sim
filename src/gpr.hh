#ifndef GPR_HH
#define GPR_HH

#include <cstdio>
#include <cmath>
#include <limits>

#include "k_func.hh"

/** Special value of time to signify an invalid measurement. */
const double gpr_invalid_t=-std::numeric_limits<double>::max();

struct measure_info {
    /** The position of the measurement. */
    double x,y,z;
    /** The time of the measurement. */
    double t;
    inline bool close(double tol,double x_,double y_,double z_) {
        return fabs(x-x_)<tol&&fabs(y-y_)<tol&&fabs(z-z_)<tol;
    }
    inline void set(double x_,double y_,double z_,double t_) {
        x=x_;y=y_;z=z_;t=t_;
    }
    inline double dis_sq(double x_,double y_,double z_) {
        double dx=x-x_,dy=y-y_,dz=z-z_;
        return dx*dx+dy*dy+dz*dz;
    }
    inline void nullify() {
        t=gpr_invalid_t;
    }
    inline bool invalid() {
        return t==gpr_invalid_t;
    }
};

class gpr {
    public:
        /** The maximum number of measurements to use in the regression. */
        int m;
        /** The memory allocation for the maximum size of the covariance
         * matrix. */
        const int mm;
        /** The counter to the last measurement in the array. */
        int v;
        /** The size of the workspace for passing to the LAPACK routines. */
        int lwork;
        /** Whether to always bypass the Woodbury formula and do a full LAPACK
         * computation each time during updates. */
        const bool full_compute;
        /** Whether the array is full or not. */
        bool full;
        /** The measurement information array. */
        measure_info* const M;
        gpr(int m_,kernel_func &KF_,double ker_tol=0,double ker_step=0,bool full_compute_=false);
        ~gpr();
        inline void reset() {
            v=0;full=false;
        }
        void add_measurement(double x,double y,double z,double t,double wx,double wy,double wz);
        inline void add_measurement(double *pos,double t,double *vel) {
            add_measurement(*pos,pos[1],pos[2],t,*vel,vel[1],vel[2]);
        }
        void compute_inverse(double *B);
        void compute_X();
        /** Performs the kernel calculations required to perform predictions.
         * This is only needed when add_measurement has been used. */
        inline void calculate() {
            compute_inverse(B);
            compute_X();
        }
        void update_measurement(double x,double y,double z,double t,double wx,double wy,double wz);
        inline void update_measurement(double *pos,double t,double *vel) {
            update_measurement(*pos,pos[1],pos[2],t,*vel,vel[1],vel[2]);
        }
        void predict(double x,double y,double z,double t,double &wx,double &wy,double &wz);
        void print(double *Q);
        double checksum();
        void diagnostics(FILE *fp);
        void print_times() {
            for(int i=0;i<(full?m:v);i++) printf("%d %g %g %g %g\n",i,M[i].x,M[i].y,M[i].z,M[i].t);
        }
   private:
        inline double min_dis(double x,double y,double z) {
            double minq=1e10,q;
            int n=full?m:v;
            for(int i=0;i<n;i++) {
                q=M[i].dis_sq(x,y,z);
                if(q<minq) minq=q;
            }
            return sqrt(minq);
        }
        /** Computes the kernel function between two measured positions.
         * \param[in] m1 a reference to the first measurement.
         * \param[in] (x,y,z) the position of the second measurement.
         * \param[in] t the time of the second measurement. */
        inline double K(measure_info &m1,double x,double y,double z,double t) {
            if(m1.invalid()) return 0;
            double dx=m1.x-x,dy=m1.y-y,dz=m1.z-z;
            return KF.eval(sqrt(dx*dx+dy*dy+dz*dz),fabs(m1.t-t));
        }
        /** Computes the kernel function between two measured positions.
         * \param[in] m1 a reference to the first measurement.
         * \param[in] m2 a reference to the second measurement. */
        inline double K(measure_info &m1,measure_info &m2) {
            if(m2.invalid()) return 0;
            return K(m1,m2.x,m2.y,m2.z,m2.t);
        }
        /** The covariance matrix. */
        double* A;
        /** The inverse of the covariance matrix. */
        double* B;
        /** The wind vector. */
        double* W;
        /** The vector used to update A^-1.*/
        double* X;
        /** The pivoting array. */
        int* ipiv;
        /** The workspace for passing to the LAPACK routines. */
        double *work;
        /** A reference to the kernel function. */
        kernel_func &KF;
        /** The tolerance used for detecting ill-conditioned kernel covariance
         * matrices, during a measurement substitution. */
        const double ker_det_tol;
        /** The tolerance used for detecting ill-conditioned kernel covariance
         * matrices, during a measurement addition. */
        const double ker_det_tol2;
        /** Counters for the different cases in covariance matrix update. */
        long co[4];
};

#endif
