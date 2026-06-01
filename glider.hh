#ifndef GLIDER_HH
#define GLIDER_HH

#include <cstdio>
#include <cmath>

#include "turb_fluid.hh"
#include <gsl/gsl_rng.h>

struct glider {
    /** The glider x position. */
    double rx;
    /** The glider y position. */
    double ry;
    /** The glider z position. */
    double rz;
    /** The glider x velocity. */
    double ux;
    /** The glider y velocity. */
    double uy;
    /** The glider z velocity. */
    double uz;
    /** The current banking position. */
    short bank;
    /** Creates the glider class, leaving the class members uninitialized. */
    glider() {}
    /** Creates the glider class, copying the position, velocity, and banking
     * angle from another.
     * \param[in] p the glider glass to copy from. */
    glider(glider &p) : rx(p.rx), ry(p.ry), rz(p.rz), ux(p.ux), uy(p.uy), uz(p.uz), bank(p.bank) {}
    /** Creates the glider class, initializing the members as a forward Euler step
     * from two others.
     * \param[in] p the base class to copy from.
     * \param[in] dt the timestep.
     * \param[in] q the class with the rate of changes. */
    glider(glider &p,double dt,glider &q) : rx(p.rx+dt*q.rx), ry(p.ry+dt*q.ry),
        rz(p.rz+dt*q.rz), ux(p.ux+dt*q.ux), uy(p.uy+dt*q.uy), uz(p.uz+dt*q.uz), bank(p.bank) {}
    /** Initializes the members of the glider class.
     * \param[in] (rx_,ry_,rz_) the glider positon.
     * \param[in] (ux_,uy_,uz_) the glider velocity.
     * \param[in] bank_ the banking angle step. */
    inline void init(double rx_,double ry_,double rz_,double ux_,double uy_,double uz_,short bank_) {
        rx=rx_;
        ry=ry_;
        rz=rz_;
        ux=ux_;
        uy=uy_;
        uz=uz_;
        bank=bank_;
    }
    /** Copies the members of the glider class from another.
     * \param[in] p the glider class to copy from. */
    inline void operator=(glider &p) {
        rx=p.rx;
        ry=p.ry;
        rz=p.rz;
        ux=p.ux;
        uy=p.uy;
        uz=p.uz;
        bank=p.bank;
    }
    /** Updates the members of class according to a forward Euler step.
     * \param[in] dt the timestep.
     * \param[in] q the class with the rate of changes. */
    inline void update(double dt,glider &q) {
        rx+=dt*q.rx;
        ry+=dt*q.ry;
        rz+=dt*q.rz;
        ux+=dt*q.ux;
        uy+=dt*q.uy;
        uz+=dt*q.uz;
    }
    /** Initializes the class members as a forward Euler step from two others.
     * \param[in] p the base class to copy from.
     * \param[in] dt the timestep.
     * \param[in] q the class with the rate of changes. */
    inline void euler(glider &p,double dt,glider &q) {
        rx=p.rx+dt*q.rx;
        ry=p.ry+dt*q.ry;
        rz=p.rz+dt*q.rz;
        ux=p.ux+dt*q.ux;
        uy=p.uy+dt*q.uy;
        uz=p.uz+dt*q.uz;
        bank=p.bank;
    }
    /** Updates the members of the class according to an improved Euler step.
     * \param[in] dt the timestep.
     * \param[in] (q,r) the two classes with the two intermediate steps
     *                  (k1 & k2) in the improved Euler method. */
    inline void improv_e(double dt,glider &q,glider &r) {
        double hdt=0.5*dt;
        rx+=hdt*(q.rx+r.rx);
        ry+=hdt*(q.ry+r.ry);
        rz+=hdt*(q.rz+r.rz);
        ux+=hdt*(q.ux+r.ux);
        uy+=hdt*(q.uy+r.uy);
        uz+=hdt*(q.uz+r.uz);
    }
    /** Copies the glider position to a given memory location.
     * \param[in] p a pointer to the memory location to copy to. */
    inline void copy_pos(double *p) {
        *p=rx;
        p[1]=ry;
        p[2]=rz;
    }
    /** Copies the glider velocity to a given memory location.
     * \param[in] p a pointer to the memory location to copy to. */
    inline void copy_vel(double *p) {
        *p=ux;
        p[1]=uy;
        p[2]=uz;
    }
    inline void output(double time,FILE *fp) {
        fprintf(fp,"%g %g %g %g %g %g %g %d\n",time,rx,ry,rz,ux,uy,uz,bank);
    }
    inline void write_pos(FILE *fp) {
        fprintf(fp," %g %g %g %d",rx,ry,rz,bank);
    }
    inline void write_min_binary(FILE *fp) {

        // Convert the glider position to single precision
        float fl[3];
        *fl=rx;fl[1]=ry;fl[2]=rz;
        fwrite(fl,sizeof(float),3,fp);

        // Write the bank angle
        fwrite(&bank,sizeof(short),1,fp);
    }
    inline void write_binary(FILE *fp,double *w) {

        // Convert the glider position to single precision
        float fl[6];
        *fl=rx;fl[1]=ry;fl[2]=rz;
        fl[3]=*w;fl[4]=w[1];fl[5]=w[2];
        fwrite(fl,sizeof(float),6,fp);

        // Write the bank angle
        fwrite(&bank,sizeof(short),1,fp);
    }
    inline double energy(double *w) {
        double vx=ux-*w,vy=uy-w[1],vz=uz-w[2];
        return 0.5*(vx*vx+vy*vy+vz*vz)+rz;
    }
};

class glider_model {
    public:
        /** The dimensionless lift coefficient. */
        const double c_L;
        /** The dimensionless drag coefficient. */
        const double c_D;
        /** The 2-norm of the (lift,drag) vector. */
        const double s;
        /** The steady-state horizonal velocity for straight flight with zero
         * banking angle. */
        const double uh0;
        /** The steady-state vertical velocity for straight flight with zero
         * banking angle. */
        const double uz0;
        /** The minimum bank range step. */
        const short br_min;
        /** The maximum bank range step. */
        const short br_max;
        /** The minimum change in bank range step. */
        const short bc_min;
        /** The maximum change in bank range step. */
        const short bc_max;
        /** The total number of bank angles. */
        const short rtot;
        /** The total number of possible changes. */
        const short ctot;
        /** The step size in the banking angle. */
        const double bstep;
        glider_model(double c_L_,double c_D_,short br_min_,short br_max_,short bc_min_,short bc_max_,double bstep_);
        /** The class destructor frees the dynamically allocated memory. */
        ~glider_model() {
            delete [] mu_tab;
            delete [] v_act;
        }
        void euler(double dt,glider &p,double *w);
        void euler(int n,double dt,glider *p,double *w);
        void ie_step1(double dt,glider &p,glider &q,glider &r,double *w);
        void ie_step1(int n,double dt,glider *p,glider *g_ie,double *w);
        void ie_step2(double dt,glider &p,glider &q,glider &r,double *w);
        void ie_step2(int n,double dt,glider *p,glider *g_ie,double *w);
        void ff(int n,glider *p,glider *q,double *w);
        void random_move(glider &p,gsl_rng *rng);
        /** Computes the number of valid moves that glider can currently make,
         * and provides a pointer to a list of valid moves.
         * \param[in] p the current glider state.
         * \param[out] v_act_ a pointer to the list of valid moves.
         * \return The total number of valid moves. */
        inline short valid_actions(glider &p,short *&v_act_) {
            short b_lo=p.bank+bc_min,b_hi=p.bank+bc_max;
            if(b_lo<br_min) b_lo=br_min;
            if(b_hi>br_max) b_hi=br_max;
            v_act_=v_act0+b_lo;
            return b_hi-b_lo+1;
        }
        inline double steady_speed() {
            return 1./sqrt(s);
        }
        double d_energy(glider &g,double *w,double *dw);
        double d_energy_component(glider &g,double *w,double *dw,double &alt,double &wacc);
        double max_timestep();
    private:
        /** The table of cosine and sine bank angles. */
        double* const mu_tab;
        /** A pointer to the zero cosine and sine angles. */
        double* const mu0;
        /** A table for communicating the possible valid actions to the MCTS
         * class. */
        short* const v_act;
        /** A pointer to the zero move in the valid action array. */
        short* const v_act0;
        inline void ff(glider &p,glider &q,double *w);
};

double two_norm(glider &p,glider &q);

#endif
