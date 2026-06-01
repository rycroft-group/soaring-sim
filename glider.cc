#include "glider.hh"

#include <ctime>
#include <gsl/gsl_randist.h>

/** Sets up the class containing all of the parameters for
 * efficiently integrating a glider.
 * \param[in] c_L_ the dimensionless lift coefficient.
 * \param[in] c_D_ the dimensionless drag coefficient.
 * \param[in] (br_min_,br_max_) the range of banking angle steps.
 * \param[in] (bc_min_,bc_max_) the range of banking angle changes.
 * \param[in] bstep_ the banking angle step in degrees. */
glider_model::glider_model(double c_L_,double c_D_,short br_min_,short br_max_,short bc_min_,short bc_max_,double bstep_)
    : c_L(c_L_), c_D(c_D_), s(sqrt(c_D*c_D+c_L*c_L)), uh0(c_L/(s*sqrt(s))),
    uz0(-c_D/(s*sqrt(s))), br_min(br_min_), br_max(br_max_), bc_min(bc_min_),
    bc_max(bc_max_), rtot(br_max-br_min+1), ctot(bc_max-bc_min+1),
    bstep(bstep_), mu_tab(new double[2*rtot]), mu0(mu_tab-2*br_min),
    v_act(new short[rtot]), v_act0(v_act-br_min) {

    // Set up the cosine and sine values of all possible banking angles, as
    // well as the table for communicating valid moves to the MCTS class.
    double *pp=mu_tab,fac=(M_PI/180.)*bstep;short *vp=v_act;
    for(short k=br_min;k<=br_max;k++) {
        *(pp++)=cos(fac*k);
        *(pp++)=sin(fac*k);
        *(vp++)=k;
    }
}

/** Updates the glider according to an Euler step.
 * \param[in] dt the timestep to use.
 * \param[in,out] p the glider position and velocity to update.
 * \param[in] w the three-component wind vector at the glider's position. */
void glider_model::euler(double dt,glider &p,double *w) {
    glider q;
    ff(p,q,w);
    p.update(dt,q);
}

/** Updates a array of gliders according to an Euler step.
 * \param[in] n the length of the array.
 * \param[in] dt the timestep to use.
 * \param[in] p the array of glider positions and velocities.
 * \param[in] w the array of wind vectors at the gliders' positions. */
void glider_model::euler(int n,double dt,glider *p,double *w) {
    glider q;
    for(int i=0;i<n;i++) {
        ff(p[i],q,w+3*i);
        p[i].update(dt,q);
    }
}

/** Updates the glider using the improved Euler method.
 * \param[in] dt the timestep to use.
 * \param[in,out] p the glider position and velocity to update.
 * \param[in] q the results of the first RK step.
 * \param[in,out] r the glider position after applying q. */
void glider_model::ie_step1(double dt,glider &p,glider &q,glider &r,double *w) {
    ff(p,q,w);
    r.euler(p,dt,q);
}

/** Updates an array of gliders using the improved Euler method.
 * \param[in] n the length of the array.
 * \param[in] dt the timestep to use.
 * \param[i]] p the array of glider positions and velocities to update.
 * \param[in] q an array of the results of the first RK step.
 * \param[in] r an array of the glider states after applying q.
 * \param[in] w the wind vector at the glider's position. */
void glider_model::ie_step1(int n,double dt,glider *p,glider *g_ie,double *w) {
    for(int i=0;i<n;i++) {
        ff(p[i],g_ie[2*i+1],w+3*i);
        g_ie[2*i].euler(p[i],dt,g_ie[2*i+1]);
    }
}

/** Updates the glider using the improved Euler method.
 * \param[in] dt the timestep to use.
 * \param[in,out] p the glider position and velocity to update.
 * \param[in] q the results of the first RK step.
 * \param[in,out] r the glider position after applying q.
 * \param[in] w the wind at position given in r. */
void glider_model::ie_step2(double dt,glider &p,glider &q,glider &r,double *w) {
    ff(r,r,w);
    p.improv_e(dt,q,r);
}

/** Updates an array of gliders using the improved Euler method.
 * \param[in] n the length of the array.
 * \param[in] dt the timestep to use.
 * \param[in] p the array of glider positions and velocities to update.
 * \param[in] g_ie information about the first step in the method.
 * \param[in] w an array of wind vectors at the first step locations. */
void glider_model::ie_step2(int n,double dt,glider *p,glider *g_ie,double *w) {
    for(int i=0;i<n;i++) {
        ff(g_ie[2*i],g_ie[2*i],w+3*i);
        p[i].improv_e(dt,g_ie[2*i],g_ie[2*i+1]);
    }
}

/** Computes the change in glider state for an array of gliders.
 * \param[in] n the length of the array.
 * \param[in] p the array of glider positions and velocities.
 * \param[in] q the array of changes.
 * \param[in] w the array of wind vectors at the gliders' positions. */
void glider_model::ff(int n,glider *p,glider *q,double *w) {
    for(int i=0;i<n;i++) ff(p[i],q[i],w+3*i);
}

/** Computes the change in glider state for a given glider.
 * \param[in] p the current glider state.
 * \param[in] q the change in the glider state.
 * \param[in] w the wind vector at the glider's position. */
inline void glider_model::ff(glider &p,glider &q,double *w) {

    // Compute glider's velocity relative to the wind
    double vx=p.ux-*w,vy=p.uy-w[1],vz=p.uz-w[2];

    // Compute the change in the glider's position
    q.rx=p.ux;
    q.ry=p.uy;
    q.rz=p.uz;

    // Compute the change in the glider's velocity, looking up the
    // trigonometric factors of bank angle in the table
    double vh2=vx*vx+vy*vy,v=sqrt(vh2+vz*vz),vh=sqrt(vh2),
           cos_mu=mu0[2*p.bank],vzcm=vz*cos_mu,
           vsm=v*mu0[2*p.bank+1],fac=c_L/vh;
    q.ux=-v*(fac*(vx*vzcm+vy*vsm)+c_D*vx);
    q.uy=-v*(fac*(vy*vzcm-vx*vsm)+c_D*vy);
    q.uz=v*(c_L*vh*cos_mu-c_D*vz)-1;
}

/** Calculates the maximum allowable timestep for integration, in order to
 * ensure stability of the ODE system around a steady solution
 * u=(c_L,0,-c_D)/s^1.5. */
double glider_model::max_timestep() {
    double q=c_D/s;
    return q<sqrt(8)/3.?1.5*q/sqrt(s):4/3./(sqrt(s)*(q+sqrt(q*q-8/9.)));
}

/** Applies a random move to the glider bank angle.
 * \param[in] rng a pointer to the GSL random number generator to use. */
void glider_model::random_move(glider &p,gsl_rng* rng) {
    int e=p.bank+bc_min+gsl_rng_uniform_int(rng,ctot);
    if(e<br_min) e=br_min;
    else if(e>br_max) e=br_max;
    p.bank=e;
}

/** Computes the Euclidean norm of the difference between two glider states.
 * \param[in] (p,q) the two glider states.
 * \return The norm. */
double two_norm(glider &p,glider &q) {
    double drx=p.rx-q.rx,dry=p.ry-q.ry,drz=p.rz-q.rz,
           dux=p.ux-q.ux,duy=p.uy-q.uy,duz=p.uz-q.uz;
    return sqrt(drx*drx+dry*dry+drz*drz
               +dux*dux+duy*duy+duz*duz);
}

/** Calculates the expected change in energy of the glider.
 * \param[in] g a reference to the glider to consider.
 * \param[in] w a pointer to the wind vector for the glider.
 * \param[in] dw a pointer to the change in wind vector along the glider's
 *               trajectory.
 * \return The change in energy. */
double glider_model::d_energy(glider &g,double *w,double *dw) {
    double vx=g.ux-*w,vy=g.uy-w[1],vz=g.uz-w[2],
           vsq=vx*vx+vy*vy+vz*vz;
    return -c_D*vsq*sqrt(vsq)+w[2]-vx*dw[0]-vy*dw[1]-vz*dw[2];
}

/** Calculates the three components of the expected change in energy of the
 * glider.
 * \param[in] g a reference to the glider to consider.
 * \param[in] w a pointer to the wind vector for the glider.
 * \param[in] dw a pointer to the change in wind vector along the glider's
 *               trajectory.
 * \param[out] alt the component from altitude.
 * \param[out] wacc the component due to wind acceleration.
 * \return The component from drag. */
double glider_model::d_energy_component(glider &g,double *w,double *dw,double &alt,double &wacc) {
    double vx=g.ux-*w,vy=g.uy-w[1],vz=g.uz-w[2],
           vsq=vx*vx+vy*vy+vz*vz;
    alt=w[2];
    wacc=-vx*dw[0]-vy*dw[1]-vz*dw[2];
    return -c_D*vsq*sqrt(vsq);
}
