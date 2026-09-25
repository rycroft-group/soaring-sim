#include "tf_grid_mr.hh"

#include "common.hh"
#include <cstdio>
#include <cmath>

/** Initializes the three-dimensional turbulent fluid generator with additional
 * FFTW-based routines for evaluating the expected fluid velocity on a grid at
 * some time interval in the future, using Hermite interpolation.
 * \param[in] (m,n,o) the dimensions of the grid.
 * \param[in] (ax_,bx_) the lower and upper x-coordinate simulation bounds.
 * \param[in] (ay_,by_) the lower and upper y-coordinate simulation bounds.
 * \param[in] (az_,bz_) the lower and upper z-coordinate simulation bounds.
 * \param[in] C_ the constant controlling mode timescales.
 * \param[in] alpha_ the constant controlling mode energy scales.
 * \param[in] h_segs_ the number of segments to divide the time interval into.
 * \param[in] Tmr_ the duration of the time interval. */
turb_fluid_grid_mr::turb_fluid_grid_mr(int m_,int n_,int o_,double ax_,double bx_,double ay_,double by_,double az_,double bz_,double Cinv_,double alpha_,int h_segs_,double Tmr_,unsigned long seed)
    : turb_fluid_grid(m_,n_,o_,ax_,bx_,ay_,by_,az_,bz_,Cinv_,alpha_,seed), h_segs(h_segs_),
    Tmr(Tmr_), rho(Tmr/(h_segs*h_segs)), irho(1./rho), umr(new double*[2*h_segs+1]) {

    // Allocate memory for storing the fluid velocity and its derivative at the
    // interpolation control points.
    for(int l=0;l<2*h_segs+1;l++)
        umr[l]=(double*)fftw_malloc(sizeof(double)*3*mno);
}

/** Deallocates the memory used for the storing the interpolation information.
 */
turb_fluid_grid_mr::~turb_fluid_grid_mr() {
    for(int l=2*h_segs;l>=0;l--) fftw_free(umr[l]);
    delete [] umr;
}

/** Calculates the fields for interpolating the mean-reverting velocity field.
 */
void turb_fluid_grid_mr::calc_mr_fields() {
    for(int l=1;l<2*h_segs+2;l++) {
        double ld=double(l/2),T=rho*ld*ld,b1=T*Cinv,b2=-2*ld*rho*Cinv;
#pragma omp for
        for(int k=0;k<o;k++) {
            fftw_complex *kp=kk+3*n*fftm*k,*kcp=kcopy+3*n*fftm*k;
            double w=sqr(facz*(k>o/2?o-k:k));
            for(int j=0;j<n;j++) {
                double ww=w+sqr(facy*(j>n/2?n-j:j));
                for(int i=0;i<fftm;i++,kp+=3,kcp+=3) {

                    // Calculate the scale factor to apply to this mode. Skip
                    // if this mode is not used.
                    int fm=f_mode(i,j,k);
                    if(fm==0) {
                        kcp[0][0]=kcp[0][1]=0.;
                        kcp[1][0]=kcp[1][1]=0.;
                        kcp[2][0]=kcp[2][1]=0.;
                        continue;
                    }

                    // For l even then calculate the velocity at an
                    // interpolation control point. For l odd then calculate
                    // the derivative of velocity.
                    double www=pow(ww+sqr(facx*i),1/3.),
                           c=l&1?b2*www*exp(-b1*www):exp(-b1*www);
                    kcp[0][0]=c*kp[0][0];
                    kcp[0][1]=c*kp[0][1];
                    kcp[1][0]=c*kp[1][0];
                    kcp[1][1]=c*kp[1][1];
                    kcp[2][0]=c*kp[2][0];
                    kcp[2][1]=c*kp[2][1];
                }
            }
        }

        // Perform an FFT to convert the mode information into the velocity
        // information
        fftw_execute_dft_c2r(plan_bck,kcopy,umr[l-1]);
    }
}

/** Evaluates the expected fluid velocity on the grid at a later point in time.
 * The routine used linear interpolation in space and Hermite interpolation in time.
 * \param[in] T the time to consider.
 * \param[in] (x,y,z) the position at which to interpolate.
 * \param[out] (ux,uy,uz) the velocity vector. */
void turb_fluid_grid_mr::lin_interp_mr(double T,double x,double y,double z,double &ux,double &uy,double &uz) {

    // Determine which box of the fluid grid we are in
    int i,j,k;
    grid_remap(x,y,z,i,j,k);

    // Calculate the fractional position within the box
    double gx=1-x,gy=1-y,gz=1-z,s=sqrt(T*irho);;

    // Calculate the displacements to the upper indices
    int xd=3*(i==m-1?1-m:1),
        yd=3*m*(j==n-1?1-n:1),
        zd=3*m*n*(k==o-1?1-o:1),
        dis=3*(i+m*(j+n*k)),

    // Calculate the interpolation interval that the time lies within
        is=int(s);
    if(is<0) fatal_error("T out of range",1);
    else if(is>=h_segs) is=h_segs-1;
    s-=is;

    // Calculate the interpolant, looping over the four Hermite basis
    // functions
    ux=uy=uz=0;
    for(int q=-1;q<3;q++) {

        // Calculate a pointer to the field to use. This must take into account
        // that the velocity field at t=0 is stored separately in the parent
        // turb_fluid_grid class.
        int ind=2*is+q;
        double *up=(ind==-1?uu:umr[2*is+q])+dis,
               *uq=up+zd,fac;

        // Compute the Hermite basis function contribution
        switch(q) {
            case -1: fac=s*s*(2*s-3)+1;break;
            case 0: fac=(s-1)*(s-1)*s;break;
            case 1: fac=s*s*(3-2*s);break;
            case 2: fac=-s*s*(1-s);
        }
        double mz=fac*z,mgz=fac*gz;

        // Perform trilinear interpolation in space
        ux+=mgz*(gy*(gx*up[0]+x*up[xd])+y*(gx*up[yd]+x*up[xd+yd]))
          +mz*(gy*(gx*uq[0]+x*uq[xd])+y*(gx*uq[yd]+x*uq[xd+yd]));
        up++;uq++;
        uy+=mgz*(gy*(gx*up[0]+x*up[xd])+y*(gx*up[yd]+x*up[xd+yd]))
          +mz*(gy*(gx*uq[0]+x*uq[xd])+y*(gx*uq[yd]+x*uq[xd+yd]));
        up++;uq++;
        uz+=mgz*(gy*(gx*up[0]+x*up[xd])+y*(gx*up[yd]+x*up[xd+yd]))
          +mz*(gy*(gx*uq[0]+x*uq[xd])+y*(gx*uq[yd]+x*uq[xd+yd]));
    }
}

/** Evaluates the expected fluid velocity on the grid at a later point in time.
 * The routine used tricubic interpolation in space and Hermite interpolation
 * in time.
 * \param[in] T the time to consider.
 * \param[in] (x,y,z) the position at which to interpolate.
 * \param[out] (ux,uy,uz) the velocity vector. */
void turb_fluid_grid_mr::cub_interp_mr(double T,double x,double y,double z,double &ux,double &uy,double &uz) {

    // Determine which box of the fluid grid that the point is within
    int i,j,k;
    grid_remap(x,y,z,i,j,k);

    // Calculate the interpolation interval that the time lies within
    double w=sqrt(T*irho);;
    int iw=int(w);
    if(iw<0) fatal_error("T out of range",1);
    else if(iw>=h_segs) iw=h_segs-1;
    w-=iw;

    // Compute coefficient tables for the cubic interpolant in the three
    // coordinates. In addition, compute the memory stride lengths, which must
    // account for the periodicity of the domain.
    double s[12];
    int q[12],dis=3*(i+m*(j+n*k));
    bic_basis(x,i,m,s,q,3);
    bic_basis(y,j,n,s+4,q+4,3*m);
    bic_basis(z,k,o,s+8,q+8,3*m*n);

    // Calculate the interpolant, looping over the four Hermite basis
    // functions
    ux=uy=uz=0;
    double vx,vy,vz;
    for(int o=-1;o<3;o++) {
        int ind=2*iw+o;
        double *up=(ind==-1?uu:umr[2*iw+o])+dis,fac;
        vx=vy=vz=0;
        for(int l=0;l<4;l++) interp_slice(vx,vy,vz,s,q,up+q[8+l],s[8+l]);

        // Compute the Hermite basis function contribution
        switch(o) {
            case -1: fac=w*w*(2*w-3)+1;break;
            case 0: fac=(w-1)*(w-1)*w;break;
            case 1: fac=w*w*(3-2*w);break;
            case 2: fac=-w*w*(1-w);
        }

        // Add scaled contribution the velocity interpolation
        ux+=fac*vx;
        uy+=fac*vy;
        uz+=fac*vz;
    }
}
