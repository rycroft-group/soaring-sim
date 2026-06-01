#include "turb_fluid.hh"
#include "common.hh"

#include <cmath>
#include <cstdlib>

#include <gsl/gsl_randist.h>

#ifdef _OPENMP
#include "omp.h"
#endif

/** Initializes the three-dimensional turbulent fluid generator.
 * \param[in] (m,n,o) the dimensions of the grid.
 * \param[in] (ax_,bx_) the lower and upper x-coordinate simulation bounds.
 * \param[in] (ay_,by_) the lower and upper y-coordinate simulation bounds.
 * \param[in] (az_,bz_) the lower and upper z-coordinate simulation bounds.
 * \param[in] Cinv_ the reciprocal of the constant controlling mode timescales.
 * \param[in] alpha_ the constant controlling mode energy scales. */
turb_fluid::turb_fluid(int m_,int n_,int o_,double ax_,double bx_,double ay_,double by_,double az_,double bz_,double Cinv_,double alpha_,unsigned long seed)
    : m(m_), n(n_), o(o_), fftm((m>>1)+1), mn(m*n), mno(mn*o), ax(ax_),
    bx(bx_), ay(ay_), by(by_), az(az_), bz(bz_), dx((bx-ax)/m), dy((by-ay)/n),
    dz((bz-az)/o), facx(2*M_PI/(bx-ax)),
    facy(2*M_PI/(by-ay)), facz(2*M_PI/(bz-az)), Cinv(Cinv_), alpha(alpha_),
    fnor(sqrt(1./(48*M_PI)*alpha*facx*facy*facz)),
    kk((fftw_complex*)fftw_malloc(sizeof(fftw_complex)*3*fftm*n*o)),
    htab(new double[2*(fftm+n)+1]), ftab(htab+1), fslots(1),
#ifdef _OPENMP
    nt(omp_get_max_threads()),
#else
    nt(1),
#endif
    rtab(new double*[nt]), rng(new gsl_rng*[nt]) {
#pragma omp parallel
    {
        int t=thread_num();
        rtab[t]=new double[24];
        rng[t]=gsl_rng_alloc(gsl_rng_taus2);
        gsl_rng_set(rng[t],t+seed);
    }
}

/** The class destructor frees the GSL random number generators and the
 * dynamically allocated memory. */
turb_fluid::~turb_fluid() {
    for(int i=0;i<nt;i++) {
        gsl_rng_free(rng[i]);
        delete [] rtab[i];
    }
    delete [] rng;
    delete [] rtab;
    delete [] htab;
    fftw_free(kk);
}

/** Initializes all of the Fourier modes to be zero. */
void turb_fluid::init_zero() {
    for(int i=0;i<3*fftm*n*o;i++) kk[i][0]=kk[i][1]=0;
}

/** Evaluates the velocity at a given position.
 * \param[in] (x,y,z) the position to consider.
 * \param[out] (ux,uy,uz) the velocity. */
void turb_fluid::vel(double x,double y,double z,double &ux,double &uy,double &uz) {
    double hx=(x-ax)*facx,hy=(y-ay)*facy,hz=(z-az)*facz,*ytab=ftab+2*fftm,*fp=ftab;

    // Assemble the table of Fourier coefficients in the x direction,
    // taking into account the symmetry
    *(fp++)=1.;*(fp++)=0.;
    for(int i=1;i<fftm-1;i++) {
        *(fp++)=2*cos(i*hx);
        *(fp++)=2*sin(i*hx);
    }
    if(m&1) {
        *(fp++)=2*cos((fftm-1)*hx);
        *(fp++)=2*sin((fftm-1)*hx);
    } else {
        *(fp++)=cos((fftm-1)*hx);
        *(fp++)=sin((fftm-1)*hx);
    }

    // Assemble the table of Fourier coefficients in the y direction
    *(fp++)=1.;*(fp++)=0.;
    for(int j=1;j<n;j++) {
        int uj=j>n/2?j-n:j;
        *(fp++)=cos(uj*hy);
        *(fp++)=sin(uj*hy);
    }

    // Loop over all Fourier modes and calculate their contribution to the
    // velocity
    ux=uy=uz=0;
#pragma omp parallel for reduction(+:ux) reduction(+:uy) reduction(+:uz)
    for(int k=0;k<o;k++) {
        int uk=k>o/2?k-o:k;
        fftw_complex *kp=kk+3*n*fftm*k;
        double sx=0,jx=0,sy=0,jy=0,sz=0,jz=0;
        for(double *yp=ytab;yp<ytab+2*n;yp+=2) {
            double rx=0,ix=0,ry=0,iy=0,rz=0,iz=0;
            for(double *xp=ftab;xp<ftab+2*fftm;xp+=2,kp+=3) {
                rx+=kp[0][0]*xp[0]-kp[0][1]*xp[1];
                ix+=kp[0][0]*xp[1]+kp[0][1]*xp[0];
                ry+=kp[1][0]*xp[0]-kp[1][1]*xp[1];
                iy+=kp[1][0]*xp[1]+kp[1][1]*xp[0];
                rz+=kp[2][0]*xp[0]-kp[2][1]*xp[1];
                iz+=kp[2][0]*xp[1]+kp[2][1]*xp[0];
            }
            sx+=rx*yp[0]-ix*yp[1];
            jx+=rx*yp[1]+ix*yp[0];
            sy+=ry*yp[0]-iy*yp[1];
            jy+=ry*yp[1]+iy*yp[0];
            sz+=rz*yp[0]-iz*yp[1];
            jz+=rz*yp[1]+iz*yp[0];
        }
        double zcos=cos(uk*hz),zsin=sin(uk*hz);
        ux+=sx*zcos-jx*zsin;
        uy+=sy*zcos-jy*zsin;
        uz+=sz*zcos-jz*zsin;
    }
}

/** Evaluates the velocity at multiple positions. This function assumes that
 * the function allocate_vel_table has been called beforehand, to ensure that
 * the number of table slots is greater or equal to the number of positions.
 * \param[in] q the number of positions.
 * \param[in] pos a pointer to the positions.
 * \param[in] vel a pointer to the velocities. */
void turb_fluid::vel_multi(int q,double *pos,double *vel) {
    fourier_table(q,pos);

    // Loop over all Fourier modes and calculate their contribution to the
    // velocity
    for(int i=0;i<3*q;i++) vel[i]=0;
#pragma omp parallel
    {
        double *r=rtab[thread_num()],*s=r+6*q;
#pragma omp for
        for(int k=0;k<o;k++) {
            int uk=k>o/2?k-o:k;
            fftw_complex *kp=kk+3*n*fftm*k;
            double *yp=ftab+2*fftm*q;
            for(int l=0;l<6*q;l++) s[l]=0;
            for(int j=0;j<n;j++) {

                // Compute contributions from a single x row of gridpoints
                for(int l=0;l<6*q;l++) r[l]=0;
                double *xp=ftab;
                for(fftw_complex *ke=kp+3*fftm;kp<ke;kp+=3)
                    for(double *rp=r;rp<s;xp+=2,rp+=6) xfill(rp,kp,xp);

                // Sum the row contributions to the xy-slice
                for(double *rp=r,*sp=s;rp<s;rp+=6,sp+=6,yp+=2) yfill(sp,rp,yp);
            }

            // Sum the xy-slice contributions into the final result
            for(double *sp=s,*hp=htab,*vp=vel;hp<htab+q;sp+=6,hp++,vp+=3)
                zfill(vp,sp,cos(*hp*uk),sin(*hp*uk));
        }
    }
}

/** Evaluates the velocity at multiple positions. This function assumes that
 * the function allocate_vel_table has been called beforehand, to ensure that
 * the number of table slots is greater or equal to the number of positions.
 * \param[in] q the number of positions.
 * \param[in] pos a pointer to the positions.
 * \param[in] vel a pointer to the velocities. */
void turb_fluid::vel_dot_multi(int q,double *pos,double *dir,double *vel,double *dvel) {
    fourier_table(q,pos);

    // Loop over all Fourier modes and calculate their contribution to the
    // velocity
    for(int i=0;i<3*q;i++) vel[i]=dvel[i]=0;
#pragma omp parallel
    {
        double *r=rtab[thread_num()],*s=r+12*q;
#pragma omp for
        for(int k=0;k<o;k++) {
            int uk=k>o/2?k-o:k;
            fftw_complex *kp=kk+3*n*fftm*k;
            double *yp=ftab+2*fftm*q,kz=facz*uk,w=kz*kz;
            for(int l=0;l<12*q;l++) s[l]=0;
            for(int j=0;j<n;j++) {

                // Compute contributions from a single x row of gridpoints
                for(int l=0;l<12*q;l++) r[l]=0;
                double *xp=ftab,ky=facy*(j>n/2?j-n:j),ww=w+ky*ky,kx=0;
                for(fftw_complex *ke=kp+3*fftm;kp<ke;kp+=3,kx+=facx)
                    for(double *rp=r,*dirp=dir;rp<s;xp+=2,rp+=12,dirp+=3) {
                        fftw_complex dkp[3];
                        xfill(rp,kp,xp);
                        double re=-pow(ww+kx*kx,1./3.)*Cinv,
                               im=kx*dirp[0]+ky*dirp[1]+kz*dirp[2];
                        dkp[0][0]=re*kp[0][0]-im*kp[0][1];
                        dkp[0][1]=re*kp[0][1]+im*kp[0][0];
                        dkp[1][0]=re*kp[1][0]-im*kp[1][1];
                        dkp[1][1]=re*kp[1][1]+im*kp[1][0];
                        dkp[2][0]=re*kp[2][0]-im*kp[2][1];
                        dkp[2][1]=re*kp[2][1]+im*kp[2][0];
                        xfill(rp+6,dkp,xp);
                }

                // Sum the row contributions to the xy-slice
                for(double *rp=r,*sp=s;rp<s;rp+=12,sp+=12,yp+=2) {
                    yfill(sp,rp,yp);
                    yfill(sp+6,rp+6,yp);
                }
            }

            // Sum the xy-slice contributions into the final result
            for(double *sp=s,*hp=htab,*vp=vel,*dvp=dvel;hp<htab+q;sp+=12,hp++,vp+=3,dvp+=3) {
                double zcos=cos(*hp*uk),zsin=sin(*hp*uk);
                zfill(vp,sp,zcos,zsin);
                zfill(dvp,sp+6,zcos,zsin);
            }
        }
    }
}

void turb_fluid::fourier_table(int q,double *pos) {

    // Assemble the table of Fourier coefficients in the x direction,
    // taking into account the symmetry
    double *fp=ftab;
    for(int l=0;l<q;l++) {
        *(fp++)=1.;*(fp++)=0.;
        htab[l]=(pos[3*l]-ax)*facx;
    }
    for(int i=1;i<fftm-1;i++) {
        for(int l=0;l<q;l++) {
            *(fp++)=2*cos(i*htab[l]);
            *(fp++)=2*sin(i*htab[l]);
        }
    }
    if(m&1) {
        for(int l=0;l<q;l++) {
            *(fp++)=2*cos((fftm-1)*htab[l]);
            *(fp++)=2*sin((fftm-1)*htab[l]);
        }
    } else {
        for(int l=0;l<q;l++) {
            *(fp++)=cos((fftm-1)*htab[l]);
            *(fp++)=sin((fftm-1)*htab[l]);
        }
    }

    // Assemble the table of Fourier coefficients in the y direction
    for(int l=0;l<q;l++) {
        *(fp++)=1.;*(fp++)=0.;
        htab[l]=(pos[3*l+1]-ay)*facy;
    }
    for(int j=1;j<n;j++) {
        int uj=j>n/2?j-n:j;
        for(int l=0;l<q;l++) {
            *(fp++)=cos(uj*htab[l]);
            *(fp++)=sin(uj*htab[l]);
        }
    }
    for(int l=0;l<q;l++) htab[l]=(pos[3*l+2]-az)*facz;
}

/** Perform a stochastic update to the Fourier modes. If mode=0, then the
 * Ornstein-Uhlenbeck equations at each mode are updated. If mode=1, then the
 * Fourier modes are initialized in steady state.
 * \param[in] dt the timestep to use for the stochastic update. This is ignored
 *               if mode=1.*/
template<int mode>
void turb_fluid::update_random(double dt) {
    double b1=dt*Cinv,b2=fnor*sqrt(2*b1);
#pragma omp parallel
    {
        gsl_rng *r=rng[thread_num()];

        // Loop over all of the Fourier modes
#pragma omp for
        for(int k=0;k<o;k++) {
            fftw_complex *kp=kk+3*n*fftm*k;
            double w=sqr(facz*(k>o/2?k-o:k));
            for(int j=0;j<n;j++) {
                double ww=w+sqr(facy*(j>n/2?j-n:j));
                for(int i=0;i<fftm;i++,kp+=3) {

                    // Diagnostic line to check on the (0,0,0) mode
                    //if(i==0&&j==0&&k==0)
                    //    printf("%g %g %g %g %g %g\n",kp[0][0],kp[0][1],kp[1][0],kp[2][1],kp[2][0],kp[2][1]);

                    // Calculate the scale factor to apply to this mode. Skip
                    // if this mode is not used.
                    int fm=f_mode(i,j,k);
                    if(fm==0) {
                        kp[0][0]=kp[0][1]=0.;
                        kp[1][0]=kp[1][1]=0.;
                        kp[2][0]=kp[2][1]=0.;
                        continue;
                    }
                    double www=ww+sqr(facx*i);
                    if(mode==0) {

                        // Apply the stochastic update for the
                        // Ornstein-Uhlenbeck process
                        double c=b1*pow(www,1./3.),
                               sig=fm*b2/sqrt(sqrt(www)*www);
                        kp[0][0]+=-kp[0][0]*c+gsl_ran_gaussian_ziggurat(r,sig);
                        kp[0][1]+=-kp[0][1]*c+gsl_ran_gaussian_ziggurat(r,sig);
                        kp[1][0]+=-kp[1][0]*c+gsl_ran_gaussian_ziggurat(r,sig);
                        kp[1][1]+=-kp[1][1]*c+gsl_ran_gaussian_ziggurat(r,sig);
                        kp[2][0]+=-kp[2][0]*c+gsl_ran_gaussian_ziggurat(r,sig);
                        kp[2][1]+=-kp[2][1]*c+gsl_ran_gaussian_ziggurat(r,sig);
                    } else {

                        // Initialize the modes to match the steady state
                        // distribution
                        double sig=fm*fnor*pow(www,-11/12.);
                        kp[0][0]=gsl_ran_gaussian_ziggurat(r,sig);
                        kp[0][1]=gsl_ran_gaussian_ziggurat(r,sig);
                        kp[1][0]=gsl_ran_gaussian_ziggurat(r,sig);
                        kp[1][1]=gsl_ran_gaussian_ziggurat(r,sig);
                        kp[2][0]=gsl_ran_gaussian_ziggurat(r,sig);
                        kp[2][1]=gsl_ran_gaussian_ziggurat(r,sig);
                    }
                }
            }
        }
    }
}

/** Perform a stochastic update to the Fourier modes. If mode=0, then the
 * Ornstein-Uhlenbeck equations at each mode are updated. If mode=1, then the
 * Fourier modes are initialized in steady state.
 * \param[in] dt the timestep to use for the stochastic update. This is ignored
 *               if mode=1.*/
double turb_fluid::mean_rms() {
    double s=0;
#pragma omp parallel for reduction(+:s)
    for(int k=0;k<o;k++) {
        double w=sqr(facz*(k>o/2?k-o:k));
        for(int j=0;j<n;j++) {
            double ww=w+sqr(facy*(j>n/2?j-n:j));
            for(int i=0;i<fftm;i++) if(f_mode(i,j,k)!=0)
                s+=pow(ww+sqr(facx*i),-11/6.);
        }
    }
    return fnor*sqrt(12*s);
}

/** Updates the modes to their expected value after some duration has passed,
 * based on mean reversion in the Ornstein-Uhlenbeck model.
 * \param[in] T the duration to consider. */
void turb_fluid::mean_revert(double T) {
    double b1=T*Cinv;

    // Loop over all of the Fourier modes
#pragma omp for
    for(int k=0;k<o;k++) {
        fftw_complex *kp=kk+3*n*fftm*k;
        double w=sqr(facz*(k>o/2?k-o:k));
        for(int j=0;j<n;j++) {
            double ww=w+sqr(facy*(j>n/2?j-n:j));
            for(int i=0;i<fftm;i++,kp+=3) {

                // Calculate the scale factor to apply to this mode. Skip
                // if this mode is not used.
                int fm=f_mode(i,j,k);
                if(fm==0) {
                    kp[0][0]=kp[0][1]=0.;
                    kp[1][0]=kp[1][1]=0.;
                    kp[2][0]=kp[2][1]=0.;
                    continue;
                }
                double c=exp(-b1*pow(ww+sqr(facx*i),1./3.));
                kp[0][0]*=c;
                kp[0][1]*=c;
                kp[1][0]*=c;
                kp[1][1]*=c;
                kp[2][0]*=c;
                kp[2][1]*=c;
            }
        }
    }
}

/** Computes a histogram of the energy as a function of the wave number
 * magnitude.
 * \param[in] hi a pointer in which to store the histogram.
 * \param[in] nbin the number of bins of the histogram.
 * \param[out] hmax the maximum wave number magnitude, setting the upper range
 *                  of the histogram bins. */
void turb_fluid::histogram(double *hi,int nbin,double &hmax) {

    // Compute the maximum wave number magnitude and the normalizing
    // factor for binning the contributions. Clear the histogram bins.
    hmax=sqrt(sqr(facx*0.5*m)+sqr(facy*0.5*n)+sqr(facz*0.5*o));
    double dsp=nbin/hmax;
    for(int i=0;i<nbin;i++) hi[i]=0;

    // Loop over the Fourier modes
#pragma omp parallel for
    for(int k=0;k<o;k++) {
        fftw_complex *kp=kk+3*n*fftm*k;
        double w=sqr(facz*(k>o/2?k-o:k));
        for(int j=0;j<n;j++) {
            int ww=w+sqr(facy*(j>n/2?j-n:j));
            for(int i=0;i<fftm;i++,kp+=3) {

                // Calculate the scale factor to apply to this mode. Skip
                // if this mode is not used.
                int fm=f_mode(i,j,k);
                if(fm==0) continue;

                // Bin the contribution from this Fourier mode
                int b=static_cast<int>(dsp*sqrt(ww+sqr(facx*i)));
                if(b<0||b>=nbin) continue;
                double c=(i>0&&2*i<m?4:1)*complex_msq(kp);
#pragma omp atomic
                hi[b]+=c;
            }
        }
    }

    // Normalize the histogram results
    for(int i=0;i<nbin;i++) hi[i]*=dsp;
}

/** Compute the maximum timestep based on how fast the relaxation term in the
 * Ornstein-Uhlenbeck equation for the maximum wave number can be integrated.
 * \return The timestep. */
double turb_fluid::est_max_timestep() {
    return pow(0.25*(sqr(facx*m)+sqr(facy*n)+sqr(facz*o)),-1/3.)/Cinv;
}

/** Checks that the Fourier coefficient tables have a given number of slots,
 * and if not, extends them.
 * \param[in] slots the number of slots. */
void turb_fluid::allocate_vel_table(int slots,bool extended) {
    if(slots>fslots) {
        fslots=slots;
        delete [] htab;
        htab=new double[(2*(fftm+n)+1)*slots];
        ftab=htab+slots;
#pragma omp parallel
        {
            int t=thread_num();
            delete [] rtab[t];
            rtab[t]=new double[(extended?24:12)*slots];
        }
    }
}

/** Calculates the velocity statistics based on randomly sampling the velocity
 * field throughout the domain.
 * \param[in] nsamp tha number of samples to take.
 * \param[out] (ubar,vbar,wbar) the components of the mean velocity vector
 * \param[out] (urms,vrms,wrms) the root mean squared of the components of the
 *                              velocity vector. */
void turb_fluid::sample_vel_stats(int nsamp,double &ubar,double &vbar,double &wbar,double &urms,double &vrms,double &wrms) {
    double *samp=new double[6*nsamp],*svel=samp+3*nsamp;
    gsl_rng *r=rng[0];

    // Create random sample points
    for(double *sp=samp;sp<samp+3*nsamp;sp+=3) {
        *sp=gsl_ran_flat(r,ax,bx);
        sp[1]=gsl_ran_flat(r,ay,by);
        sp[2]=gsl_ran_flat(r,az,bz);
    }

    // Calculate the velocities at the random sample points
    allocate_vel_table(nsamp);
    vel_multi(nsamp,samp,svel);

    // Compute the velocity statistics from the random samples
    vel_stats_internal(svel,nsamp,ubar,vbar,wbar,urms,vrms,wrms);
    delete [] samp;
}

/** Calculates the velocity statistics of samples in an array.
 * \param[in] up a pointer to the start of the array.
 * \param[in] nsamp the number of samples in the array (arranged as (ux,uy,uz)
 *                  triplets).
 * \param[out] (ubar,vbar,wbar) the components of the mean velocity vector
 * \param[out] (urms,vrms,wrms) the root mean squared of the components of the
 *                              velocity vector. */
void turb_fluid::vel_stats_internal(double *up,int nsamp,double &ubar,double &vbar,double &wbar,double &urms,double &vrms,double &wrms) {
    double su=0,sv=0,sw=0,suu=0,svv=0,sww=0,fac=1./nsamp;
    for(double *sp=up;sp<up+3*nsamp;sp+=3) {
        su+=*sp;
        sv+=sp[1];
        sw+=sp[2];
        suu+=*sp*(*sp);
        svv+=sp[1]*sp[1];
        sww+=sp[2]*sp[2];
    }

    // Calculate mean and variance of the velocity components
    ubar=su*fac;
    vbar=sv*fac;
    wbar=sw*fac;
    urms=sqrt(suu*fac);
    vrms=sqrt(svv*fac);
    wrms=sqrt(sww*fac);
}

/** Creates a number of random samples of velocity.
 * \param[in] samp an array in which to store the sample positions and velocities. */
void turb_fluid::correl_init(double *samp,int nsamp) {

    // Create the random sample points
#pragma omp parallel
    {
        gsl_rng *r=rng[thread_num()];

#pragma omp for
        for(double *sp=samp;sp<samp+3*nsamp;sp+=3) {
            //*sp=dx*gsl_rng_uniform_int(r,m);
            //sp[1]=dy*gsl_rng_uniform_int(r,n);
            //sp[2]=dz*gsl_rng_uniform_int(r,o);
            *sp=gsl_ran_flat(r,ax,bx);
            sp[1]=gsl_ran_flat(r,ay,by);
            sp[2]=gsl_ran_flat(r,az,bz);
        }
    }

    // Calculate the velocities at the random sample points
    allocate_vel_table(nsamp);
    vel_multi(nsamp,samp,samp+3*nsamp);
}

/** Computes the correlation function.
 * \param[in] w an array for storing the correlation function.
 * \param[in] nbin the number of bins for the correlation function.
 * \param[in] mrad the maximum separation to calculate the correlation function to.
 * \param[in] samp an array of sample positions and velocities at a previous
 *                 (or current) time.
 * \param[in] nsamp the number of samples of velocity to do in each bin. */
void turb_fluid::correl_function(double *w,int nbin,double mrad,double *samp,int nsamp) {
    double h=mrad/(nbin-1),*samp2=new double[6*nsamp],
           *svel=samp+3*nsamp,*svel2=samp2+3*nsamp;

    for(int i=0;i<nbin;i++) {
        gsl_rng *r=rng[0];
        double saa=0,sab=0,sbb=0,dx,dy,dz,rmin=h*i,rlen;

        for(int j=0;j<3*nsamp;j+=3) {

            // Compute a random vector with length in [rmin,rmin+h], and
            // hence compute a second position that is displaced by this
            // vector from the first
            gsl_ran_dir_3d(r,&dx,&dy,&dz);
            rlen=rmin+h*gsl_rng_uniform(r);
            samp2[j]=samp[j]+rlen*dx;
            samp2[j+1]=samp[j+1]+rlen*dy;
            samp2[j+2]=samp[j+2]+rlen*dz;
        }

        // Compute velocities at the sample positions
        vel_multi(nsamp,samp2,svel2);

        // Compute the sums required for the correlation coefficient
        for(int j=0;j<3*nsamp;j+=3) {
            double *s1=svel+j,*s2=svel2+j;
            saa+=*s1*(*s1)+s1[1]*s1[1]+s1[2]*s1[2];
            sab+=*s1*(*s2)+s1[1]*s2[1]+s1[2]*s2[2];
            sbb+=*s2*(*s2)+s2[1]*s2[1]+s2[2]*s2[2];
        }

        // Store the correlation coefficient
        w[i]=sab/sqrt(saa*sbb);
    }

    // Remove temporary array for computing velocity samples
    delete [] samp2;
}

/** Saves the complete state of the modes and parameters to an open file
 * handle.
 * \param[in] fp the file handle to write to. */
void turb_fluid::save(FILE *fp) {

    // Write the mode sizes (m,n,o). In addition, include the derived
    // quantities fftm, mn, and mno, which can be used to check file integrity.
    fwrite(&m,sizeof(int),6,fp);

    // Write the coordinate ranges
    fwrite(&ax,sizeof(double),6,fp);

    // Write the wind parameters
    fwrite(&Cinv,sizeof(double),2,fp);

    // Write the mode information
    fwrite(kk,sizeof(fftw_complex),3*fftm*n*o,fp);
}

// Explicit instantiation
template void turb_fluid::update_random<0>(double);
template void turb_fluid::update_random<1>(double);
