#include "tf_grid.hh"

#include "common.hh"
#include <cstdio>

/** Initializes the three-dimensional turbulent fluid generator with additional
 * FFTW-based routines for evaluating the fluid on a grid.
 * \param[in] (m_,n_,o_) the dimensions of the grid.
 * \param[in] (ax_,bx_) the lower and upper x-coordinate simulation bounds.
 * \param[in] (ay_,by_) the lower and upper y-coordinate simulation bounds.
 * \param[in] (az_,bz_) the lower and upper z-coordinate simulation bounds.
 * \param[in] Cinv_ the reciprocal of the constant controlling mode timescales.
 * \param[in] alpha_ the constant controlling mode energy scales.
 * \param[in] seed a base seed for initializing the GSL random number
 *                 generators. */
turb_fluid_grid::turb_fluid_grid(int m_,int n_,int o_,double ax_,double bx_,double ay_,double by_,double az_,double bz_,double Cinv_,double alpha_,unsigned long seed)
    : turb_fluid(m_,n_,o_,ax_,bx_,ay_,by_,az_,bz_,Cinv_,alpha_,seed),
    xsp(1./dx), ysp(1./dy), zsp(1./dz),
    ksize(sizeof(fftw_complex)*3*fftm*n*o),
    kcopy((fftw_complex*)fftw_malloc(ksize)),
    uu((double*)fftw_malloc(sizeof(double)*3*mno)),
    plan_bck(fftw_plan_many_dft_c2r(3,(const int*) this,3,kcopy,NULL,3,1,uu,NULL,3,1,FFTW_MEASURE)),
    a_c(turb_fluid_grid_a_c_default) {}

/** The class destructor destroys the FFTW plans and frees the dynamically
 * allocated memory. */
turb_fluid_grid::~turb_fluid_grid() {
    fftw_destroy_plan(plan_bck);
    fftw_free(uu);
    fftw_free(kcopy);
}

/** Outputs a cross-section of the velocity field.
 * \param[in] filename the name of the file the write to.
 * \param[in] uu a pointer to the (0,0) element of the cross section to write.
 * \param[in] a1 the lower limit in the first dimension.
 * \param[in] d1 the step size in the first dimension.
 * \param[in] a2 the lower limit in the second dimension.
 * \param[in] d2 the step size in the second dimension.
 * \param[in] l1 the number of grid points to write in the first dimension.
 * \param[in] l2 the number of grid points to write in the second dimension.
 * \param[in] s1 the memory shift length in the first dimension.
 * \param[in] s2 the memory shift length in the second dimension. */
void turb_fluid_grid::output_cross_section(const char* filename,double *uu,double a1,double d1,double a2,double d2,int l1,int l2,int s1,int s2) {

    // Open the output file and create temporary memory
    FILE *outf=safe_fopen(filename,"wb");
    float *buf=new float[l1+1];

    // Output the first line of the file
    int i,j;
    float *bp=buf+1,*be=bp+l1;
    *buf=l1;
    for(i=0;i<l1;i++) *(bp++)=a1+d1*i;
    fwrite(buf,sizeof(float),l1+1,outf);

    // Output the field values to the file
    double *ur=uu;
    for(j=0;j<l2;j++,ur+=s2) {
        double *up=ur;
        *buf=a2+d2*j;bp=buf+1;
        while(bp<be) {*(bp++)=*up;up+=s1;}
        fwrite(buf,sizeof(float),l1+1,outf);
    }

    // Close the file and free the dynamically allocated memory
    fclose(outf);
    delete [] buf;
}

/** Computes a full discrete Fourier transform of the velocity field, using
 * direct evaluation at every point. This is a computationally expensive
 * routine, and is used for diagnostic purposes to check it agrees with the
 * FFTW output. */
void turb_fluid_grid::full_dft() {
    double *up=uu;
    for(int k=0;k<o;k++) {
        double z=az+k*dz;
        for(int j=0;j<n;j++) {
            double y=ay+j*dy;
            for(int i=0;i<o;i++,up+=3)
                vel(ax+i*dx,y,z,*up,up[1],up[2]);
        }
    }
}

/** Performs trilinear interpolation of the fluid velocity on the grid.
 * \param[in] (x,y,z) the position at which to interpolate.
 * \param[out] (ux,uy,uz) the velocity vector. */
void turb_fluid_grid::lin_interp(double x,double y,double z,double &ux,double &uy,double &uz) {

    // Compute the box in the primary domain that this point is within, and the
    // fractional position of the point within the box
    int i,j,k;
    grid_remap(x,y,z,i,j,k);

    // Calculate the displacements to the upper indices
    int xd=3*(i==m-1?1-m:1),
        yd=3*m*(j==n-1?1-n:1),
        zd=3*m*n*(k==o-1?1-o:1);

    // Calculate pointer to the lower corner of the box
    double *up=uu+3*(i+m*(j+n*k)),gx=1-x,gy=1-y,gz=1-z,
           *uq=up+zd;

    // Calculate the linear interpolants
    ux=gz*(gy*(gx*up[0]+x*up[xd])+y*(gx*up[yd]+x*up[xd+yd]))
      +z*(gy*(gx*uq[0]+x*uq[xd])+y*(gx*uq[yd]+x*uq[xd+yd]));
    up++;uq++;
    uy=gz*(gy*(gx*up[0]+x*up[xd])+y*(gx*up[yd]+x*up[xd+yd]))
      +z*(gy*(gx*uq[0]+x*uq[xd])+y*(gx*uq[yd]+x*uq[xd+yd]));
    up++;uq++;
    uz=gz*(gy*(gx*up[0]+x*up[xd])+y*(gx*up[yd]+x*up[xd+yd]))
      +z*(gy*(gx*uq[0]+x*uq[xd])+y*(gx*uq[yd]+x*uq[xd+yd]));
}

/** Performs interpolation of the fluid velocity on the grid using a four-point
 * stencil in each direction (either tricubic or Lanczos2).
 * \param[in] cub true for tricubic interpolation, false for Lanczos2
 *                interpolation.
 * \param[in] (x,y,z) the position at which to interpolate.
 * \param[out] (ux,uy,uz) the velocity vector. */
void turb_fluid_grid::four_pt_interp(bool cub,double x,double y,double z,double &ux,double &uy,double &uz) {

    // Determine which box of the fluid grid we are in
    int i,j,k;
    grid_remap(x,y,z,i,j,k);

    // Compute coefficient tables for the four-point interpolant in the three
    // coordinates. In addition, compute the memory stride lengths, which must
    // account for the periodicity of the domain.
    double s[12];
    int q[12];
    cub?setup_cubic_basis(x,y,z,s)
       :setup_lanczos2_basis(x,y,z,s);
    setup_memory_strides(i,j,k,q);

    // Compute the tricubic interpolant value using the 4x4x4 grid of nearby
    // values
    double *up=uu+3*(i+m*(j+n*k));
    ux=uy=uz=0;
    for(int l=0;l<4;l++) interp_slice(ux,uy,uz,s,q,up+q[8+l],s[8+l]);
}

/** Sums up contributions to a tricubic interpolant of the velocity in a 4x4
 * slice of points.
 * \param[in,out] (ux,uy,uz) the velocity to add these contributions to.
 * \param[in] s a pointer to the interpolant coefficients in the x and y
 *              directions.
 * \param[in] q a pointer to the memory stride lengths.
 * \param[in] up a pointer to the base gridpoint for interpolation.
 * \param[in] sf a scale factor to apply to the results. */
void turb_fluid_grid::interp_slice(double &ux,double &uy,double &uz,double *s,int *q,double *up,double sf) {
    double vx=0,vy=0,vz=0;
    interp_line(vx,vy,vz,s,q,up+q[4],s[4]);
    interp_line(vx,vy,vz,s,q,up+q[5],s[5]);
    interp_line(vx,vy,vz,s,q,up+q[6],s[6]);
    interp_line(vx,vy,vz,s,q,up+q[7],s[7]);
    ux+=sf*vx;
    uy+=sf*vy;
    uz+=sf*vz;
}

/** Sums up contributions to a tricubic interpolant of the velocity in a
 * 4-gridpoint line in the x direction.
 * \param[in,out] (vx,vy,vz) the velocity to add these contributions to.
 * \param[in] s a pointer to the interpolant coefficients in the x direction.
 * \param[in] q a pointer to the memory stride lengths.
 * \param[in] up a pointer to the base gridpoint for interpolation.
 * \param[in] sf a scale factor to apply to the results. */
inline void turb_fluid_grid::interp_line(double &vx,double &vy,double &vz,double *s,int *q,double *up,double sf) {
    vx+=sf*(*s*up[*q]+s[1]*up[q[1]]+s[2]*up[q[2]]+s[3]*up[q[3]]);
    vy+=sf*(*s*up[1+*q]+s[1]*up[1+q[1]]+s[2]*up[1+q[2]]+s[3]*up[1+q[3]]);
    vz+=sf*(*s*up[2+*q]+s[1]*up[2+q[1]]+s[2]*up[2+q[2]]+s[3]*up[2+q[3]]);
}

/** Sets up the cubic basis coefficients in each of the three coordinate
 * directions.
 * \param[in] (x,y,z) the fractional position of a point within a grid box at
 *                    which to interpolate.
 * \param[out] s an array of length 12 for storing the coefficients (4 each for
 *              the three coordinate directions). */
void turb_fluid_grid::setup_cubic_basis(double x,double y,double z,double *s) {
    cubic_coeffs(x,s);
    cubic_coeffs(y,s+4);
    cubic_coeffs(z,s+8);
}

/** Sets up the Lanczos2 basis coefficients in each of the three coordinate
 * directions.
 * \param[in] (x,y,z) the fractional position of a point within a grid box at
 *                    which to interpolate.
 * \param[out] s an array of length 12 for storing the coefficients (4 each for
 *              the three coordinate directions). */
void turb_fluid_grid::setup_lanczos2_basis(double x,double y,double z,double *s) {
    lanczos2_coeffs(x,s);
    lanczos2_coeffs(y,s+4);
    lanczos2_coeffs(z,s+8);
}

/** Calculates the mean of the square components of the velocity field
 * using the grid-based values.
 * \param[out] (mux,muy,muz) the mean of the square components of the velocity.
 */
void turb_fluid_grid::mean_squared(double &mux,double &muy,double &muz) {
    mux=muy=muz=0;
    for(double *up=uu;up<uu+3*mno;up+=3) {
        mux+=*up*(*up);
        muy+=up[1]*up[1];
        muz+=up[2]*up[2];
    }
    double fac=1./mno;
    mux*=fac;muy*=fac;muz*=fac;
}
