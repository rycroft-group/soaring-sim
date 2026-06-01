#include "k_func.hh"

#include "common.hh"
#include <cstdio>

/** A function to link to the gsl_function class in order to perform the kernel
 * integration with (r,t) dependence using the GSL integration routines.
 * \param[in] x the function argument.
 * \param[in] params a pointer to additional arguments. */
double kern_rt_integrand(double x,void *params) {
    double *p=reinterpret_cast<double*>(params),
           r=*p,t=p[1],Cinv=p[2],xr=x*r,
           tauinv=Cinv*pow(x,(2./3));
    return pow(x,(-5./3))*(xr==0?1:sin(xr)/(xr))*exp(-t*tauinv);
}

/** A function to link to the gsl_function class in order to perform the kernel
 * integration with r dependence using the GSL integration routines.
 * \param[in] x the function argument.
 * \param[in] params a pointer that is ignored. */
double kern_r_integrand(double x,void *params) {
    return pow(x,(-8./3))*sin(x);
}

/** Sets up the kernel evaluation function.
 * \param[in] N_ the number of turbulent fluid modes in one direction.
 * \param[in] B_ The turbulent fluid box size.
 * \param[in] rs_ The size of the grid in the r direction.
 * \param[in] rmax_ The maximum extent of the interpolation in the r direction.
 * \param[in] bl_ The size of the interpolation grid. */
kernel_func::kernel_func(int N_,double B_,int rs_,double rmax_,int blsize)
    : N(N_), rs(rs_), miss(0), rmax(rmax_), B(B_), kmin(sqrt(3)*M_PI/B),
    kmax(kmin*N), anor((2/3.)/(pow(kmin,(-2./3))-pow(kmax,(-2./3)))),
    rfac(rmax/((rs-1)*(rs-1))), irfac((rs-1)/sqrt(rmax)),
    bl(new double[blsize]) {}

/** Sets up the kernel evaluation function.
 * \param[in] N_ the number of turbulent fluid modes in one direction.
 * \param[in] B_ The turbulent fluid box size.
 * \param[in] Cinv_ The reciprocal of the turbulent fluid mode timescale.
 * \param[in] (rs_,ts_) The bilinear interpolation grid dimensions.
 * \param[in] (rmax_,tmax_) The maximum extent of the bilinear interpolation
 *                          grid. */
kernel_rt::kernel_rt(int N_,double B_,double Cinv_,int rs_,int ts_,double rmax_,double tmax_)
    : kernel_func(N_,B_,rs_,rmax_,rs_*ts_), ts(ts_), tmax(tmax_), Cinv(Cinv_),
    tfac(tmax/((ts-1)*(ts-1))), itfac((ts-1)/sqrt(tmax)) {

    // Create temporary table of r values for bilinear interpolation
    double *rv=new double[rs];
    for(int i=0;i<rs;i++) rv[i]=rstretch(i);

    // Fill in the bilinear interpolation table entries
#pragma omp parallel
    {

        // Set up the GSL integration workspace and parameters to pass to the
        // integrand function
        double params[3],res,err;
        gsl_integration_workspace *giw=gsl_integration_workspace_alloc(kernel_func_ws_size);
        gsl_function F;
        F.function=&kern_rt_integrand;
        F.params=params;
        params[2]=Cinv;

        // Loop over the bilinear interpolation grid, using different threads
        // to process different rows
#pragma omp for
        for(int j=0;j<ts;j++) {
            double *p=bl+rs*j,t=tstretch(j);
            params[1]=t;
            for(int i=0;i<rs;i++) {

                // Perform the adaptive integration, and scale the result
                *params=rv[i];
                gsl_integration_qags(&F,kmin,kmax,kernel_func_integ_tol,0,kernel_func_ws_size,giw,&res,&err);
                *(p++)=anor*res;
            }
        }
        gsl_integration_workspace_free(giw);
    }

    // Remove the temporary table
    delete [] rv;
}

/** Sets up the kernel evaluation function.
 * \param[in] N_ the number of turbulent fluid modes in one direction.
 * \param[in] B_ The turbulent fluid box size.
 * \param[in] C_ The turbulent fluid mode timescale.
 * \param[in] (rs_,ts_) The bilinear interpolation grid dimensions.
 * \param[in] (rmax_,tmax_) The maximum extent of the bilinear interpolation
 *                          grid. */
kernel_r::kernel_r(int N_,double B_,int rs_,double rmax_)
    : kernel_func(N_,B_,rs_,rmax_,rs_) {

    // Set up GSL integration
    gsl_integration_workspace *giw=gsl_integration_workspace_alloc(kernel_func_ws_size);
    gsl_function F;
    F.function=&kern_r_integrand;

    // Loop over the different r values to fill in the interpolation table
    *bl=1;
    for(int i=1;i<rs;i++) {
        double r=rstretch(i),res,err;

        // Perform integration. No parameters are needed since they can all be
        // scaled out of the integrand.
        gsl_integration_qags(&F,kmin*r,kmax*r,kernel_func_integ_tol,0,kernel_func_ws_size,giw,&res,&err);
        bl[i]=anor*pow(r,2/3.)*res;
    }
    gsl_integration_workspace_free(giw);
}

/** The class destructor frees the dynamically allocated memory. */
kernel_func::~kernel_func() {
    delete [] bl;
}

/** Outputs the table of pre-computed kernel function values.
 * \param[in] filename the name of the file to write to. */
void kernel_rt::output_table(const char* filename) {
    int i,j;

    // Open file and print error if there is a problem
    FILE *outf=safe_fopen(filename,"wb");
    if(outf==NULL) {
        fputs("Can't open file\n",stderr);
        exit(1);
    }

    // Allocate memory and write the header file
    float *fbuf=new float[rs+1],*fp=fbuf;
    *(fp++)=rs;for(i=0;i<rs;i++) *(fp++)=rstretch(i);
    fwrite(fbuf,sizeof(float),rs+1,outf);

    // Write field entries line-by-line
    double *p=bl;
    for(j=0;j<ts;j++) {

        // Write header entry
        fp=fbuf;*(fp++)=tstretch(j);

        // Write a horizontal line to the buffer
        for(i=0;i<rs;i++) *(fp++)=*(p++);
        fwrite(fbuf,sizeof(float),rs+1,outf);
    }

    // Remove temporary memory and close file
    delete [] fbuf;
    fclose(outf);
}

/** Evaluates the kernel function using the previously computed table.
 * \param[in] r the separation.
 * \param[in] t the time difference.
 * \return The kernel function. */
double kernel_rt::eval(double r,double t) {

    // Perform the nonlinear transformation on the input values
    double qr=irstretch(r),
           qt=itstretch(t);

    // Identify the grid square that the result is in
    int i=static_cast<int>(qr),j=static_cast<int>(qt);

    // If i and j are negative, then this is an error
    if(i<0||j<0) {
        fprintf(stderr,"Out of bounds: (r,t)=(%g,%g)\n",r,t);
        exit(1);
    }

    // If i and j are too large, then just return K=0 with a warning
    if(i>rs-1||j>ts-1) {
        if(kernel_func_warn_too_large)
            printf("# Out of range: %d %d (r,t)=(%g,%g) [performing integral]\n",i,j,r,t);
#pragma omp atomic
        miss++;
        return integrate(r,t);
    }

    // Return the bilinear interpolation of values in the table
    double *p=bl+i+rs*j;
    qr-=i;qt-=j;
    return (*p*(1-qr)+p[1]*qr)*(1-qt)+(p[rs]*(1-qr)+p[rs+1]*qr)*qt;
}

/** Evaluates the kernel function using the previously computed table.
 * \param[in] r the separation.
 * \param[in] t the time difference.
 * \return The kernel function. */
double kernel_r::eval(double r,double t) {

    // Perform the nonlinear transformation on the input value
    double qr=irstretch(r);

    // Identify the grid square that the result is in. If i is negative, then
    // this is an error.
    int i=static_cast<int>(qr);
    if(i<0) {
        fprintf(stderr,"Out of bounds: r=%g\n",r);
        exit(1);
    }

    // If i is too large, then call
    if(i>rs-1) {
        if(kernel_func_warn_too_large)
            printf("# Out of range: i=%d r=%g [performing integral]\n",i,r);
#pragma omp atomic
        miss++;
        return integrate(r);
    }

    // Return the bilinear interpolation of values in the table
    qr-=i;
    return bl[i]*(1-qr)+bl[i+1]*qr;
}

/** Calculates the kernel function using adaptive integration. This routine
 * designed for testing purposes, since allocates and deallocates GSL
 * integration workspace each time it is called.
 * \param[in] r the separation.
 * \param[in] t the time difference.
 * \return The kernel function. */
double kernel_rt::integrate(double r,double t) {

    // Set up GSL workspace and function
    double res,err,params[3];
    gsl_integration_workspace *giw=gsl_integration_workspace_alloc(kernel_func_ws_size);
    gsl_function F;
    F.function=&kern_rt_integrand;

    // Set up parameters to pass to integrand function
    F.params=params;
    *params=r;
    params[1]=t;
    params[2]=Cinv;

    // Perform the adaptive integration and free the workspace
    gsl_integration_qags(&F,kmin,kmax,kernel_func_integ_tol,0,kernel_func_ws_size,giw,&res,&err);
    gsl_integration_workspace_free(giw);
    return anor*res;
}

/** Calculates the kernel function using adaptive integration. This routine
 * designed for testing purposes, since allocates and deallocates GSL
 * integration workspace each time it is called.
 * \param[in] r the separation.
 * \return The kernel function. */
double kernel_r::integrate(double r) {
    if(r==0) return 1;

    // Set up GSL workspace and function
    double res,err;
    gsl_integration_workspace *giw=gsl_integration_workspace_alloc(kernel_func_ws_size);
    gsl_function F;
    F.function=&kern_r_integrand;

    // Perform the adaptive integration and free the workspace
    gsl_integration_qags(&F,kmin*r,kmax*r,kernel_func_integ_tol,0,kernel_func_ws_size,giw,&res,&err);
    gsl_integration_workspace_free(giw);
    return anor*pow(r,2/3.)*res;
}

/** Outputs the table of pre-computed kernel function values.
 * \param[in] filename the name of the file to write to. */
void kernel_r::output_table(const char* filename) {

    // Open file and print error if there is a problem
    FILE *outf=safe_fopen(filename,"w");
    if(outf==NULL) {
        fputs("Can't open file\n",stderr);
        exit(1);
    }

    // Output the table values and close the file
    for(int i=0;i<rs;i++) fprintf(outf,"%.10g %.10g\n",rstretch(i),bl[i]);
    fclose(outf);
}
