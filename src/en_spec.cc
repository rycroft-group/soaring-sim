#include "en_spec.hh"

/** Initializes the hyperbolic arcsine mapping function.
 * \param[in] kmin the minimum wavenumber in the turbulent fluid class.
 * \param[in] kgeo the geometric mean of the smallest wavenumbers along
 * each coordinate direction.
 * \param[in] kmax the maximum wavenumber in the turbulent fluid class. */
void en_spec_param::init(double kmin,double kgeo,double kmax) {

    // Initialize parameters for the nonlinear scaling function and bin size
    r=ltl_k*kgeo;ir=1./r;
    s=(1./log(10))*dbin;is=1./s;

    // Check that the lower bin range is greater than zero, because otherwise
    // the bin choices may give unexpected results
    flo=s*asinh(kmin*ir)-0.5;
    if(flo<=0) {
        fputs("Lower bin range in energy spectrum is negative; adjust parameters\n",stderr);
        exit(1);
    }

    // Compute number of bins to use, and allocate memory
    nbin=int(f(kmax))+1;

    // Compute the normalizing factors to apply to the energy spectrum
    if(Enor!=NULL) delete [] Enor;
    Enor=new double[nbin];
    double k1=f_inv(0),k2;
    for(int i=0;i<nbin;i++) {
        k2=f_inv(i+1);
        Enor[i]=1/(k2-k1);
        k1=k2;
    }
}
