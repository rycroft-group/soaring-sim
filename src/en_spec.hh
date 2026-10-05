#ifndef EN_SPEC_HH
#define EN_SPEC_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

/** A data structure for binning information in the mode energy spectrum. */
struct en_spec_data {
    /** The number of evaluations of E(k) that were made in this bin. */
    long n;
    /** The wavevector k. */
    double k;
    /** The energy E(k). */
    double E;
    /** Sets the counters in the data structure to zero. */
    inline void clear() {n=0;k=0;E=0;}
    /** Adds the contribution of a (k,E(k)) pair to the counters.
     * \param[in] (k_,E_) the data to add. */
    inline void contrib(double k_,double E_) {
#pragma omp atomic
        n++;
        double lk=log(k_);
#pragma omp atomic
        k+=lk;
#pragma omp atomic
        E+=E_;
    }
    /** Processes the previous counted data to obtain the final (k,E(k)) for
     * analysis.
     * \param[in] nor the normalization factor to apply to the E measurement. */
    inline void process(double nor) {if(n>0) {k=exp(k/n);E*=nor;}}
};

/** Class for representing the parameters and nonlinear mapping functions used
 * to compute a mode energy spectrum. */
class en_spec_param {
    public:
        /** The number of bins per decade for large wavenumbers. */
        const double dbin;
        /** The geometric mean of mode number at which to position
         * the linear-to-logarithmic bin transition. */
        const double ltl_k;
        /** The number of bins to use. */
        int nbin;
        /** The class constructor sets that parameters that control
         * how many bins to use in the energy spectrum calculation.
         * \param[in] dbin_ the number of bins per decade for large wavenumbers.
         * \param[in] ltl_k_ the geometric mean of mode number at which to position
         *                  the linear-to-logarithmic bin transition. */
        en_spec_param(double dbin_,double ltl_k_) : dbin(dbin_), ltl_k(ltl_k_),
            Enor(NULL) {}
        /** The destructor frees the dynamically allocated normalization table,
         * if in use. */
        ~en_spec_param() {
            if(Enor!=NULL) delete [] Enor;
        }
        void init(double kmin,double kgeo,double kmax);
        /** Clears an energy spectrum array. */
        void clear(en_spec_data *ed) {
            for(int i=0;i<nbin;i++) ed[i].clear();
        }
        /** Bins a (k,E(k)) into an energy spectrum array.
         * \param[in] ed a pointer to the energy spectrum array.
         * \param[in] (k,E) the data to bin. */
        inline void bin(en_spec_data *ed,double k,double E) {
            int i=int(f(k));
            if(i<0||i>=nbin) {
                fprintf(stderr,"Error: k=%g, E=%g, i=%d\n",k,E,i);
                exit(1);
            }
            ed[i].contrib(k,E);
        }
        /** Evaluates hyperbolic arcsine mapping function.
         * \param[in] x the function argument.
         * \return The function result. */
        inline double f(double x) {
            return s*asinh(x*ir)-flo;
        }
        /** Evaluates the inverse mapping function.
         * \param[in] y the inverse function argument.
         * \return The inverse function result. */
        inline double f_inv(double y) {
            return r*sinh((y+flo)*is);
        }
        /** Completes the computation of the energy spectrum, by normalizing
         * E(k) contributions, and computing the geometric mean wavenumber k
         * per bin.
         * \param[in] ed the array of energy spectrum data.
         * \return The integral of E(k). */
        inline double finalize(en_spec_data *ed) {
            double E_integ=0;
            for(int i=0;i<nbin;i++) {
                E_integ+=ed[i].E;
                ed[i].process(Enor[i]);
            }
            return E_integ;
        }
    private:
        /** The first scaling factor used in the asinh mapping function. */
        double r;
        /** The reciprocal of the first scaling factor used in the asinh
         * mapping function. */
        double ir;
        /** The second scaling factor used in the asinh mapping function. */
        double s;
        /** The reciprocal of the second scaling factor used in the asinh
         * mapping function. */
        double is;
        /** The offset to apply to the scaling function, based on the lowest
         * value that f can evaluate to. */
        double flo;
        /** The array of normalization factors to apply to the energy spectrum
         * data. */
        double *Enor;
};

#endif
