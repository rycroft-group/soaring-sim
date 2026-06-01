#ifndef STATS_HH
#define STATS_HH

#include <cmath>
#include <limits>

/** A data structure for storing basic statistics about a collection of values.
 */
struct cli_stats {
    /** The first moment of the values. */
    double fmom;
    /** The second moment of the values. */
    double smom;
    /** The minimum of the values. */
    double mi;
    /** The maximum of the values. */
    double ma;
    /** Initializes the statistical counters before processing the values. */
    void init() {
        fmom=smom=0;
        mi=std::numeric_limits<double>::max();
        ma=-std::numeric_limits<double>::max();
    }
    /** Calculates the mean and standard deviation of the values.
     * \param[in] nor a normalizing factor, equal to the reciprocal of the
     *                number of values.
     * \param[out] mu the mean.
     * \param[out] sig the standard deviation. */
    inline void mu_sig(double nor,double &mu,double &sig) {
        mu=fmom*nor;
        sig=sqrt(smom*nor-mu*mu);
    }
    /** Updates the statistics to account for a given value.
     * \param[in] x the value. */
    inline void contrib(double x) {
        fmom+=x;
        smom+=x*x;
        if(x<mi) mi=x;
        if(x>ma) ma=x;
    }
};

/** A data structure for storing basic statistics about a collection of values.
 */
struct mti_stats {
    /** The total number of measurements. */
    long n;
    /** The first moment of the values. */
    long fmom;
    /** The second moment of the values. */
    long smom;
    /** The minimum of the values. */
    int mi;
    /** The maximum of the values. */
    int ma;
    /** Initializes the statistical counters before processing the values. */
    void init() {
        n=fmom=smom=0;
        mi=std::numeric_limits<int>::max();
        ma=0;
    }
    /** Calculates the mean and standard deviation of the values.
     * \param[in] nor a normalizing factor, equal to the reciprocal of the
     *                number of values.
     * \param[out] mu the mean.
     * \param[out] sig the standard deviation. */
    inline void mu_sig(double &mu,double &sig) {
        double nor=1./n;
        mu=fmom*nor;
        sig=sqrt(smom*nor-mu*mu);
    }
    /** Updates the statistics to account for a given value.
     * \param[in] x the value. */
    inline void contrib(double x) {
        n++;
        fmom+=x;
        smom+=x*x;
        if(x<mi) mi=x;
        if(x>ma) ma=x;
    }
    /** Adds the contributions of another instance of this class to this one.
     * \param[in] m a reference of the class to consider. */
    inline void add(mti_stats &m) {
        n+=m.n;
        fmom+=m.fmom;
        smom+=m.smom;
        if(m.mi<mi) mi=m.mi;
        if(m.ma>ma) ma=m.ma;
    }
};

#endif
