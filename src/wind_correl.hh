#ifndef WIND_CORREL_HH
#define WIND_CORREL_HH

#include <cstdio>
#include <cmath>

#include "glider.hh"
#include "tf_grid.hh"

/** The total number of timepoints for measuring the wind correlation. */
const int wc_pts=7;

/** A data structure for computing Pearson's correlation coefficient between
 * two different sets of values. */
struct p_cor_coeff {
    /** The number of samples that have been taken. */
    long n;
    /** The sum of the x values. */
    double sx;
    /** The sum of the y values. */
    double sy;
    /** The sum of the squares of x values. */
    double sxx;
    /** The sum of the products of the x and y values. */
    double sxy;
    /** The sum of the squares of y values. */
    double syy;
    /** The container constructor leaves the internal variables uninitialized.
     */
    p_cor_coeff() {}
    /** Initializes the internal variables in preparation for measuring (x,y)
     * pairs. */
    inline void init() {
        n=0;
        sx=sy=sxx=sxy=syy=0;
    }
    /** Adds a measurement pair to the calculation.
     * \param[in] (x,y) the pair of measured values. */
    inline void add(double x,double y) {
        n++;
        sx+=x;
        sy+=y;
        sxx+=x*x;
        sxy+=x*y;
        syy+=y*y;
    }
    /** Combines the variables of another p_cor_coeff class into this one.
     * \param[in] pc a reference to the other class. */
    inline void combine(p_cor_coeff &pc) {
        n+=pc.n;
        sx+=pc.sx;
        sy+=pc.sy;
        sxx+=pc.sxx;
        sxy+=pc.sxy;
        syy+=pc.syy;
    }
    /** Computes Pearson's correlation coefficient. */
    inline double coeff() {
        double ninv=1./n,ex=sx*ninv,ey=sy*ninv;
        return (sxy*ninv-ex*ey)/sqrt((sxx*ninv-ex*ex)*(syy*ninv-ey*ey));
    }
    void diagnostic(const char* s,FILE *fp=stdout);
};

/** A data structure for computing the correlation between the actual wind
 * field and the GPR prediction. */
class wind_correl {
    public:
        /** The number of gliders per trial. */
        const int gpt;
        /** The number of MCTS trials to sample. */
        const int n_mcts;
        /** The total number of measurements for a given MCTS instance. */
        const int totm;
        /** A pointer to the class with the turbulent wind field. */
        turb_fluid_grid* const tf;
        /** A table of glider positions during an MCTS computation. */
        double* const ptab;
        /** A table of vertical wind velocities predicted using GPR during an
         * MCTS computation. */
        double* const wtab;
        /** A table for computing actual vertical wind velocities. */
        double* const mztab;
        /** Counters for the number of */
        int* const ctab;
        wind_correl(int n_mcts_,int gpt_,turb_fluid_grid* tf_);
        ~wind_correl();
        inline void reset_counters() {
            for(int i=0;i<wc_pts*gpt;i++) ctab[i]=0;
        }
        inline void measure(int j,int step,int id,glider &g_,double wz) {

            // Determine whether this step needs to be recorded or not
            int co;
            if(j==0) {
                if(step==0) co=0;
                else if(step==1) co=4;
                else if(step==2) co=5;
                else if(step==5) co=6;
                else return;
            } else if(step==0) {
                if(j==tenth) co=1;
                else if(j==fifth) co=2;
                else if(j==half) co=3;
                else return;
            } else return;

            // Store the measurement
            measure_internal(co,id,g_,wz);
        }
        /** Sets the markers used to find the timepoints corresponding
         * to 0.1, 0.2, and 0.5 of a control duration.
         * \param[in] substeps the total number of integration substeps in a
         *                     control duration. This should always be a
         *                     multiple of 10. */
        inline void set_substep_markers(int substeps) {
            if(substeps%10!=0) fatal_error("Substep problem",1);
            tenth=substeps/10;
            fifth=2*tenth;
            half=5*tenth;
        }
        void wind_computation();
        void output(FILE *fp);
    private:
        void measure_internal(int co,int id,glider &g_,double wz);
        /** The correlation coefficient classes. */
        p_cor_coeff wc[wc_pts];
        /** The number of substeps corresponding to a tenth of a step. */
        int tenth;
        /** The number of substeps corresponding to a fifth of a step. */
        int fifth;
        /** The number of substeps corresponding to half of a step. */
        int half;
};

#endif
