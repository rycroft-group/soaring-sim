#ifndef REINF_LEARN_HH
#define REINF_LEARN_HH

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <limits>

// When choosing the best action, it is possible that some actions may have
// equal value. This constant sets a tolerance for when the actions should be
// treated as equivalent.
const double reinf_learn_equiv_factor=10.*std::numeric_limits<double>::epsilon();

class reinf_learn {
    public:
        /** The total number of states. */
        const int ns;
        /** The total number of actions. */
        const int na;
        /** The total number of (state,action) pairs. */
        const int nsa;
        /** The learning rate. */
        const double alpha;
        /** The discount factor. */
        const double gamma;
        /** The expected value of each (state,action) pair. */
        double* const Q;
        /** The frequency visiting a state. */
        int* const Freq;
        const char *filename;
        reinf_learn(int ns_,int na_,double alpha_,double gamma_);
        ~reinf_learn();
        void reset_Q(double Qinit=0);
        void trained_Q();
        void shuffle_Q();
        void straight_Q(double directness=1.1);
        void straight_Q2();
        void load_Q();
        double max_Q(int s);
        double min_Q(int s);
        void print_Q();
        int pick_best_with_random(int s,double eps);
        int pick_best_with_prob(int s,double eps);
        int pick_best_uct(int s,double eps,int n_step);
        void update_uct();
        int pick_best(int s);
        void write_Q(const char* qfile);
        void write_Freq(const char* file);
        inline void update(int s,int a,int new_s,double rwd,int k) {
            Q[na*s+a]+=alpha*(rwd+gamma*max_Q(new_s)-Q[na*s+a]);}
        // inline void update(int s,int a,int new_s,double rwd,int k) {
        //     Q[na*s+a]+=(rwd+gamma*max_Q(new_s)-Q[na*s+a])*(1.0/k);} // use learning rate=1/t
    private:
        /** A temporary array for storing the best actions
         * for a given state. */
        int *ba;
        /** A pointer to the GSL random number generator. */
        gsl_rng *rng;
};

#endif
