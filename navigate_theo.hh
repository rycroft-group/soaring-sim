#ifndef NAVIGATE_HH
#define NAVIGATE_HH

#include <cstdio>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "reinf_learn.hh"

// Threshold on the concentration gradient for when to treat it as being zero
const double navigate_zero_grad_threshold=1e-6;
const double navigate_zgt_sq=navigate_zero_grad_threshold
                            *navigate_zero_grad_threshold;

class navigate_theo {
    public:
        /** The horizontal grid size. */
        const int m;
        /** The vertical grid size. */
        const int n;
        /** The number of blocks */
        const int n_block;
        /** The size of a subspce */
        const int block_size;
        /** The number of frames in dynamic concentration field. */
        const int n_frame;
        /** The total number of grid points. */
        const int mn;
        /** The total number of points(x,y,t). */
        const int mnf;
        /** The horizontal grid spacing. */
        const double dx;
        /** The vertical grid spacing. */
        const double dy;
        /** The output directory filename. */
        const char *filename;
        /** The first odor concentration. */
        double* const c1;
        /** The first reward field. */
        double* const r1;
        /** The first odor concentration. */
        double* const c2;
        /** The first reward field. */
        double* const r2;
        /** The first odor concentration. */
        double* const c3;
        /** The first reward field. */
        double* const r3;
        navigate_theo(int m_,int n_,int n_block_,int n_frame,double sx,double sy,double alpha,double gamma,const char *filename);
        ~navigate_theo();
        void init_conc_and_reward();
        void init_conc_and_reward2();
        void init_conc_and_reward3();
        void init_conc_and_reward4();
        double* read_cfile(const char* cfilename);
        int state(int i,int j,int k,int port);
        double training_episode(int i,int j,int steps,double eps,FILE *fp,int port);
        double test_episode(int i,int j,int steps,double eps,FILE *fp,int port);
        void grad_descent(int episodes,int steps,int port);
        void train(int episodes,int steps,double eps,int freq=0);
        int discrete_conc_grad(int i,int j,int k,int port);
        int weak_odor_cue(int i,int j,int k,int port);
        int strong_odor_cue(int i,int j,int k,int port);
        void write_fields();
        inline void print_Q() {rl.print_Q();}
        inline void write_Q(const char *qfile) {rl.write_Q(qfile);}
        inline void write_Freq(const char *qfile) {rl.write_Freq(qfile);}
       // inline void write_Q(FILE *fp) {rl.write_Q(FILE *fp);}
    private:
        void output_field(const int mode);
        inline bool out_of_bounds(int i,int j,int a);
        inline void move(int &i,int &j,int a);
        void set_conc_extrema();
        /** A buffer for assembling output filenames. */
        char* const buf;
        /** The minimum of the concentration field. */
        double cmin;
        /** The maximum of the concentration field. */
        double cmax;
        /** A pointer to the GSL random number generator. */
        gsl_rng *rng;
        /** The reinforcement learning class. */
        reinf_learn rl;
};

#endif
