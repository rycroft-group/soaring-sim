#ifndef FILEINFO_HH
#define FILEINFO_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>

#include "common.hh"
#include "glider.hh"

/** The size of the temporary buffer for parsing the input file. */
const int fileinfo_buf_size=512;

/** The padding for the buffer for assembling output filenames. */
const int fileinfo_fbuf_pad_size=256;

/** The type of model used for predicting the wind field. */
enum wind_model {
    wm_unset, wm_gpr, wm_full_linear, wm_full_cubic
};

/** The planning model type. */
enum planning_type {
    pt_unset, pt_zero, pt_random, pt_mcts
};

/** \brief A class for parsing the drop impact parameters from a text
 * configuration file. */
class fileinfo {
    public:
        /** The number of trials. */
        int num_trials;
        /** The number of gliders per trial. */
        int gpt;
        /** The number of Fourier modes in the x direction. */
        int nx;
        /** The number of Fourier modes in the y direction. */
        int ny;
        /** The number of Fourier modes in the z direction. */
        int nz;
        /** The number of memory entries in Gaussian process regression. */
        int gprm;
        /** The number of memory measurements to take per control event. */
        int gs_pc;
        /** The number of outputs to take per control event. */
        int out_pc;
        /** The number of MCTS trials to sample. */
        int n_mcts;
        /** The depth of the MCTS trials, corresponding to the number of
         * actions taken. */
        int mcts_depth;
        /** The MCTS exploration factor, specified as a dimensionless constant. */
        float mcts_ex_fac;
        /** The number of Hermite interpolation segments, used for interpolating
         * forward in time using the full information model. */
        int h_segs;
        /** The number of r control points in the GPR kernel map. */
        int gpr_rs;
        /** The number of t control points in the GPR kernel map. */
        int gpr_ts;
        /** The number of control durations between outputting all of the MCTS
         * paths. */
        int path_interval;
        /** For full information gliders, whether to bypass the mean-reversion
         * prediction of the wind field. */
        bool bypass_mr;
        /** Whether to always bypass the Woodbury formula in Gaussian process
         * regression and do a full LAPACK computation each time during
         * updates. */
        bool gpr_full_compute;
        /** Flags controlling the types of file output. */
        unsigned int fflags;
        /** The seed for the random number generator. */
        unsigned long base_seed;
        /** The length scale to convert from simulation units to physical
         * units. */
        double l_phys;
        /** The time scale to convert from simulation units to physical
         * units. */
        double t_phys;
        /** The velocity scale to convert from simulation units to physical
         * units. */
        double v_phys;
        /** The gravitational acceleration scale to convert from simulation
         * units to physical units. */
        double g_phys;
        /** The viscosity of air in physical units. */
        double nu_phys;
        /** The size of the turbulent fluid box in the x direction. */
        double lx;
        /** The size of the turbulent fluid box in the y direction. */
        double ly;
        /** The size of the turbulent fluid box in the z direction. */
        double lz;
        /** The reciprocal of the constant controlling the mode timescales in
         * the turbulent wind field. */
        double w_Cinv;
        /** The constant controlling the mode energy scales in the turbulent
         * wind field. */
        double w_alpha;
        /** The root mean squared (RMS) wind speed. */
        double w_rms;
        /** The padding factor in the maximum distance calculation for the GPR class. */
        double gpr_dist_pad;
        /** The tolerance factor for detecting duplicate locations in the GPR class. */
        double gpr_ker_tol;
        /** The timestep padding factor for the glider. */
        double gl_pad;
        /** The timestep padding factor for the turbulent fluid. */
        double tf_pad;
        /** The duration of a control step. */
        double c_dur;
        /** The simulation duration. */
        double duration;
        /** The time cutoff before which to take no measurements of wind
         * correlation. */
        double wic_tcut;
        /** The time cutoff before which to take no measurements of the
         * MCTS tree info. */
        double mti_tcut;
        /** The type of model used for predicting the wind field. */
        wind_model wmodel;
        /** The numerical integration type. */
        integration_type itype;
        /** The planning model type, used by the glider for deciding how to
         * bank. */
        planning_type ptype;
        /** A pointer to the class containing glider model parameters. */
        glider_model *gm;
        /** The output directory filename. */
        char *filename;
        /** The buffer for assembling output filenames. */
        char *fbuf;
        /** A vector of information about fluid cross-sections to output. */
        std::vector<int> w_cs;
        fileinfo(const char* infile);
        /** The class destructor frees the dynamically allocated memory. */
        ~fileinfo() {
            delete gm;
            delete [] fbuf;
            delete [] filename;
        }
        void select_timestep(bool verbose=false);
        void print_info(FILE *fp=stdout);
        /** Prints information about all of the constants stored within the
         * class.
         * \param[in] filename the name of the file to write to. */
        inline void print_info(const char* filename) {
            FILE *fp=safe_fopen(filename,"w");
            print_info(fp);
            fclose(fp);
        }
    protected:
        /** Returns whether or not an interpolation method is used for wind
         * prediction.
         * \return True if the linear or cubic interpolation method is used,
         * flase otherwise. */
        inline bool wm_interpolation() {
            return wmodel==wm_full_linear||wmodel==wm_full_cubic;
        }
        /** Returns whether the wind field is frozen.
         * \return True if frozen, false otherwise. */
        inline bool frozen() {
            return w_Cinv==0;
        }
        /** Returns whether or not the simulation will need to use the
         * mean-reversion velocity field prediction class.
         * \return True if the turb_fluid_grid_mr class is needed, false
         * otherwise. */
        inline bool mr_pred_used() {
            return !frozen()&&wm_interpolation();
        }
        /** The integration timestep. */
        double dt;
        /** The duration of an integration interval. */
        double i_dur;
        /** The duration between output snapshots. */
        double out_dur;
        /** The number of steps in an integration interval. */
        int steps;
        /** The total number of output snapshots (not including the zeroth
         * snapshot). */
        int snaps;
        /** The number of integration intervals. */
        int intervals;
        /** The number of integration intervals per control event. */
        int i_pc;
        /** The number of integration intervals between GPR samples. */
        int gs_freq;
        /** The number of integration intervals between saving output. */
        int out_freq;
    private:
        /** Finds the next token in a string and interprets it as a double
         * precision floating point number. If none is availble, it gives an
         * error message.
         * \param[in] ln the current line number. */
        inline double next_double(int ln) {
            return atof(next_token(ln));
        }
        /** Finds the next token in a string, interprets it as a double
         * precision floating point number, and checks that there are no
         * subsequent values.
         * \param[in] ln the current line number. */
        inline double final_double(int ln) {
            double temp=next_double(ln);
            check_no_more(ln);
            return temp;
        }
        /** Finds the next token in a string, interprets it as an integer, and
         * checks that there are no subsequent values.
         * \param[in] ln the current line number. */
        inline int final_int(int ln) {
            int temp=atoi(next_token(ln));
            check_no_more(ln);
            return temp;
        }
        /** Tests to see if two strings are equal.
         * \param[in] p1 a pointer to the first string.
         * \param[in] p2 a pointer to the second string.
         * \return True if they are equal, false otherwise. */
        inline bool se(const char *p1,const char *p2) {
            return strcmp(p1,p2)==0;
        }
        int gcd(int a,int b);
        /** Computes the least common multiple of two
         * integers.
         * \param[in] (a,b) the two integers.
         * \return The least common multiple. */
        inline int lcm(int a,int b) {
            return a/gcd(a,b)*b;
        }
        char* next_token(int ln);
        void check_no_more(int ln);
        void check_invalid(double val,const char *p);
        void print_range(FILE *fp,int rmin,int rmax);
        inline double max_tf_timestep();
        inline void calculate_wind_param(bool set_rms);
        inline double sqr(double x) {return x*x;}
        /** Returns a string describing the current wind model. */
        inline const char* s_wind_model() {
            switch(wmodel) {
                case wm_gpr: return "GPR";
                case wm_full_linear: return "Full info, linear";
                case wm_full_cubic: return "Full info, cubic";
                default: break;
            }
            return "Unset";
        }
        /** Returns a string describing the current integration type. */
        inline const char* s_integration_type() {
            switch(itype) {
                case it_euler: return "Euler";
                case it_improv_e: return "Improved Euler";
                default: break;
            }
            return "Unset";
        }
        /** Returns a string describing the current planning type. */
        inline const char* s_planning_type() {
            switch(ptype) {
                case pt_zero: return "Zero";
                case pt_random: return "Random";
                case pt_mcts: return "MCTS";
                default: break;
            }
            return "Unset";
        }
        /** Computes the scaling factor to apply to a mode, to take
         * into account symmetries.
         * \param[in] (i,j,k) the mode index to consider.
         * \return The scaling factor, which can be zero, one, or two. */
        inline int f_mode(int i,int j,int k) {
            return i!=0?1:(2*k>=nz?(2*k==nz?1:(2*j==ny?1:0))
                                  :(k>0?(2*j>=ny?(2*j==ny?1:2):2)
                                       :(2*j>=ny?(2*j==ny?1:0):(j==0?0:2))));
        }
};

#endif
