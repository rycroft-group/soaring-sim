#ifndef SOARING_SIM_HH
#define SOARING_SIM_HH

#include <cstdio>
#include <gsl/gsl_rng.h>

#include "fileinfo.hh"
#include "glider.hh"
#include "tf_grid.hh"
#include "tf_grid_mr.hh"
#include "k_func.hh"
#include "gpr.hh"
#include "mcts.hh"
#include "wind_correl.hh"
#include "path_store.hh"
#include "stats.hh"

class soaring_sim : public fileinfo {
    public:
        /** The number of OpenMP threads. */
        const int nt;
        /** The number of OpenMP threads to use for the Monte Carlo tree
         * search. It is equal to the minimum of the total threads available,
         * and the total number of gliders. */
        const int mctst;
        /** The current fluid simulation number, used to set snpashot output
         * filenames. */
        int fsim;
        /** The current time. */
        double time;
        /** A boolean value controlling whether to output the MCTS paths. */
        bool path_output;
        /** A pointer to the class with the turbulent wind field. */
        turb_fluid_grid* const tf;
        /** The array of glider states. */
        glider* const g;
        /** A second array of glider states, needed when using the improved
         * Euler method for numerical integration. */
        glider* const g_ie;
        /** A pointer to the kernel function used with Gaussian process
         * regression (GPR). */
        kernel_func* KF;
        /** An array of pointers to the Gaussian process regression (GPR) classes
         * for predicting the wind field. */
        gpr** const mem;
        /** An array of pointers to the Monte Carlo tree search (MCTS) classes
         * for planning the glider banking. */
        mcts** const mc;
        soaring_sim(const char* fileinfo);
        ~soaring_sim();
        void run();
        void init_gliders();
        void simulate_gliders();
        void plan();
        void plan_mcts();
        inline void play(glider &g_,int id,int step,short action,bool diag) {
            diag?play_internal<true>(g_,id,step,action)
                :play_internal<false>(g_,id,step,action);
        }
        void output_start();
        void snapshots_start();
        void write_snapshots(int s);
        void snapshots_end();
        void output_end();
        float base_score(glider &g_,int id);
        float score(glider &g_,int id,float bs);
        /** Calculates the valid actions for a given glider state.
         * \param[in] g_ the current glider state.
         * \param[out] v_act a pointer to an array containing the valid
         *                   actions, as short integers.
         * \return The total number of valid actions. */
        inline short valid_actions(glider &g_,short*& v_act) {
            return gm->valid_actions(g_,v_act);
        }
    private:
        /** Computes the wind velocities at the current glider positions, if
         * they aren't already up to date. */
        inline void glider_velocities() {
            if(flag_current&1) return;
            for(int j=0;j<gpt;j++) g[j].copy_pos(pos+3*j);
            tf->vel_multi(gpt,pos,wnd);
            flag_current|=1;
        }
        /** Computes the grid-based representation of the turbulent wind field,
         * if it isn't already up to date. */
        inline void grid_compute() {
            if(flag_current&2) return;
            tf->transform();
            flag_current|=2;
        }
        /** Adds a measurement to the GPR classes for each glider at the
         * current time and position.
         * \param[in] time the current simulation time. */
        inline void gpr_measurement(double time) {
            glider_velocities();
#pragma omp parallel for
            for(int j=0;j<gpt;j++) mem[j]->update_measurement(pos+3*j,time,wnd+3*j);
        }
        /** Predicts the wind field at a given glider position, using the
         * currently selected prediction model.
         * \param[in] t_ the time into the future for prediction.
         * \param[in] id the ID of this glider.
         * \param[in] g_ the glider state.
         * \param[in] w a pointer to array for writing the predicted wind
         *              components. */
        inline void predict_wind(double t,int id,glider &g_,double *w) {
            switch(wmodel) {
                case wm_gpr: mem[id]->predict(g_.rx,g_.ry,g_.rz,time+t,*w,w[1],w[2]);break;
                case wm_full_linear:
                    bypass_mr||frozen()?tf->lin_interp(g_.rx,g_.ry,g_.rz,*w,w[1],w[2])
                          :((turb_fluid_grid_mr*)tf)->lin_interp_mr(t,g_.rx,g_.ry,g_.rz,*w,w[1],w[2]);
                    break;
                case wm_full_cubic:
                    bypass_mr||frozen()?tf->cub_interp(g_.rx,g_.ry,g_.rz,*w,w[1],w[2])
                          :((turb_fluid_grid_mr*)tf)->cub_interp_mr(t,g_.rx,g_.ry,g_.rz,*w,w[1],w[2]);
                    break;
                default: break;
            }
        }
        /** Runs an MCTS planning step, enabling the output of the MCTS paths.
         */
        inline void plan_with_paths() {
            path_output=true;
            plan();
            printf(" (Path)");
            for(path_info **p=pth;p<pth+gpt;p++) (*p)->output_and_reset(path_file);
            path_output=false;
        }
        template<bool diag>
        void play_internal(glider &g_,int id,int step,short action);
        /** A flag indicating whether the glider velocity array is up to date
         * (in the 1 bit) and whether the grid-based wind fieled is up to date
         * (in the 2 bit). */
        unsigned int flag_current;
        /** An array for holding the current glider positions, used to compute
         * the wind velocities. */
        double* const pos;
        /** An array holding the wind velocities. */
        double* const wnd;
        /** The wind correlation data structure. */
        wind_correl *wic;
        /** The path output data structure. */
        path_info **pth;
        /** A pointer to a GSL random number generator. */
        gsl_rng *rng;
        /** A file handle for writing glider position information. */
        FILE *glider_xyz;
        /** A file handle for writing glider position information in a minimal
         * binary format. */
        FILE *glider_mb;
        /** A file handle for writing glider position information in a larger
         * binary format (with wind data). */
        FILE *glider_lb;
        /** A file handle for writing the MCTS paths in binary format. */
        FILE *path_file;
        /** A file handle for writing the GPR checksums. */
        FILE *checksums_file;
        /** Initial glider heights. */
        double* ginit;
        /** The climb statistics. */
        cli_stats* cstats;
        /** The energy component statistics. */
        cli_stats* enstats;
};

#endif
