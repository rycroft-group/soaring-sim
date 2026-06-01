#include <cstring>
#include <cstdio>

#include <gsl/gsl_randist.h>

#include "soaring_sim.hh"
#include "tf_grid_mr.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

/** Initializes the soaring simulation by reading parameters from a
 * configuration file, and allocating memory for the gliders, prediction
 * models, and turbulent flow simulation.
 * \param[in] filename the name of the configuration file to read from. */
soaring_sim::soaring_sim(const char* filename) : fileinfo(filename),
    nt(omp_get_max_threads()), mctst(gpt<nt?gpt:nt), path_output(false),
    tf(mr_pred_used()?(turb_fluid_grid*)
       new turb_fluid_grid_mr(nx,ny,nz,0,lx,0,ly,0,lz,w_Cinv,
                              w_alpha,h_segs,c_dur*mcts_depth,base_seed+1)
      :new turb_fluid_grid(nx,ny,nz,0,lx,0,ly,0,lz,w_Cinv,w_alpha,base_seed+1)),
    g(new glider[gpt]), g_ie(itype==it_improv_e?new glider[2*gpt]:NULL),
    mem(wmodel==wm_gpr?new gpr*[gpt]:NULL),
    mc(ptype==pt_mcts?new mcts*[mctst]:NULL),
    pos(new double[6*gpt]), wnd(pos+3*gpt),
    wic(fflags&256?new wind_correl(gpt,n_mcts,tf):NULL),
    pth(fflags&1024?new path_info*[gpt]:NULL),
    rng(gsl_rng_alloc(gsl_rng_taus2)) {
    gsl_rng_set(rng,base_seed);

    // If MCTS is in use, allocate memory for the MCTS calculations
    if(ptype==pt_mcts) {
        unsigned long sbase=1+omp_get_max_threads()+base_seed;
#pragma omp parallel num_threads(mctst)
        {

            // Initialize the MCTS class using a default score value based on
            // steady-state falling in zero wind
            int t=omp_get_thread_num();
            mc[t]=new mcts(gm->rtot,mcts_depth,mcts_ex_fac,fflags&512,sbase+t);

            // Initialize path data storage if needed
            if(fflags&1024) {
                int n=(1+mcts_depth*out_pc)*n_mcts;
#pragma omp for
                for(int j=0;j<gpt;j++) pth[j]=new path_info(n);
            }
        }
    }

    // If the GPR calculations are in use, allocate memory for them
    if(wmodel==wm_gpr) {

        // Compute the maximum time and distance that the kernel function will
        // need to be evaluated for, to set up the bilinear interpolation table
        if(ptype!=pt_mcts) fatal_error("Need MCTS enabled for GPR model",1);
        double spd=gm->steady_speed(),
               tmax=c_dur*(mcts_depth+double(gprm)/gs_pc),
               rmax=tmax*gpr_dist_pad*spd,
               ker_step=spd*c_dur/gs_pc;

        // Allocate the kernel table. Add a small padding to the maximum time,
        // to ensure rounding errors don't cause a table miss when t=tmax.
        const double pad=(1+16*std::numeric_limits<double>::epsilon());
        KF=frozen()?(kernel_func*) new kernel_r(nx,lx,gpr_rs,rmax)
                   :(kernel_func*) new kernel_rt(nx,lx,w_Cinv,gpr_rs,gpr_ts,rmax,tmax*pad);
#pragma omp parallel for
        for(int j=0;j<gpt;j++) mem[j]=new gpr(gprm,*KF,gpr_ker_tol,ker_step,gpr_full_compute);
    }

    // Allocate memory in the turbulent fluid class for performing simultaneous
    // velocity calculations for all of the gliders
    tf->allocate_vel_table(gpt,fflags&1);
}

/** The class destructor frees the dynamically allocated memory. */
soaring_sim::~soaring_sim() {

    // Deallocate the GPR instances if needed
    if(wmodel==wm_gpr) {
        for(int j=gpt-1;j>=0;j--) delete mem[j];
        delete KF;
        delete [] mem;
    }

    // Deallocate the MCTS instances if needed
    if(ptype==pt_mcts) {
        for(int j=0;j<mctst;j++) delete mc[j];
        delete [] mc;
    }

    // Deallocate random number generator
    gsl_rng_free(rng);

    // Delete path info data storage
    if(pth!=NULL) {
        for(int j=gpt-1;j>=0;j--) delete pth[j];
        delete [] pth;
    }

    // Delete the wind correlation computation class
    if(wic!=NULL) delete wic;

    // Delete glider and wind arrays
    delete [] pos;
    if(g_ie!=NULL) delete [] g_ie;
    delete [] g;

    // Deallocate the turbulent wind field
    if(mr_pred_used()) delete (turb_fluid_grid_mr*) tf;
    else delete tf;
}

/** Runs the simulation, creating a number of instances of the turbulent
 * fluid, and modeling the gliders moving through it. */
void soaring_sim::run() {

    // If the wind correlation computation is enabled, then set the substep
    // markers
    if(wic!=NULL) wic->set_substep_markers(steps*i_pc);

    // Do any required initial setup for output
    output_start();

    // Loop over the different turbulent fluid simulations
    for(fsim=0;fsim<num_trials;fsim++) {

        // Initialize the gliders at random positions
        init_gliders();

        // Set up the turbulent fluid class in steady state
        tf->init_steady_state();
        flag_current=0;

        // Simulate the gliders in the turbulent fluid
        snapshots_start();
        simulate_gliders();
        snapshots_end();
    }

    // Output any summary statistics, and clean up any temporary memory used
    // for output
    output_end();
}

/** Initializes the gliders at a random positions, and with random directions.
 */
void soaring_sim::init_gliders() {
    double uh0=gm->uh0,uz0=gm->uz0;
    for(int i=0;i<gpt;i++) {
        double vx,vy;
        gsl_ran_dir_2d(rng,&vx,&vy);
        g[i].init(lx*gsl_rng_uniform(rng),
                  ly*gsl_rng_uniform(rng),
                  lz*gsl_rng_uniform(rng),
                  vx*uh0,vy*uh0,uz0,gsl_rng_uniform_int(rng,gm->rtot)+gm->br_min);
    }
}

/** Simulates the gliders moving in the turbulent wind field and using the
 * prediction model to control their bank angle to maximize energy gain. */
void soaring_sim::simulate_gliders() {

    // Make an initial measurement of the wind at each glider's location
    // if GPR is being used
    if(wmodel==wm_gpr) {
        for(int j=0;j<gpt;j++) mem[j]->reset();
        gpr_measurement(0);
    }
    time=0;
    plan();
    write_snapshots(0);

    for(int i=1;i<=intervals;i++) {
        for(int j=0;j<steps;j++) {

            // Compute the wind velocities at the glider positions, and use them to
            // integrate the gliders forward. For the Euler method, this will
            // complete the integration; for the improved Euler method, this will
            // be a preliminary step.
            glider_velocities();
            itype==it_improv_e?gm->ie_step1(gpt,dt,g,g_ie,wnd)
                              :gm->euler(gpt,dt,g,wnd);

            // Perform a stochastic integration step on the turbulent fluid modes
            if(!frozen()) tf->step_forward(dt);
            flag_current=0;

            // For the improved Euler method, compute the wind velocities with the
            // preliminary step applied, and use those to complete the integration
            // step
            if(itype==it_improv_e) {
                for(int k=0;k<gpt;k++) g_ie[2*k].copy_pos(pos+3*k);
                tf->vel_multi(gpt,pos,wnd);
                gm->ie_step2(gpt,dt,g,g_ie,wnd);
            }
        }

        // If GPR is being used, make a measurement of the wind at each
        // glider's location
        printf("# Interval %d",i);
        time=i_dur*i;
        if(wmodel==wm_gpr&&(i%gs_freq==0)) {printf(" (Measure)");gpr_measurement(time);}

        // Perform a planning calculation to update the gliders' bank angles
        if(i%i_pc==0) {
            printf(" (Plan)");
            fflags&1024&&i%(i_pc*path_interval)==0?plan_with_paths():plan();
        }

        // Write any request snapshot files
        if(i%out_freq==0) {puts(" (Write)");write_snapshots(i/out_freq);}
        else putchar('\n');
    }
}

/** Updates the gliders' banking angles based on planning using MCTS.
 * \param[in] path_output whether to output the paths (if enabled). */
void soaring_sim::plan_mcts() {

    // Determine whether diagnostic routines for wind correlation or MCTS tree
    // info need to be considered
    unsigned short dflags=fflags&512&&time>=mti_tcut?2:0;
    if(fflags&256&&time>=wic_tcut) {
        dflags|=1;
        wic->reset_counters();
    }
    if(path_output) dflags|=1;

#pragma omp parallel num_threads(mctst)
    {
        // Each thread has a unique MCTS instance to use, which can be reset
        // and reused to plan for multiple gliders
        mcts *mcp=mc[omp_get_thread_num()];
        mcp->dflags=dflags;
#pragma omp for
        for(int j=0;j<gpt;j++) {

            // Run the MCTS by simulating many paths to build a decision tree
            // weighted toward favorable paths
            dflags?mcp->simulate_paths<true>(*this,g[j],j,n_mcts)
                  :mcp->simulate_paths<false>(*this,g[j],j,n_mcts);

            // Update the glider's banking angle based on the best move from
            // the MCTS
            g[j].bank=mcp->pick_best();
        }
    }

    // If the wind correlation computation is enabled, then compute the wind at
    // all of the glider locations in MCTS, and add their contributions to the
    // correlation functions
    if(fflags&256&&time>=wic_tcut) wic->wind_computation();
}

/** Takes a given glider position and integrates it forward over one
 * control duration with given move.
 * \param[in,out] g_ the glider to update.
 * \param[in] id the ID number of the glider.
 * \param[in] step the step number of the play, used to calculate the current
 *                 time for predicting the wind field.
 * \param[in] action the action to apply to the glider. */
template<bool diag>
void soaring_sim::play_internal(glider &g_,int id,int step,short action) {
    double w[3],t=step*c_dur;
    int sof=steps*out_freq;
    g_.bank=action;
    if(itype==it_euler) {

        // Integrate by forward Euler
        for(int j=0;j<steps*i_pc;j++) {
            predict_wind(t,id,g_,w);

            // If needed, perform diagnostic computations for wind correlations
            // and path output
            if(diag) {
                if(wic!=NULL) wic->measure(j,step,id,g_,w[2]);
                if(path_output&&j%sof==0) pth[id]->store(g_,w);
            }
            gm->euler(dt,g_,w);
            t+=dt;
        }
    } else {
        glider g2,g3;

        // Integrate by improved Euler
        for(int j=0;j<steps*i_pc;j++) {
            predict_wind(t,id,g_,w);

            // If needed, perform diagnostic computations for wind correlations
            // and path output
            if(diag) {
                if(wic!=NULL) wic->measure(j,step,id,g_,w[2]);
                if(path_output&&j%sof==0) pth[id]->store(g_,w);
            }
            gm->ie_step1(dt,g_,g2,g3,w);
            t+=dt;
            predict_wind(t,id,g3,w);
            gm->ie_step2(dt,g_,g2,g3,w);
        }
    }

    // If path output is enabled, and this is the last step of a path, then do
    // an extra dummy wind prediction for the endpoint of the path
    if(diag&&path_output&&step==mcts_depth-1) {
        predict_wind(t,id,g_,w);
        pth[id]->store(g_,w);
    }
}

/** Computes a baseline score for the glider at the start of an MCTS playout.
 * This is later used to compute the delta change in energy.
 * \param[in] g_ the glider state.
 * \param[in] id the glider ID number.
 * \return The baseline score. */
float soaring_sim::base_score(glider &g_,int id) {
    double w[3];
    predict_wind(0,id,g_,w);
    return static_cast<float>(g_.energy(w));
}

/** Computes the score for an MCTS playout.
 * \param[in] g_ the glider state at the end of the playout
 * \param[in] id the glider ID number.
 * \param[in] bs the basline score at the start of the playout.
 * \return The score, computed at the change in glider energy from the start
 * until now. */
float soaring_sim::score(glider &g_,int id,float bs) {
    double w[3];
    predict_wind(c_dur*mcts_depth,id,g_,w);
    return static_cast<float>(g_.energy(w))-bs;
}

/** Updates the gliders' banking angles based on the selected planning
 * strategy.
 * \param[in] path_output whether to output the MCTS paths (if enabled). */
void soaring_sim::plan() {
    if(ptype==pt_mcts) {
        if(wm_interpolation()) {
            grid_compute();
            if(!frozen()) ((turb_fluid_grid_mr*)tf)->calc_mr_fields();
        }
        plan_mcts();
    } else if(ptype==pt_random) for(int j=0;j<gpt;j++) gm->random_move(g[j],rng);
}

// Explicit instantiation
#include "mcts_sp.cc"
template void mcts::simulate_paths<true,soaring_sim,glider>(soaring_sim&,glider&,int,int);
template void mcts::simulate_paths<false,soaring_sim,glider>(soaring_sim&,glider&,int,int);
template void soaring_sim::play_internal<true>(glider&,int,int,short);
template void soaring_sim::play_internal<false>(glider&,int,int,short);
