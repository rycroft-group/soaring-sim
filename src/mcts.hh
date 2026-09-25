#ifndef MCTS_HH
#define MCTS_HH

#include <gsl/gsl_rng.h>
#include <cmath>

#include "common.hh"
#include "node.hh"
#include "stats.hh"

/** The initial memory allocation for the tree. */
const int mcts_init_tree_mem=4096;

/** The maximum memory allocation for the tree. */
const int mcts_max_tree_mem=16777216;

class mcts {
    public:
        /** The maximum number of children per node. Equal to maximum number of
         * actions. */
        const short m;
        /** Flags controlling diagnostic output. 1: call diagnostic play
         * function, 2: collect MCTS tree info. */
        unsigned short dflags;
        /** The depth of the Monte Carlo tree search, corresponding to the total
         * number of moves to make in a path. */
        const int depth;
        /** The exploration factor, specified as a dimensionless constant. */
        const float ex_fac;
        /** The current tree size. */
        int s;
        /** The current tree memory. */
        int mem;
        /** The total number of simulations at the root node. */
        int root_n;
        /** The total number of valid actions at the root node. */
        short root_nact;
        /** An array for storing the counts of depths in the MCTS. */
        int* const ti;
        /** An array for storing MCTS depth statistics. */
        mti_stats* ms;
        /** The tree. */
        node* t;
        mcts(short m_,int depth_,float ex_fac_,bool tdiag,unsigned long seed=1);
        ~mcts();
        template<bool diag,class game,class state>
        void simulate_paths(game &ga,state &st,int id,int num);
        void create_leaf(short nact,short* v_act);
        int select_child(int d,int n,short nact);
        short pick_best();
        void print_tree();
    private:
        /** The GSL random number generator. */
        gsl_rng *rng;
        /** The weight array for selecting children. */
        float* const wei;
        /** The moves taken in a simulated playout. */
        short* const actions;
        void add_tree_memory();
};

#endif
