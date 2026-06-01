#ifndef GAME_MCTS_HH
#define GAME_MCTS_HH

/** The initial memory allocation for the tree. */
const int mcts_init_tree_mem=4096;

/** The maximum memory allocation for the tree. */
const int mcts_max_tree_mem=16777216;

#include <gsl/gsl_rng.h>
#include "game.hh"
#include <cmath>

struct node {
    /** The number of visits. */
    int n;
    /** The child node. */
    int c;
    /** The expected win in one trajactory. */
    float w;
     /** The expected increased/ decreased win in one trajactory. */
    float dw;
    /** The expected win in one step. */
    float w_step;
    inline void init() {
        n=0;c=0;w=0.;dw=0.;
    }
    inline float weight(float fac,bool opponent=false) {
//	printf("n %d w %g dw %g \n",n,w/n,dw);
        return n>0?(opponent?n-w:w)/n+fac/sqrt(n)+exp(0.05/(abs(dw-w/n)+0.05)):0.5;
    }
    inline void update(float w_) {
        n++;
        w+=w_;
        dw=w_;
    }
    inline void update_step(float w_){
        w_step+=w_;
    }
};

class game_mcts {
    public:
        /** The number of children per node. Equal to number of actions. */
        const int m;
	/** The searching depth of the tree. */
        const int depth;
        /** The current tree size. */
        int s;
        /** The current tree memory. */
        int mem;
        /** The total number of simulations at the root node. */
        int rt;
        /** The tree. */
        node* t;
        game_mcts(int m_,int depth_,game_glider &ttt_);
        ~game_mcts();
        void simulate_path(int p,int time);
        void simulate_path_glider(int p,int time,bool output,FILE *fp);
        void create_leaf();
        int select_child(int d,int nt,int p);
        int pick_best(int d);
        void print_tree();
        inline void reset() {
            rt=s=0;create_leaf();
        }
    private:
        /** The GSL random number generator. */
        gsl_rng *rng;
        /** The weight array for selecting children. */
        float *wei;
        /** A reference to a tic-tac-toe game. */
        game_glider &ttt;
        void add_tree_memory();
};

#endif
