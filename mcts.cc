#include "mcts.hh"

#include <cstring>
#include <limits>

/** Initializes the MCTS class, setting up the memory for the action tree.
 * \param[in] m_ the maximum number of actions.
 * \param[in] depth_ the maximum depth of the search tree.
 * \param[in] seed a seed for the GSL random number generator.
 * \param[in] ex_fac_ an exploration factor, specified as a dimensionless
 *                    number.
 * \param[in] tdiag whether to set up memory for the MCTS tree info. */
mcts::mcts(short m_,int depth_,float ex_fac_,bool tdiag,unsigned long seed)
    : m(m_), depth(depth_), ex_fac(ex_fac_), mem(mcts_init_tree_mem),
    ti(tdiag?new int[depth+1]:NULL), ms(tdiag?new mti_stats[depth+1]:NULL),
    t(new node[mem]), rng(gsl_rng_alloc(gsl_rng_taus2)), wei(new float[m]),
    actions(new short[depth]) {
    gsl_rng_set(rng,seed);

    // Clear the MCTS tree info statistics classes
    if(tdiag) for(int o=0;o<=depth;o++) ms[o].init();
}

/** The class destructor frees the dynamically allocated memory. */
mcts::~mcts() {
    delete [] actions;
    delete [] wei;
    gsl_rng_free(rng);
    delete [] t;

    // Free the MCTS tree info arrays if they were allocated
    if(ms!=NULL) {
        delete [] ms;
        delete [] ti;
    }
}

/** Select a child in tree based on weighting.
 * \param[in] d parent node index.
 * \param[in] n number of simulations at the parent node.
 * \param[in] nact the number of valid actions. */
int mcts::select_child(int d,int n,short nact) {
    short k,nz=0;
    float tw=0,wei_min=std::numeric_limits<float>::max(),
               wei_max=std::numeric_limits<float>::min();

    for(k=0;k<nact;k++) {
        if(t[d+k].n>0) {
            wei[k]=t[d+k].weight();
            if(wei[k]<wei_min) wei_min=wei[k];
            if(wei[k]>wei_max) wei_max=wei[k];
            tw+=wei[k];
        } else nz++;
    }

    // Choose a weighted random child
    short nc=nact-nz;
    if(nc<=1) return gsl_rng_uniform_int(rng,nact);
    float fac=wei_max-wei_min;
    if(fac==0) return gsl_rng_uniform_int(rng,nact);
    fac=1./fac;

    float extra=ex_fac/sqrt(n),ex1=-fac*wei_min+extra,ex2=0.5+extra,
          r=gsl_rng_uniform(rng)*(fac*tw+nc*ex1+nz*ex2);
    k=0;
    while(k<nact-1) {
        r-=t[d+k].n>0?fac*wei[k]+ex1:0.5+ex2;
        if(r<0) return k;
        k++;
    }
    return nact-1;
}

/** Prints information about all of the nodes in the tree. */
void mcts::print_tree() {
    printf("Root n: %d\nRoot nact: %hd\n",root_n,root_nact);
    for(int i=0;i<s;i++) {
        if(t[i].n==0) printf("Node %d: n=0, act=%hd\n",i,t[i].act);
        else{
            printf("Node %d: n=%d, c=%d, w=%g, act=%hd",
                    i,t[i].n,t[i].c,t[i].weight(),t[i].act);
            if(t[i].c>0) printf(", nact=%hd\n",t[i].nact);
            else putchar('\n');
        }
    }
}

/** Find the best action at the root of the tree to play.
 * \return The action. */
short mcts::pick_best() {
    short k0=0;
    while(true) {
        if(t[k0].n>0) break;
        k0++;
        if(k0>=root_nact) return t[gsl_rng_uniform_int(rng,root_nact)].act;
    }
    float max_wei=t[k0].weight(),w;
    for(short k=k0+1;k<root_nact;k++) {
        if(t[k].n>0) {
            w=t[k].weight();
            if(w>max_wei) k0=k,max_wei=w;
        }
    }
    return t[k0].act;
}

/** Allocates memory for a leaf, including space for new branches for each
 * valid action.
 * \param[in] nact the number of valid actions at this leaf. */
void mcts::create_leaf(short nact,short* v_act) {
    if(s+nact>mem) add_tree_memory();
    for(short k=0;k<nact;k++) t[s+k].init(v_act[k]);
    s+=nact;
}

/** Doubles the memory in the tree array. */
void mcts::add_tree_memory() {

    // Check if the absolute threshold on tree memory allocation has been
    // exceeded
    if(mem>=mcts_max_tree_mem)
        fatal_error("Maximum tree memory allocation reached",1);

    // Allocate a new tree array with double the size
    int nmem=mem<<1;
    node *nt=new node[nmem];

    // Copy the contents of the existing array into this array
    memcpy(nt,t,mem*sizeof(node));

    // Delete the old array and update the pointer to the new array
    delete [] t;
    t=nt;mem=nmem;
}
