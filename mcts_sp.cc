#include "mcts.hh"

/** Performs a Monte Carlo tree search (MCTS) by simulating a number of
 * playouts, taking random actions that are biased toward more promising
 * outcomes.
 * \param[in] game a reference to game class for simulating the game state, and
 *                 evaluating scores and valid moves.
 * \param[in] state a reference to a class containing the current, complete
 *                  game state.
 * \param[in] id the numerical ID of this particular game state.
 * \param[in] num the total number of playouts to simulate. */
template<bool diag,class game,class state>
void mcts::simulate_paths(game &ga,state &st,int id,int num) {
    int i,k,d,n,o;
    short nact,*v_act;
    float bs=ga.base_score(st,id),win;
    bool pdiag=diag&&dflags&1,
         tdiag=diag&&dflags&2;
    root_n=0;
    s=0;

    // Clear the tree info array
    if(tdiag) for(o=0;o<=depth;o++) ti[o]=0;

    // Compute the number of valid actions at the tree root, and create a leaf
    // with space for branches ot all of these actions
    root_nact=ga.valid_actions(st,v_act);
    create_leaf(root_nact,v_act);

    for(o=0;o<num;o++) {

        // Make a copy of the current state, to perform a simulated playout.
        // Initialize constants for traversing the tree.
        state st_(st);
        n=root_n;
        nact=root_nact;
        d=0;
        i=0;
        while(i<depth) {

            // Choose a move, record it, and play the move on the board
            actions[i]=(k=select_child(d,n,nact));
            ga.play(st_,id,i++,t[d+k].act,pdiag);

            // Check if the new node on the tree has a child or not
            n=t[d+k].n;
            if(t[d+k].c==0) {
                if(i==depth) break;

                // There is no child. Create a new leaf here.
                t[d+k].c=s;
                t[d+k].nact=(nact=ga.valid_actions(st_,v_act));
                create_leaf(nact,v_act);

                // Take a random move from the new leaf
                k=gsl_rng_uniform_int(rng,nact);
                actions[i]=k;
                ga.play(st_,id,i++,v_act[k],pdiag);

                // Simulate a random playout
                for(int l=i;l<depth;l++) {
                    nact=ga.valid_actions(st_,v_act);
                    k=gsl_rng_uniform_int(rng,nact);
                    ga.play(st_,id,l,v_act[k],pdiag);
                }
                break;
            } else {

                // There is a child
                nact=t[d+k].nact;
                d=t[d+k].c;
            }
        }

        // Store the length of the tree traverse
        if(tdiag) ti[i]++;

        // Backpropagation
        d=0;
        win=ga.score(st_,id,bs);
        for(int j=0;j<i;j++) {
            t[d+actions[j]].update(win);
            d=t[d+actions[j]].c;
        }
        root_n++;
    }

    // Store the MCTS tree info if needed
    if(tdiag) for(o=0;o<=depth;o++) ms[o].contrib(ti[o]);
}
