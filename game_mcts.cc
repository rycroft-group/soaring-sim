#include "game_mcts.hh"
#include "common.hh"

#include <cstring>

game_mcts::game_mcts(int m_,int depth_,game_glider &ttt_) : m(m_),depth(depth_), s(0), mem(mcts_init_tree_mem),
    rt(0), t(new node[mem]), rng(gsl_rng_alloc(gsl_rng_taus2)),
    wei(new float[m]), ttt(ttt_) {
    create_leaf();
}

/** The class destructor frees the dynamically allocated memory. */
game_mcts::~game_mcts() {
    delete [] wei;
    gsl_rng_free(rng);
    delete [] t;
}

/** A single simulation by playing till game over.
 * \param[in] p player id (1 or 2). */
void game_mcts::simulate_path(int p,int time) {
    int nt=rt;
    int moves[256],i=0;
    float win;

    int d=0;
    while(true) {

        // Choose a move, record it, and play the move on the board
     int k=select_child(d,nt,p);
        moves[i++]=k;
        if(ttt.play(k,p,time)==1) {
            win=p==1?1:0;break;
        }
        if(ttt.full()) {
            win=0.5;break;
        }
        p^=3;
        nt=t[d+k].n;
        // Check if the new node on the tree has a child or not
        if(t[d+k].c==0) {

            // There is no child. Create a new leaf here.
            t[d+k].c=s;
            create_leaf();

            // Pick a random move to make
            int l=i++;
            while(true) {
                do {
                    k=gsl_rng_uniform_int(rng,m);
                } while(!ttt.valid_move(k));
        moves[l++]=k;
                if(ttt.play(k,p,time)==1) {
                    win=p==1?1:0;break;
                }
                if(ttt.full()) {
                    win=0.5;break;
                }
                p^=3;
            }
            // for(int o=i;o<l;o++) ttt.board[moves[o]]=0;
            for(int o=l-1;o>i-1;o--) {
//          if(time==9999 || time==0) ttt.print();
            ttt.remove(moves[o]);
        }
            break;
        } else {

            // There is a child
            d=t[d+k].c;
        }
    }

    // Backpropagation
    d=0;
    for(int j=i-1;j>-1;j--) {
        // ttt.board[moves[j]]=0;
//  if(time==9999 || time==0) ttt.print();
        ttt.remove(moves[j]);

    }
    for(int j=0;j<i;j++){
    t[d+moves[j]].update(win);
        d=t[d+moves[j]].c;
    }
    rt++;
}

/** A single simulation for glider.
 * \param[in] p player id (1 or 2).
 * \param[in] time the current time frame */
void game_mcts::simulate_path_glider(int p,int time,bool output,FILE *fp) {
    int nt=rt;
    int moves[depth+1],i=0;
    double win=0;
    double win_step;

    double rx=ttt.rx;
    double ry=ttt.ry;
    double rz=ttt.rz;

    // Get initial glider velocity and wind velocity
    double ux0=ttt.ux,uy0=ttt.uy,uz0=ttt.uz;
    double wx0=ttt.wx,wy0=ttt.wy,wz0=ttt.wz;
    double vx0=ttt.vx,vy0=ttt.vy,vz0=ttt.vz;
    double bank_=ttt.bank;

    int d=0;
    while(i<depth){

   	 //printf("simulate step %d\n",i);
        // Choose a move, record it, and play the move on the board
        int k=select_child(d,nt,p);
	ttt.bank+=ttt.mu_f(k);
	double temp=0;
        for(int j=0;j<50;j++){
            moves[i++]=k;
	   // printf("simulate step %d\n",j);
           printf("SP %d %d %d\n",i,j,time);
            win_step=ttt.simulate_play(k,p,time+i);
            if(output) ttt.output(i+1,k,fp);
           if(std::isnan(win_step)){printf("nan win %g k %d i %d \n",win_step,k,i);exit(0);}
	       win+=win_step;
	       temp+=win_step;
        }
	    double vx1=ttt.vx,vy1=ttt.vy,vz1=ttt.vz;
            double wx1=ttt.wx,wy1=ttt.wy,wz1=ttt.wz;
            temp+=(vx0*wx0+vy0*wy0+vz0*wz0-vx1*wx1-vy1*wy1-vz1*wz1);
            t[d].update_step(temp/50.);
         // puts("pass simulate play");
       	nt=t[d+k].n;
        // Check if the new node on the tree has a child or not
        if(t[d+k].c==0) {

            // There is no child. Create a new leaf here.
            t[d+k].c=s;
            create_leaf();

            // Pick a random move to make
            int l=i++;
            while(l<depth) {
                do {
                    k=gsl_rng_uniform_int(rng,m);
                } while(!ttt.valid_move(k));
		          ttt.bank+=ttt.mu_f(k);
                for(int j=0;j<50;j++){
                   moves[l++]=k;
//		   printf("l %d\n",l);
                    win_step=ttt.simulate_play(k,p,time+l);
                    if(output) ttt.output(l+1,k,fp);
                    win+=win_step;
                  //  double vx1=ttt.vx,vy1=ttt.vy,vz1=ttt.vz;
                   // double wx1=ttt.wx,wy1=ttt.wy,wz1=ttt.wz;
                   // t[d+k].update_step(win_step+(vx0*wx0+vy0*wy0+vz0*wz0-vx1*wx1-vy1*wy1-vz1*wz1));
                }
            }
            break;
        } else {

            // There is a child
            d=t[d+k].c;
        }
    }

    // Get final glider velocity and wind velocity
    double vx1=ttt.vx,vy1=ttt.vy,vz1=ttt.vz;
    double wx1=ttt.wx,wy1=ttt.wy,wz1=ttt.wz;

    // Calculate the total energy gain in the simulation
    win+=(vx0*wx0+vy0*wy0+vz0*wz0-vx1*wx1-vy1*wy1-vz1*wz1);

    // Backpropagation
    d=0;
    double avg_win=win;
    // printf("final win %g time %d i %d\n",avg_win,time,i);
    for(int j=0;j<i-1;j++){
   	    t[d+moves[j]].update(avg_win);
        d=t[d+moves[j]].c;
    }
    rt++;

    // Reset to the start position
    ttt.reset(rx,ry,rz,ux0,uy0,uz0,vx0,vy0,vz0,bank_);
}

/** Select a child in tree based on weighting.
 * \param[in] d parent node index.
 * \param[in] nt number of simulations at the parent node. */
int game_mcts::select_child(int d,int nt,int p) {
    int k;
    float fac=sqrt(2*log(nt)),tw=0;

    // Calculate total weight of children
    // Check valid child (could not be parents )
    for(k=0;k<m;k++) {
        // w[k]=t[d].weight(fac);
        if(ttt.valid_move(k)){
            wei[k]=t[d+k].weight(fac,p==2);
            tw+=wei[k];
        }
        else wei[k]=0;
    }

    // Choose a weighted random child
    float r=gsl_rng_uniform(rng)*tw;
    k=0;
    while(k<m-1) {
        r-=wei[k];
        if(r<0) return k;
        k++;
    }
    return m-1;
}

void game_mcts::print_tree() {
    for(int i=0;i<s;i++) {
        printf("Node %d: n=%d, w=%g, c=%d\n",i,t[i].n,t[i].weight(0),t[i].c);
        if(i%9==8) putchar('\n');
    }
}

/** Find the max weight child.
 * \param[in] d index of parent node
 * \return the best child. */
int game_mcts::pick_best(int d) {
    int k,idx=0;
    float w,max_w=t[d].weight(0);

    //printf("k %d weight %g\n",k,max_w);
    for(k=1;k<m;k++) {
        w=t[d+k].weight(0);
        if(w>max_w) idx=k,max_w=w;
//	printf("k %d weight %g\n",k,w);
    }
    return idx;
}

void game_mcts::create_leaf() {
    if(s+m>mem) add_tree_memory();
    for(int k=0;k<m;k++) t[s+k].init();
    s+=m;
}

/** Doubles the memory in the tree array. */
void game_mcts::add_tree_memory() {

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
