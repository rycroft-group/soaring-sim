#include "common.hh"
#include "reinf_learn.hh"

#include <cstdio>
#include <cstdlib>
#include <algorithm>

/** Constructs the reinforcement learning class.
 * \param[in] ns_ the number of states.
 * \param[in] na_ the number of actions.
 * \param[in] alpha_ the learning rate.
 * \param[in] gamma_ the discount factor. */
reinf_learn::reinf_learn(int ns_,int na_,double alpha_,double gamma_)
    : ns(ns_), na(na_), nsa(ns*na), alpha(alpha_), gamma(gamma_),
    Q(new double[nsa]), Freq(new int[nsa]), ba(new int[na]),
    rng(gsl_rng_alloc(gsl_rng_taus2)) {
    reset_Q();
}

/** The class destructor frees the dynamically allocated memory, and frees the
 * random number generator. */
reinf_learn::~reinf_learn() {
    gsl_rng_free(rng);
    delete [] ba;
    delete [] Q;
    delete [] Freq;
}

/** Resets the value array.
 * \param[in] Qinit the initial value for each (state,array) pair. */
void reinf_learn::reset_Q(double Qinit) {

    // shuffle_Q();
    // straight_Q(1.1);
    // trained_Q();
    // load_Q();
    // straight_Q2();

    for(int i=0;i<nsa;i++) {
        // Q[i]=gsl_rng_uniform(rng);
        Q[i]=gsl_ran_gaussian(rng,1.0)+1;
        // Q[i]=0;
        Freq[i]=0;
    }
}

void reinf_learn::trained_Q(){

    for(int i=0;i<nsa;i++) {
        // Q[i]=gsl_rng_uniform(rng);
        Q[i]=0;
        Freq[i]=0;
    }

    FILE *fp=fopen("nt1.odr/Final_Q_values_111.bin","rb");
    if (fp==NULL){
        fputs("Error opening file\n",stderr);
        exit(1);
    }

    for(int i=0;i<nsa;i++) {fscanf(fp,"%lf\n",&Q[i]);}
    fclose(fp);

    double aa;
    fp=fopen("nt1.odr/Final_Q_values_222.bin","rb");
    if (fp==NULL){
        fputs("Error opening file\n",stderr);
        exit(1);
    }

    for(int i=0;i<nsa;i++) {fscanf(fp,"%lf\n",&aa); Q[i]+=aa;}
    fclose(fp);

    fp=fopen("nt1.odr/Final_Q_values_333.bin","rb");
    if (fp==NULL){
        fputs("Error opening file\n",stderr);
        exit(1);
    }

    for(int i=0;i<nsa;i++) {fscanf(fp,"%lf\n",&aa); Q[i]+=aa;}
    fclose(fp);

}

void reinf_learn::shuffle_Q(){
    trained_Q();
    std::random_shuffle(&Q[0],&Q[nsa-1]);
}

void reinf_learn::load_Q(){

    for(int i=0;i<nsa;i++) {
        Q[i]=0;
        Freq[i]=0;
    }

    FILE *fp=fopen("nt1.odr/Direct_Q_values.bin","rb");
    if (fp==NULL){
        fputs("Error opening file\n",stderr);
        exit(1);
    }

    for(int i=0;i<nsa;i++) {fscanf(fp,"%lf\n",&Q[i]);}
    fclose(fp);

}

void reinf_learn::straight_Q(double directness){

    for(int i=0;i<nsa;i++) {
        Q[i]=0;
        Freq[i]=0;
    }

    // assign direction for rows
    for(int i=0;i<20;i++){
        double pp=gsl_rng_uniform(rng);
        if(pp<directness){
            double prob=gsl_rng_uniform(rng);
            int a=prob>0.5?1:3;
            for(int j=0;j<20;j++){
                for(int k=0;k<5;k++){
                    Q[na*(k+5*(i+20*j))+a]=1;
                }
            }
        }
    }

    // assgin direction for columns
    for(int j=0;j<20;j++){
        double pp=gsl_rng_uniform(rng);
        if(pp<directness){
            double prob=gsl_rng_uniform(rng);
            int a=prob>0.5?0:2;
            for(int i=0;i<20;i++){
                for(int k=0;k<5;k++){
                    Q[na*(k+5*(i+20*j))+a]=1;
                }
            }
        }
    }

}

void reinf_learn::straight_Q2(){

    for(int i=0;i<nsa;i++) {
        Q[i]=0;
        Freq[i]=0;
    }

    int rol_col_idx[40];
    for(int i=0;i<40;i++){
        rol_col_idx[i]=i;
    }
    std::random_shuffle(&rol_col_idx[0],&rol_col_idx[39]);

    for(int i=0;i<40;i++){
        printf("random shuffle index %3d\n",rol_col_idx[i]);
    }

    for(int t=0;t<40;t++){
        int idx=rol_col_idx[t];
        if(idx<20){
            double prob=gsl_rng_uniform(rng);
            int a=prob>0.5?1:3;
            for(int j=0;j<20;j++){
                for(int k=0;k<5;k++){
                    for(int act=0;act<5;act++ ){
                        Q[na*(k+5*(idx+20*j))+act]=0;
                    }
                    Q[na*(k+5*(idx+20*j))+a]=1;
                }
            }

        }

        else{

            double prob=gsl_rng_uniform(rng);
            int a=prob>0.5?0:2;
            for(int i=0;i<20;i++){
                for(int k=0;k<5;k++){
                    for(int act=0;act<5;act++ ){
                        Q[na*(k+5*(i+20*(idx-20)))+act]=0;
                    }
                    Q[na*(k+5*(i+20*(idx-20)))+a]=1;
                }
            }
        }
    }

}

/** Calculates the maximum value over all actions for a given state.
 * \param[in] s the state.
 * \return The maximum value. */
double reinf_learn::max_Q(int s) {
    double mQ=Q[na*s];
    for(int a=1;a<na;a++) if(Q[na*s+a]>mQ) mQ=Q[na*s+a];
    return mQ;
}

double reinf_learn::min_Q(int s) {
    double mQ=Q[na*s];
    for(int a=1;a<na;a++) if(Q[na*s+a]<mQ) mQ=Q[na*s+a];
    return mQ;
}

/** Prints the Q values. */
void reinf_learn::print_Q() {
    for(int s=0;s<ns;s++) {
        printf("# s=%3d:",s);
        for(int a=0;a<na-1;a++) printf(" %8g",Q[na*s+a]);
        printf(" %8g\n",Q[na*s+(na-1)]);

    }
    /*
    puts("Frequency");
    for(int s=0;s<ns;s++) {

        printf("# s=%3d:",s);
        for(int a=0;a<na-1;a++) printf(" %8d",Freq[na*s+a]);
        printf(" %8d\n",Freq[na*s+(na-1)]);
    }*/

}

/** Chooses a set of actions An={a1,a2,...,an} for given starting state,with a probability epsilon
 * of picking a random action instead.
 * \param[in] s the state.
 * \param[in] n_step number of steps to explore.
 * \return the chosen action set. */
int reinf_learn::pick_best_uct(int s,double eps,int n_step){

}

void reinf_learn::update_uct()
/** Chooses the best action for a given state, with a probability epsilon
 * of picking a random action instead.
 * \param[in] s the state.
 * \return The chosen action. */
int reinf_learn::pick_best_with_random(int s,double eps) {
    double b=gsl_rng_uniform(rng);
    Freq[s]++;
    return b<eps?static_cast<int>(na*b/eps)
                :pick_best(s);
}

/** Chooses the best action for a given state.
 * \param[in] s the state.
 * \return The chosen action. */
int reinf_learn::pick_best(int s) {

    // Find all of the actions that are within a small tolerance of the best
    // action
    double mQ=max_Q(s);
    mQ*=mQ>0?1-reinf_learn_equiv_factor
            :1+reinf_learn_equiv_factor;
    int k=0;
    for(int a=0;a<na;a++) if(Q[na*s+a]>=mQ) ba[k++]=a;

    // Check for the case when no action is found, and give an error
    if(k==0) fatal_error("Error finding best action",1);

    // Randomly choose from the best actions
    return ba[k==1?0:gsl_rng_uniform_int(rng,k)];
}

/** Chooses the action with probability calculated from Q[s,a] for a given state.
 * \param[in] s the state.
 * \return The chosen action. */
int reinf_learn::pick_best_with_prob(int s,double eps) {

    double minQ=min_Q(s);
    double maxQ=max_Q(s);

    // Check if all Q[s,a]=0 for given state s
    if (minQ==0 and maxQ==0) return pick_best_with_random(s,eps);

    // Calculate the probability according to Q[s,a]
    double p [5];
    for(int a=0;a<na;a++) p[a]=Q[na*s+a]+abs(minQ)+0.0001;
    double sum_p=0;
    for(int a=0;a<na;a++) sum_p+=p[a];
    for(int a=0;a<na;a++) p[a]=p[a]/sum_p;

    // Check which action to choose
    double accumulate_p [5];
    accumulate_p[0]=0;
    for(int a=0;a<na-1;a++) accumulate_p[a+1]=accumulate_p[a]+p[a];
    double b=gsl_rng_uniform(rng);
    double ma=0;
    for(int a=0;a<na;a++) if(b>accumulate_p[a]) ma=a;

    return ma;

}

/** Writes the Q values. */
void reinf_learn::write_Q(const char* filename) {

    FILE *fp=fopen(filename,"w");
    if (fp==NULL){
        fputs("Error opening file\n",stderr);
        exit(1);
    }
    for(int i=0;i<nsa;i++) fprintf(fp,"%.10f\n",Q[i]);
    fclose(fp);

}

/** Writes the frequency values. */
void reinf_learn::write_Freq(const char* filename) {

    FILE *fp=fopen(filename,"w");
    if (fp==NULL){
        fputs("Error opening file\n",stderr);
        exit(1);
    }
    for(int i=0;i<nsa;i++) fprintf(fp,"%.10d\n",Freq[i]);
    fclose(fp);

}
