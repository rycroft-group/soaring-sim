#include "wind_correl.hh"

#include <cstdlib>

/** Outputs the internal variables and Pearson's correlation coefficient to a
 * file.
 * \param[in] s a header string to write.
 * \param[in] fp the file handle to write to. */
void p_cor_coeff::diagnostic(const char* s,FILE *fp) {
    double ninv=1./n,ex=sx*ninv,ey=sy*ninv,
           varx=sxx*ninv-ex*ex,
           vary=syy*ninv-ey*ey,
           cov=sxy*ninv-ex*ey;
    fprintf(fp,"%s %ld %g %g %g %g %g %g\n",s,n,ex,ey,sqrt(varx),sqrt(vary),cov,cov/sqrt(varx*vary));
}

/** Initializes the wind correlation computution class.
 * \param[in] n_mcts_ the number of MCTS samples.
 * \param[in] gpt_ the number of gliders per trial.
 * \param[in] tf_ a pointer to the turbulent fluid class. */
wind_correl::wind_correl(int gpt_,int n_mcts_,turb_fluid_grid* tf_) :
    gpt(gpt_), n_mcts(n_mcts_), totm(wc_pts*n_mcts*gpt), tf(tf_),
    ptab(new double[6*totm]), wtab(ptab+3*totm), mztab(new double[totm]),
    ctab(new int[gpt*wc_pts]) {

    // Ensure that the turbulent fluid class has enough memory to compute the
    // wind at all of the glider locations encountered during the MCTS
    tf->allocate_vel_table(totm);

    // Initialize the correlation coefficient classes
    for(int i=0;i<wc_pts;i++) wc[i].init();
}

/** The class destructor frees the dynamically allocated memory. */
wind_correl::~wind_correl() {
    delete [] ctab;
    delete [] mztab;
    delete [] ptab;
}

/** Computes the wind at all of the glider locations encountered in MCTS, and
 * then compares them to the GPR predictions.
 * \param[in] time the current simulation time. */
void wind_correl::wind_computation() {

    // Check that the expected number of measurements were taken
    bool measure_err=false;
    for(int i=0;i<gpt;i++) if(ctab[i]!=n_mcts) {
        printf("ctab[%d]=%d   not equal %d\n",i,ctab[i],n_mcts);
        measure_err=true;
    }
    if(measure_err) exit(1);

    // Compute the actual wind field at the glider's position,
    tf->vel_multi(totm,ptab,wtab);
    for(int i=0;i<totm;i+=wc_pts)
        for(int j=0;j<wc_pts;j++) wc[j].add(mztab[i+j],wtab[3*(i+j)+2]);
}

/** Stores the glider position and GPR wind prediction for a later correlation
 * computation.
 * \param[in] co the measurement point index.
 * \param[in] id the glider ID.
 * \param[in] g_ a reference to the glider state.
 * \param[in] wz a GPR prediction for the vertical wind component at the
 *               glider's current position. */
void wind_correl::measure_internal(int co,int id,glider &g_,double wz) {

    // Add the glider position and the GPR prediction of the vertical wind
    // component to the internal tables
    int ind=(n_mcts*id+ctab[wc_pts*id+co]++)*wc_pts+co;
    double *pp=ptab+3*ind;
    *pp=g_.rx;
    pp[1]=g_.ry;
    pp[2]=g_.rz;
    mztab[ind]=wz;
}

/** Outputs the correlation information at the measurement points.
 * \param[in] fp a file handle to write to. */
void wind_correl::output(FILE *fp) {
    const char* ti[wc_pts]={"0","0.1","0.2","0.5","1","2","5"};
    for(int k=0;k<wc_pts;k++)
        wc[k].diagnostic(ti[k],fp);
}
