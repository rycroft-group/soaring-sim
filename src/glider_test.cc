#include "glider_test.hh"

/** Initializes the glider test integration class.
 * \param[in] gm_ a reference to the glider model.
 * \param[in] tf_ a pointer to a turbulent fluid class for wind computations.
 *                If set to a null pointer, then the wind field is treated as
 *                zero.
 * \param[in] bank0_ the initial bank angle in the tests.
 * \param[in] duration_ the duration of the tests.
 * \param[in] nframes_ the number of output frames to save. */
glider_test::glider_test(glider_model &gm_,turb_fluid *tf_,int bank0_,double duration_,int nframes_) :
    gm(gm_), tf(tf_), dt_max(gm.max_timestep()), duration(duration_),
    fr_int(duration/nframes_), bank0(bank0_), nframes(nframes_) {

    // Set initial glider velocity
    double nor=gm.s*sqrt(gm.s);
    vx0=gm.c_L/nor;vz0=-gm.c_D/nor;

    // If there is no turbulent fluid, then set the wind vector to zero
    *w=w[1]=w[2]=0;
    *dw=dw[1]=dw[2]=0;

    if(tf!=NULL) tf->allocate_vel_table(1,true);
}

/** Integrates the glider position forward for a specified duration.
 * \param[in] itype the integration type to use.
 * \param[in] dt_pad the padding factor to apply to the maximum timestep.
 * \param[out] dt the actual timestep used, which is adjusted to ensure that an
 *                integer number fit within the output frame duration.
 * \param[out] g the glider state at the end of the test.
 * \param[in] fp a file handle to write to; if set to NULL, then the routine
 *               does not output anything. */
void glider_test::integrate(integration_type itype,double dt_pad,double &dt,glider &g,FILE *fp) {
    glider g2,g3;

    // Choose the timestep so that an integer number fit within the output
    // frame duration
    int l=timestep_select(dt_pad,dt);

    // Initialize the glider position
    g.init(0,0,0,vx0,0,vz0,bank0);
    if(fp!=NULL) g.output(0,fp);

    // Perform the integration
    for(int k=1;k<=nframes;k++) {

        for(int i=0;i<l;i++) {
            // Compute the wind velocity at the glider's position, and integrate
            // the glider's position forward. For the Euler method, this will
            // complete the integration; for the improved Euler method, this will
            // be a preliminary step.
            if(tf!=NULL) tf->vel(g.rx,g.ry,g.rz,*w,w[1],w[2]);

            itype==it_improv_e?gm.ie_step1(dt,g,g2,g3,w)
                              :gm.euler(dt,g,w);

            // If the turbulent fluid is in use and it is not frozen, then integrate
            // it forward
            if(tf!=NULL&&tf->Cinv!=0) tf->step_forward(dt);

            // For the improved Euler method, compute the wind velocities with the
            // preliminary step applied, and use those to complete the integration
            // step
            if(itype==it_improv_e) {
                if(tf!=NULL) tf->vel(g3.rx,g3.ry,g3.rz,*w,w[1],w[2]);
                gm.ie_step2(dt,g,g2,g3,w);
            }
        }

        // Save the glider's information to the output file
        if(fp!=NULL) g.output(k*fr_int,fp);
    }
}
