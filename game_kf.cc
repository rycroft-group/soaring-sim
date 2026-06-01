#include "game.hh"
#include "game_mcts.hh"
#include "turb_fluid.hh"
#include "gpr.hh"
#include <gsl/gsl_rng.h>
#include "common.hh"
#include <cstring>
#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef _OPENMP
#include "omp.h"
#endif

// Multi glider
// Compare known field gliders under static field and dynamic field
// Only mcts is used, don't need gpr

int main(){
    double t0=wtime();
    gsl_rng *rng;
    rng=gsl_rng_alloc(gsl_rng_taus2);
   // turb_fluid tf(256,256,256,0,50,0,50,0,50,1.,10.);
   // tf.init_steady_state();
    int action_size=3;
    int wind_state=0;
    int memory=50;
    int depth=1000;
    int nmode=64;
    const int ngl=10;
    int action_scale=50;
    const double box_size=50.;
    double time_scale=128.;
    double alpha=M_PI*M_PI/(9.*(pow(sqrt(3.)*M_PI,-2/3.)-pow(sqrt(3.)*M_PI*nmode,-2/3.)));

    bool frozen=false;
    double wind_val=0.05;
    double* r=new double[3*ngl];
    double* wind=new double[3*ngl];
    char fn[100];
    char buffer[150];
    bool debug=false;

    // Compute a time step for wind forward
    const double dt_pad=0.06;
    const double gl_dt=0.01;

    for(int l=0;l<20;l++){
        printf("tf %d m %d d %d frozen %d wind_val %g time scale %g \n",l,memory,depth,frozen,wind_val,time_scale);
        turb_fluid_grid tf(nmode,nmode,nmode,0,box_size,0,box_size,0,box_size,time_scale*3.094,alpha);
        tf.init_steady_state();
        // Save the turbulent fluid field
        // Output the initial cross-section
        //tf.transform();

        gpr **gp_test=new gpr*[ngl];
        game_glider **gl_kf=new game_glider*[ngl];
            FILE** fp1=new FILE*[ngl];
            FILE** fp2=new FILE*[ngl];

        //FILE** fp=static_cast<FILE**>(malloc(sizeof(FILE*)*ngl));
        for(int o=0;o<ngl;o++){
            // Initialize gaussion progress regressor
            gp_test[o]=new gpr(memory,time_scale);
            // gp[o]->init(100,50,sqrt(3.),time_scale*4.);
            gp_test[o]->init(); // not used
            gl_kf[o]=new game_glider(1.,1/15.,wind_val,2,frozen,tf,*gp_test[o]);
            // Start at random locations
            double rx_,ry_,rz_,ux_,uy_,uz_,vx_,vy_,vz_,wx_,wy_,wz_,bank_;
            rx_=box_size*gsl_rng_uniform(rng);
            ry_=box_size*gsl_rng_uniform(rng);
            rz_=box_size*gsl_rng_uniform(rng);

            gl_kf[o]->init(rx_,ry_,rz_,1.,0.,-1./15);

            // Open file and store initial positions
            sprintf(fn,"test/dynamic_time_scale4/t%d_d%d_a%d/tf_result_%d",int(time_scale),depth,action_scale,l*ngl+o);
                mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
            sprintf(buffer,"%s/game_mcts_kf.dat",fn);
            fp1[o]=safe_fopen(buffer,"w");
                sprintf(buffer,"%s/game_mcts_kf_sim.dat",fn);
            fp2[o]=safe_fopen(buffer,"w");
            gl_kf[o]->output(0,0,fp1[o]);
        }
         // Compute a time step for wind forward
        double tf_dt=dt_pad*tf.est_max_timestep();
        double sint=action_scale*gl_kf[0]->dt;
        printf("%g\n",tf_dt);
        return 0;
        int tf_l=static_cast<int>(sint/tf_dt)+1;
        tf_dt=sint/tf_l;
            int i=1;
            while(i<50000){
                              tf.transform();
                //printf("i %d\n",i);
                if(debug) printf("i %d\n",i);
             #pragma omp parallel for schedule(dynamic)
            // Do MCTS for known tf
                for(int o=0;o<ngl;o++){
                              double rx_,ry_,rz_,ux_,uy_,uz_,vx_,vy_,vz_,wx_,wy_,wz_,bank_;
                    game_mcts mcts(action_size,depth,*gl_kf[o]);

                    // Glider move
                    mcts.reset();
                    // Store the current status
                    rx_=gl_kf[o]->rx,ry_=gl_kf[o]->ry,rz_=gl_kf[o]->rz;
                    ux_=gl_kf[o]->ux,uy_=gl_kf[o]->uy,uz_=gl_kf[o]->uz;
                    vx_=gl_kf[o]->vx,vy_=gl_kf[o]->vy,vz_=gl_kf[o]->vz;
                    wx_=gl_kf[o]->wx,wy_=gl_kf[o]->wy,wz_=gl_kf[o]->wz;
                    bank_=gl_kf[o]->bank;
                    if(debug) printf("vx %g vy %g vz %g\n",vx_,vy_,vz_);
                    // if(i%action_scale==0) g.add_measurement(rx_,ry_,rz_,0.,wx_,wy_,wz_);
                    if(debug) printf("rx %g ry %g rz %g\n",rx_,ry_,rz_);
                    // Monte Carlo tree search
                    bool out=false; // Output simulated path
                    for(int j=0;j<10000;j++) {
                //          if(debug) puts("mcts step");
                          out=false;
                        if((j%1000==0)&&(i%1000==1)) out=true;
                                    mcts.simulate_path_glider(1,i,out,fp2[o]);
                        for(int k=0;k<action_size;k++){
                            // fprintf(fp7[o],"%g\n",mcts.t[k].weight(0)/mcts.t[k].n);

                          }
                    }
                   int action=mcts.pick_best(0);
                   double max_w=mcts.t[action].w_step/(mcts.t[action].n);

                    // Restore the current status
                    gl_kf[o]->reset(rx_,ry_,rz_,ux_,uy_,uz_,vx_,vy_,vz_,bank_);

                    // action time scale = dt*10
                    gl_kf[o]->bank+=gl_kf[o]->mu_f(action);
            }
            // Known tf glider play
            for(int k=0;k<action_scale;k++,i++){
                for(int o=0;o<ngl;o++){
                    r[3*o]=gl_kf[o]->rx;
                    r[3*o+1]=gl_kf[o]->ry;
                    r[3*o+2]=gl_kf[o]->rz;
                }
                // Calculate the velocities at the random sample points
                tf.allocate_vel_table(ngl);
                tf.vel_multi(ngl,r,wind);
                for(int o=0;o<ngl;o++){
                    gl_kf[o]->wx=wind_val*wind[3*o];
                    gl_kf[o]->wy=wind_val*wind[3*o+1];
                    gl_kf[o]->wz=wind_val*wind[3*o+2];
                    gl_kf[o]->wind=3;
                    if(debug) printf("wx %g wy %g wz %g\n",wind[3*o],wind[3*o+1],wind[3*o+2]);
                    double w=gl_kf[o]->play(0,1,i);
                    // fprintf(fp2[o],"%g %g\n",max_w,w);
                    gl_kf[o]->wind=2;
                    // Store the move
                    gl_kf[o]->output(i+1,0,fp1[o]);
                }
                       tf.step_forward(gl_kf[0]->dt);
            }

            // wind forward
            if(!frozen){
                    //for(int o=0;o<tf_l;o++) tf.step_forward(tf_dt);
                    }

            }

        for(int o=0;o<ngl;o++){
            delete gl_kf[o];
            delete gp_test[o];
            if(fp1[o]!=NULL) fclose(fp1[o]);
            if(fp2[o]!=NULL) fclose(fp2[o]);
        }
        delete [] gl_kf;
        delete [] gp_test;
        delete [] fp2;
        delete [] fp1;

    }
    delete [] wind;
    delete [] r;
    double t1=wtime();
    printf("running time is %g\n",t1-t0);
}
