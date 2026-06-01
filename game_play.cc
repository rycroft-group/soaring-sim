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

int main(){

    double t0=wtime();
    gsl_rng *rng;
    rng=gsl_rng_alloc(gsl_rng_taus2);
   // turb_fluid tf(256,256,256,0,50,0,50,0,50,1.,10.);
   // tf.init_steady_state();
    int action_size=3;
    int wind_state=2;
    int memory=50;
    int depth=1000;
    double wind_val=0.05;
    bool frozen=true;

#pragma omp parallel for
    for(int l=0;l<20;l++){
        // Initialize gaussion progress regressor
        gpr g(memory);
        g.init();
        turb_fluid_grid tf(256,256,256,0,50,0,50,0,50,0.67,3.47);
        tf.init_steady_state();
        // Initialize glider
        game_glider ttt(1.,1/15.,wind_val,wind_state,frozen,tf,g);
        // Start at random locations
        double rx_,ry_,rz_,ux_,uy_,uz_,vx_,vy_,vz_,wx_,wy_,wz_,bank_;
        rx_=25+25.*(2*gsl_rng_uniform(rng)-1);
        ry_=25+25.*(2*gsl_rng_uniform(rng)-1);
        rz_=25+25.*(2*gsl_rng_uniform(rng)-1);
        //rx_=30.,ry_=30.,rz_=30.;
    //rz_=0.;
          ttt.init(rx_,ry_,rz_,1.,0.,-1./15);

        char fn[128];
        char buffer[192];
	// Open file and store initial positions
        sprintf(fn,"test_certain_wind/tf_slice/m%d_d%d/tf_result_%d",memory,depth,l);
	    mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
        sprintf(buffer,"%s/game_mcts.dat",fn);
        FILE *fp=safe_fopen(buffer,"w");
        sprintf(buffer,"%s/game_w.dat",fn);
        FILE *fp2=safe_fopen(buffer,"w");
        sprintf(buffer,"%s/game_mcts_sim.dat",fn);
        FILE *fp3=safe_fopen(buffer,"w");
        sprintf(buffer,"%s/game_mcts_policy.dat",fn);
        FILE *fp4=safe_fopen(buffer,"w");
        sprintf(buffer,"%s/game_mcts_random.dat",fn);
        FILE *fp5=safe_fopen(buffer,"w");
        sprintf(buffer,"%s/game_mcts_free.dat",fn);
        FILE *fp6=safe_fopen(buffer,"w");
        sprintf(buffer,"%s/game_var.dat",fn);
        FILE *fp7=safe_fopen(buffer,"w");
//    ttt.output(0,0,fp);

        // Start tree search
        game_mcts mcts(action_size,depth,ttt);
        float w;
        int action;
        for(int i=0;i<200;i++){
            action=gsl_rng_uniform_int(rng,action_size);
            ttt.play(action,1,i);
        }
        // ttt.rx=40.,ttt.ry=40.,ttt.rz=40.;
        ttt.output(0,0,fp);

        // Initialize a new glider using the fixed policy
        gpr g_test(50);
        g_test.init();
        game_glider gl_test(1.,1/15.,wind_val,wind_state,frozen,tf,g_test);
        rx_=ttt.rx,ry_=ttt.ry;rz_=ttt.rz;
        gl_test.init(rx_,ry_,rz_,1.,0.,-1./15);

        // Initialize a new glider using random policy
        game_glider gl_random(1.,1/15.,wind_val,wind_state,frozen,tf,g_test);
        gl_random.init(rx_,ry_,rz_,1.,0.,-1./15);

        // Initialize a new glider using no action policy
        game_glider gl_free(1.,1/15.,wind_val,wind_state,frozen,tf,g_test);
        gl_free.init(rx_,ry_,rz_,1.,0.,-1./15);

        int i=1;
        int action_scale=50;
        while(i<10000){

            printf("time frame %d\n",i);

            // Glider move
            mcts.reset();
            // Store the current status
            rx_=ttt.rx,ry_=ttt.ry,rz_=ttt.rz;
            ux_=ttt.ux,uy_=ttt.uy,uz_=ttt.uz;
            vx_=ttt.vx,vy_=ttt.vy,vz_=ttt.vz;
            wx_=ttt.wx,wy_=ttt.wy,wz_=ttt.wz;
            bank_=ttt.bank;
            // printf("vx %g vy %g vz %g\n",vx_,vy_,vz_);
            // if(i%action_scale==0) g.add_measurement(rx_,ry_,rz_,0.,wx_,wy_,wz_);
            //printf("rx %g rx_ %g ry %g ry_ %g rz %g rz_ %g\n",ttt.rx,rx_,ttt.ry,ry_,ttt.rz,rz_);
            // Monte Carlo tree search
            bool out=false; // Output simulated path
            for(int j=0;j<10000;j++) {
		out=false;
                if((j%1000==0)&&(i%1000==1)) out=true;
                mcts.simulate_path_glider(1,i,out,fp3);
                for(int k=0;k<action_size;k++){
                    fprintf(fp7,"%g\n",mcts.t[k].weight(0)/mcts.t[k].n);

                  }
            }
            action=mcts.pick_best(0);
           double max_w=mcts.t[action].w_step/(mcts.t[action].n);

        //action=0; //Test
            // Restore the current status
            ttt.reset(rx_,ry_,rz_,ux_,uy_,uz_,vx_,vy_,vz_,bank_);

            // Simulate a move using the action chose and store
            ttt.bank+=ttt.mu_f(action);
           for(int k=0;k<action_scale;k++){
              ttt.simulate_play(action,1,i);
              ttt.output(i+1,action,fp3);
            }

            // Restore the current status
            ttt.reset(rx_,ry_,rz_,ux_,uy_,uz_,vx_,vy_,vz_,bank_);
            //printf("rx %g rx_ %g ry %g ry_ %g rz %g rz_ %g\n",ttt.rx,rx_,ttt.ry,ry_,ttt.rz,rz_);
            // Actual move
        //action=gsl_rng_uniform_int(rng,action_size);
            //action=0;

            // action time scale = dt*10
            ttt.bank+=ttt.mu_f(action);
            for(int k=0;k<action_scale;k++,i++){
                w=ttt.play(action,1,i);
                fprintf(fp2,"%g %g\n",max_w,w);

                // Store the move
                ttt.output(i+1,action,fp);
            }

        // Test glider move a step
        // Choose an action use the fixed policy
            int policy_action=gl_test.pick_action_policy(i);

            // Take a move using the action chose and store
            gl_test.bank+=gl_test.mu_f(policy_action);
            for(int k=0;k<action_scale;k++){
              gl_test.play(policy_action,1,i);
              gl_test.output(i+1,policy_action,fp4);
            }

        policy_action=gsl_rng_uniform_int(rng,action_size);
    //policy_action=0;
            // Take a move using the action chose and store
            gl_random.bank+=gl_random.mu_f(policy_action);
        //    gl_random.bank=gl_random.mu_f(1)*4;
    for(int k=0;k<action_scale;k++){
              gl_random.play(policy_action,1,i);
              gl_random.output(i+1,policy_action,fp5);
        }

        policy_action=1;

            // Take a move using the action chose and store
            gl_free.bank=0.;
            for(int k=0;k<action_scale;k++){
              gl_free.play(policy_action,1,i);
              gl_free.output(i+1,policy_action,fp6);
	    }

        }
    if(fp!=NULL) fclose(fp);
    if(fp2!=NULL) fclose(fp2);
    if(fp3!=NULL) fclose(fp3);
    if(fp4!=NULL) fclose(fp4);
    if(fp5!=NULL) fclose(fp5);
    if(fp6!=NULL) fclose(fp6);
    if(fp7!=NULL) fclose(fp7);
    }
    double t1=wtime();
    printf("running time is %g\n",t1-t0);
}
