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
    int wind_state=0;
    int memory=50;
    int depth=1000;
    int nmode=128;
    const int ngl=10;
    int action_scale=50;
    const double box_size=20.;
    double time_scale=256.;

    bool frozen=false;
    double wind_val=0.05;
    double* r=new double[3*ngl];
    double* wind=new double[3*ngl];
    char fn[100];
    char buffer[150];
    bool debug=false;

    for(int l=0;l<20;l++){
        printf("tf %d m %d d %d frozen %d wind_val %g time scale %g \n",l,memory,depth,frozen,wind_val,time_scale);
        turb_fluid_grid tf(nmode,nmode,nmode,0,box_size,0,box_size,0,box_size,time_scale*3.125,3.47);
        tf.init_steady_state();
        // Save the turbulent fluid field
        // Output the initial cross-section
        tf.transform();
        /**
	sprintf(fn,"tf.odr");
        mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
        sprintf(buffer,"%s/ux.0",fn);
        tf.output_z(buffer,0,0);
        sprintf(buffer,"%s/uy.0",fn);
        tf.output_z(buffer,1,0);
        sprintf(buffer,"%s/uz.0",fn);
        tf.output_z(buffer,2,0);
        puts("# Output frame 0");

        // Integrate the modes and output additional snapshots
#pragma omp parallel
{
	#pragma omp for
	for(int k=1;k<=nmode;k++) {
            // Output a cross-section
            tf.transform();
            sprintf(buffer,"%s/ux.%d",fn,k);
            tf.output_z(buffer,0,k);
            sprintf(buffer,"%s/uy.%d",fn,k);
            tf.output_z(buffer,1,k);
            sprintf(buffer,"%s/uz.%d",fn,k);
            tf.output_z(buffer,2,k);
        }
}
        puts("# Finish Output frames");
*/
        gpr **gp=new gpr*[ngl];
        gpr **gp_test=new gpr*[ngl];
        game_glider **gl=new game_glider*[ngl];
        game_glider **gl_policy=new game_glider*[ngl];
        game_glider **gl_random=new game_glider*[ngl];
        game_glider **gl_free=new game_glider*[ngl];
        game_glider **gl_kf=new game_glider*[ngl];
    	FILE** fp=new FILE*[ngl];
    	FILE** fp3=new FILE*[ngl];
    	FILE** fp4=new FILE*[ngl];
    	FILE** fp5=new FILE*[ngl];
    	FILE** fp6=new FILE*[ngl];
        FILE** fp7=new FILE*[ngl];
	    FILE** fp8=new FILE*[ngl];

        //FILE** fp=static_cast<FILE**>(malloc(sizeof(FILE*)*ngl));
        for(int o=0;o<ngl;o++){
            // Initialize gaussion progress regressor
            gp[o]=new gpr(memory);
            gp_test[o]=new gpr(memory);
            gp[o]->init(100,50,sqrt(3.),time_scale*4.);
            gp_test[o]->init(100,50,sqrt(3.),time_scale*4.);
            gl[o]=new game_glider(1.,1/15.,wind_val,wind_state,frozen,tf,*gp[o]);
            gl_policy[o]=new game_glider(1.,1/15.,wind_val,wind_state,frozen,tf,*gp_test[o]);
            gl_random[o]=new game_glider(1.,1/15.,wind_val,wind_state,frozen,tf,*gp_test[o]);
            gl_free[o]=new game_glider(1.,1/15.,wind_val,wind_state,frozen,tf,*gp_test[o]);
            gl_kf[o]=new game_glider(1.,1/15.,wind_val,2,frozen,tf,*gp_test[o]);
            // Start at random locations
            double rx_,ry_,rz_,ux_,uy_,uz_,vx_,vy_,vz_,wx_,wy_,wz_,bank_;
            rx_=box_size*gsl_rng_uniform(rng);
            ry_=box_size*gsl_rng_uniform(rng);
            rz_=box_size*gsl_rng_uniform(rng);
            //rx_=30.,ry_=30.,rz_=30.;
        //rz_=0.;
            gl[o]->init(rx_,ry_,rz_,1.,0.,-1./15);

    	// Open file and store initial positions
            sprintf(fn,"test/dynamic_time_scale2/t%d_d%d_a%d/tf_result_%d",int(time_scale),depth,action_scale,l*ngl+o);
    	    mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
            sprintf(buffer,"%s/game_mcts.dat",fn);
            fp[o]=safe_fopen(buffer,"w");
            //sprintf(buffer,"%s/game_w.dat",fn);
            // FILE *fp2=safe_fopen(buffer,"w");
            sprintf(buffer,"%s/game_mcts_sim.dat",fn);
            fp3[o]=safe_fopen(buffer,"w");
            sprintf(buffer,"%s/game_mcts_policy.dat",fn);
            fp4[o]=safe_fopen(buffer,"w");
            sprintf(buffer,"%s/game_mcts_random.dat",fn);
            fp5[o]=safe_fopen(buffer,"w");
            sprintf(buffer,"%s/game_mcts_free.dat",fn);
            fp6[o]=safe_fopen(buffer,"w");
            sprintf(buffer,"%s/game_mcts_kf.dat",fn);
            fp7[o]=safe_fopen(buffer,"w");
	    sprintf(buffer,"%s/game_mcts_kf_sim.dat",fn);
            fp8[o]=safe_fopen(buffer,"w");
            // sprintf(buffer,"%s/game_var.dat",fn);
            // FILE *fp7=safe_fopen(buffer,"w");
    //    ttt.output(0,0,fp);
            // Start tree search
            float w;
            int action;
            for(int i=0;i<10;i++){
                action=gsl_rng_uniform_int(rng,action_size);
                gl[o]->play(action,1,i);
            }
            // Initialize a new glider using the fixed policy
            rx_=gl[o]->rx,ry_=gl[o]->ry;rz_=gl[o]->rz;
            gl_policy[o]->init(rx_,ry_,rz_,1.,0.,-1./15);

            // Initialize a new glider using random policy
            gl_random[o]->init(rx_,ry_,rz_,1.,0.,-1./15);

            // Initialize a new glider using no action policy
            gl_free[o]->init(rx_,ry_,rz_,1.,0.,-1./15);
            // Initialize a new glider using known field policy
            gl_kf[o]->init(rx_,ry_,rz_,1.,0.,-1./15);

	        gl[o]->output(0,0,fp[o]);
            gl_policy[o]->output(0,0,fp4[o]);
            gl_random[o]->output(0,0,fp5[o]);
            gl_free[o]->output(0,0,fp6[o]);
            gl_kf[o]->output(0,0,fp7[o]);
        }

            // Get current wind status and add measurement
            for(int o=0;o<ngl;o++){
                r[3*o]=gl[o]->rx;
                r[3*o+1]=gl[o]->ry;
                r[3*o+2]=gl[o]->rz;
            }
            // Calculate the velocities at the random sample points
            tf.allocate_vel_table(ngl);
            tf.vel_multi(ngl,r,wind);
            for(int o=0;o<ngl;o++){
                gl[o]->wx=wind_val*wind[3*o];
                gl[o]->wy=wind_val*wind[3*o+1];
                gl[o]->wz=wind_val*wind[3*o+2];
                gp[o]->add_measurement(r[3*o]/box_size,r[3*o+1]/box_size,r[3*o+2]/box_size,0.,wind_val*wind[3*o],wind_val*wind[3*o+1],wind_val*wind[3*o+2]);

            }
            int i=1;
            while(i<50000){
		if(debug) printf("i %d\n",i);
            // Do MCTS
            #pragma omp parallel for schedule(dynamic)
                for(int o=0;o<ngl;o++){
			        double rx_,ry_,rz_,ux_,uy_,uz_,vx_,vy_,vz_,wx_,wy_,wz_,bank_;
		            game_mcts mcts(action_size,depth,*gl[o]);

                    // Glider move
                    mcts.reset();
                    // Store the current status
                    rx_=gl[o]->rx,ry_=gl[o]->ry,rz_=gl[o]->rz;
                    ux_=gl[o]->ux,uy_=gl[o]->uy,uz_=gl[o]->uz;
                    vx_=gl[o]->vx,vy_=gl[o]->vy,vz_=gl[o]->vz;
                    wx_=gl[o]->wx,wy_=gl[o]->wy,wz_=gl[o]->wz;
                    bank_=gl[o]->bank;
                    // printf("vx %g vy %g vz %g\n",vx_,vy_,vz_);
                    // if(i%action_scale==0) g.add_measurement(rx_,ry_,rz_,0.,wx_,wy_,wz_);
                //    printf("rx %g rx_ %g ry %g ry_ %g rz %g rz_ %g\n",ttt.rx,rx_,ttt.ry,ry_,ttt.rz,rz_);
                    // Monte Carlo tree search
                    bool out=false; // Output simulated path
                    for(int j=0;j<5000;j++) {
		//	    printf("i %d j %d \n",i,j);
        		          out=false;
                        if((j%1000==0)&&(i%1000==1)) out=true;
                        mcts.simulate_path_glider(1,i,out,fp3[o]);
                        for(int k=0;k<action_size;k++){
                            // fprintf(fp7[o],"%g\n",mcts.t[k].weight(0)/mcts.t[k].n);

                          }
                    }
                   int action=mcts.pick_best(0);
                   double max_w=mcts.t[action].w_step/(mcts.t[action].n);

                //action=0; //Test
                    // Restore the current status
                    gl[o]->reset(rx_,ry_,rz_,ux_,uy_,uz_,vx_,vy_,vz_,bank_);
                    //printf("rx %g rx_ %g ry %g ry_ %g rz %g rz_ %g\n",ttt.rx,rx_,ttt.ry,ry_,ttt.rz,rz_);
                    // Actual move
                //action=gsl_rng_uniform_int(rng,action_size);
                    //action=0;

                    // action time scale = dt*10
                    gl[o]->bank+=gl[o]->mu_f(action);
            }

		if(debug) puts("finish MCTS\n");
            // Do MCTS for known tf
            //#pragma omp parallel for schedule(dynamic)
                for(int o=0;o<ngl;o++){
                  tf.transform();
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
                     //printf("vx %g vy %g vz %g\n",vx_,vy_,vz_);
                    // if(i%action_scale==0) g.add_measurement(rx_,ry_,rz_,0.,wx_,wy_,wz_);
                    //printf("rx %g rx_ %g ry %g ry_ %g rz %g rz_ %g\n",ttt.rx,rx_,ttt.ry,ry_,ttt.rz,rz_);
                    // Monte Carlo tree search
                    bool out=false; // Output simulated path
                    for(int j=0;j<5000;j++) {
                          out=false;
                        if((j%1000==0)&&(i%1000==1)) out=true;
			mcts.simulate_path_glider(1,i,out,fp8[o]);
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
            // Do play
            for(int k=0;k<action_scale;k++,i++){
                for(int o=0;o<ngl;o++){
                    r[3*o]=gl[o]->rx;
                    r[3*o+1]=gl[o]->ry;
                    r[3*o+2]=gl[o]->rz;
                }
                // Calculate the velocities at the random sample points
                tf.allocate_vel_table(ngl);
                tf.vel_multi(ngl,r,wind);
                for(int o=0;o<ngl;o++){
                    gl[o]->wx=wind_val*wind[3*o];
                    gl[o]->wy=wind_val*wind[3*o+1];
                    gl[o]->wz=wind_val*wind[3*o+2];
		            if(k==action_scale-1){
                        double time=0.;
                        if(!frozen) time=i*gl[o]->dt;
                        // printf("time %g i %d\n",time,i);
                        gp[o]->add_measurement(r[3*o]/box_size,r[3*o+1]/box_size,r[3*o+2]/box_size,time,wind_val*wind[3*o],wind_val*wind[3*o+1],wind_val*wind[3*o+2]);
                    }
		            gl[o]->wind=3;
                    double w=gl[o]->play(0,1,i);
                    // fprintf(fp2[o],"%g %g\n",max_w,w);
                    gl[o]->wind=wind_state;
                    // Store the move
                    gl[o]->output(i+1,0,fp[o]);
                }
            }
	    if(debug) puts("finish play\n");
            // Known tf glider
            for(int k=0;k<action_scale;k++){
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
                    double w=gl_kf[o]->play(0,1,i);
                    // fprintf(fp2[o],"%g %g\n",max_w,w);
                    gl_kf[o]->wind=2;
                    // Store the move
                    gl_kf[o]->output(i+1,0,fp7[o]);
                }
            }
            // Policy glider move a step
            // Choose an action use the fixed policy
            int policy_action;
            for(int o=0;o<ngl;o++){
                policy_action=gl_policy[o]->pick_action_policy(i);
                gl_policy[o]->bank=gl_policy[o]->bank+gl_policy[o]->mu_f(policy_action);
            }
            for(int k=0;k<action_scale;k++){
                for(int o=0;o<ngl;o++){
                    r[3*o]=gl_policy[o]->rx;
                    r[3*o+1]=gl_policy[o]->ry;
                    r[3*o+2]=gl_policy[o]->rz;
                }
                // Calculate the velocities at the random sample points
                tf.allocate_vel_table(ngl);
                tf.vel_multi(ngl,r,wind);
                for(int o=0;o<ngl;o++){
                    gl_policy[o]->wx=wind_val*wind[3*o];
                    gl_policy[o]->wy=wind_val*wind[3*o+1];
                    gl_policy[o]->wz=wind_val*wind[3*o+2];
                    gl_policy[o]->wind=3;
                    double w=gl_policy[o]->play(0,1,i);
                    // fprintf(fp2[o],"%g %g\n",max_w,w);
                    gl_policy[o]->wind=wind_state;
                    // Store the move
                    gl_policy[o]->output(i+1,0,fp4[o]);
                }
            }

            // Random glider
            for(int o=0;o<ngl;o++){
                policy_action=gsl_rng_uniform_int(rng,action_size);
                gl_random[o]->bank=(gl_random[o]->bank)+(gl_random[o]->mu_f(policy_action));
            }

            for(int k=0;k<action_scale;k++){
                for(int o=0;o<ngl;o++){
                    r[3*o]=gl_random[o]->rx;
                    r[3*o+1]=gl_random[o]->ry;
                    r[3*o+2]=gl_random[o]->rz;
                }
                // Calculate the velocities at the random sample points
                tf.allocate_vel_table(ngl);
                tf.vel_multi(ngl,r,wind);
                for(int o=0;o<ngl;o++){
                    gl_random[o]->wx=wind_val*wind[3*o];
                    gl_random[o]->wy=wind_val*wind[3*o+1];
                    gl_random[o]->wz=wind_val*wind[3*o+2];
                    gl_random[o]->wind=3;
                    double w=gl_random[o]->play(0,1,i);
                    // fprintf(fp2[o],"%g %g\n",max_w,w);
                    gl_random[o]->wind=wind_state;
                    // Store the move
                    gl_random[o]->output(i+1,0,fp5[o]);
                }
            }

            // Free glider
            for(int o=0;o<ngl;o++){
                gl_free[o]->bank=0.;
            }

            for(int k=0;k<action_scale;k++){
                for(int o=0;o<ngl;o++){
                    r[3*o]=gl_free[o]->rx;
                    r[3*o+1]=gl_free[o]->ry;
                    r[3*o+2]=gl_free[o]->rz;
                }
                // Calculate the velocities at the random sample points
                tf.allocate_vel_table(ngl);
                tf.vel_multi(ngl,r,wind);
                for(int o=0;o<ngl;o++){
                    gl_free[o]->wx=wind_val*wind[3*o];
                    gl_free[o]->wy=wind_val*wind[3*o+1];
                    gl_free[o]->wz=wind_val*wind[3*o+2];
                    gl_free[o]->wind=3;
                    double w=gl_free[o]->play(0,1,i);
                    // fprintf(fp2[o],"%g %g\n",max_w,w);
                    gl_free[o]->wind=wind_state;
                    // Store the move
                    gl_free[o]->output(i+1,0,fp6[o]);
                }
            }
	    // wind forward
	    if(!frozen) tf.step_forward(action_scale*gl[0]->dt);

            }

        // if(fp!=NULL) fclose(fp);
        // if(fp2!=NULL) fclose(fp2);
        // if(fp3!=NULL) fclose(fp3);
        // if(fp4!=NULL) fclose(fp4);
        // if(fp5!=NULL) fclose(fp5);
        // if(fp6!=NULL) fclose(fp6);
        // if(fp7!=NULL) fclose(fp7);

        for(int o=0;o<ngl;++o){
            delete gp[o];
            delete gl[o];
            delete gl_policy[o];
            delete gl_random[o];
            delete gl_free[o];
            delete gl_kf[o];
            if(fp[o]!=NULL) fclose(fp[o]);
            if(fp3[o]!=NULL) fclose(fp3[o]);
            if(fp4[o]!=NULL) fclose(fp4[o]);
            if(fp5[o]!=NULL) fclose(fp5[o]);
            if(fp6[o]!=NULL) fclose(fp6[o]);
            if(fp7[o]!=NULL) fclose(fp7[o]);
	        if(fp8[o]!=NULL) fclose(fp8[o]);
        }
        delete [] gp;
        delete [] gl;
        delete [] gl_policy;
        delete [] gl_random;
        delete [] gl_free;
        delete [] gl_kf;
        delete [] fp;
        delete [] fp8;
	    delete [] fp7;
        delete [] fp6;
        delete [] fp5;
        delete [] fp4;
        delete [] fp3;

    }
    delete [] wind;
    delete [] r;
    double t1=wtime();
    printf("running time is %g\n",t1-t0);
}
