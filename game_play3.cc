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

int main(int argc, char **argv){

    if(argc!=4){
        fputs("Syntax: ./game_play3 <integer> <integer> <float> \n",stderr);
        return 1;
    }

    int memory=atoi(argv[1]);
    int nmode=atoi(argv[2]);
    double time_scale=atof(argv[3]);

    // Perform some basic checks on the numbers
    if(memory<0) {
        fputs("<integer> must be positive\n",stderr);
        return 1;
    }
    if(time_scale<0) {
        fputs("<float> must be positive\n",stderr);
        return 1;
    }

    double t0=wtime();
    gsl_rng *rng;
    rng=gsl_rng_alloc(gsl_rng_taus2);
   // turb_fluid tf(256,256,256,0,50,0,50,0,50,1.,10.);
   // tf.init_steady_state();
    int action_size=3;
    int wind_state=0;
    // int memory=25;
    int depth=1000;
    //int nmode=64;
    const int ngl=10;
    int action_scale=50;
    const double box_size=50.;
    // double time_scale=1.;
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

        gpr **gp=new gpr*[ngl];
        game_glider **gl=new game_glider*[ngl];
    	FILE** fp1=new FILE*[ngl];
    	FILE** fp2=new FILE*[ngl];

        //FILE** fp=static_cast<FILE**>(malloc(sizeof(FILE*)*ngl));
        for(int o=0;o<ngl;o++){
            // Initialize gaussion progress regressor
            gp[o]=new gpr(memory,time_scale);
            gp[o]->init(100,50,sqrt(3.),50);
            //gp[o]->init();
            gl[o]=new game_glider(1.,1/15.,wind_val,wind_state,frozen,tf,*gp[o]);
            // Start at random locations
            double rx_,ry_,rz_,ux_,uy_,uz_,vx_,vy_,vz_,wx_,wy_,wz_,bank_;
            rx_=box_size*gsl_rng_uniform(rng);
            ry_=box_size*gsl_rng_uniform(rng);
            rz_=box_size*gsl_rng_uniform(rng);

            gl[o]->init(rx_,ry_,rz_,1.,0.,-1./15);

    	// Open file and store initial positions
            sprintf(fn,"new_test/memory_%d/t%d_d%d_a%d_n%d/tf_result_%d",memory,int(time_scale),depth,action_scale,nmode,l*ngl+o);
    	    mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
            sprintf(buffer,"%s/game_mcts_kf.dat",fn);
            fp1[o]=safe_fopen(buffer,"w");
	        sprintf(buffer,"%s/game_mcts_kf_sim.dat",fn);
            fp2[o]=safe_fopen(buffer,"w");
            gl[o]->output(0,0,fp1[o]);
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
	 // Compute a time step for wind forward
        double tf_dt=dt_pad*tf.est_max_timestep();
	double sint=action_scale*gl[0]->dt;
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
                    game_mcts mcts(action_size,depth,*gl[o]);

                    // Glider move
                    mcts.reset();
                    // Store the current status
                    rx_=gl[o]->rx,ry_=gl[o]->ry,rz_=gl[o]->rz;
                    ux_=gl[o]->ux,uy_=gl[o]->uy,uz_=gl[o]->uz;
                    vx_=gl[o]->vx,vy_=gl[o]->vy,vz_=gl[o]->vz;
                    wx_=gl[o]->wx,wy_=gl[o]->wy,wz_=gl[o]->wz;
                    bank_=gl[o]->bank;
                    if(debug) printf("vx %g vy %g vz %g\n",vx_,vy_,vz_);
                    // if(i%action_scale==0) g.add_measurement(rx_,ry_,rz_,0.,wx_,wy_,wz_);
                    if(debug) printf("rx %g ry %g rz %g\n",rx_,ry_,rz_);
                    // Monte Carlo tree search
                    bool out=false; // Output simulated path
                    for(int j=0;j<10000;j++) {
		//	  if(debug) puts("mcts step");
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
                    gl[o]->reset(rx_,ry_,rz_,ux_,uy_,uz_,vx_,vy_,vz_,bank_);

                    // action time scale = dt*10
                    gl[o]->bank+=gl[o]->mu_f(action);
            }
            // Known tf glider play
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
		    if(debug) printf("wx %g wy %g wz %g\n",wind[3*o],wind[3*o+1],wind[3*o+2]);
                    double w=gl[o]->play(0,1,i);
                    // fprintf(fp2[o],"%g %g\n",max_w,w);
                    gl[o]->wind=wind_state;
                    // Store the move
                    gl[o]->output(i+1,0,fp1[o]);
                }
		       tf.step_forward(gl[0]->dt);
            }

	    // wind forward
	    if(!frozen){
		    //for(int o=0;o<tf_l;o++) tf.step_forward(tf_dt);
	   	 }

            }

        for(int o=0;o<ngl;o++){
            delete gl[o];
            delete gp[o];
            if(fp1[o]!=NULL) fclose(fp1[o]);
            if(fp2[o]!=NULL) fclose(fp2[o]);
        }
        delete [] gl;
        delete [] gp;
        delete [] fp2;
	delete [] fp1;

    }
    delete [] wind;
    delete [] r;
    double t1=wtime();
    printf("running time is %g\n",t1-t0);
}
