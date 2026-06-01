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
    double time_scale=1.;
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
            gp[o]=new gpr(memory);
            // gp[o]->init(100,50,sqrt(3.),time_scale*4.);
            gp[o]->init(); // not used
            gl[o]=new game_glider(1.,1/15.,wind_val,3,frozen,tf,*gp[o]);
            // Start at random locations
            double rx_,ry_,rz_,ux_,uy_,uz_,vx_,vy_,vz_,wx_,wy_,wz_,bank_;
            rx_=box_size*gsl_rng_uniform(rng);
            ry_=box_size*gsl_rng_uniform(rng);
            rz_=box_size*gsl_rng_uniform(rng);

            gl[o]->init(rx_,ry_,rz_,1.,0.,-1./15);

    	// Open file and store initial positions
            sprintf(fn,"test/random/t%d_d%d_a%d/tf_result_%d",int(time_scale),depth,action_scale,l*ngl+o);
    	    mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
            sprintf(buffer,"%s/game_random.dat",fn);
            fp1[o]=safe_fopen(buffer,"w");
	        sprintf(buffer,"%s/game_random_sim.dat",fn);
            fp2[o]=safe_fopen(buffer,"w");
            gl[o]->output(0,0,fp1[o]);
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
            // Take a random step
                for(int o=0;o<ngl;o++){

                   int action=gsl_rng_uniform_int(rng,action_size);
                    gl[o]->bank+=gl[o]->mu_f(action);
            }
            // Random glider play
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
		    if(debug) printf("wx %g wy %g wz %g\n",wind[3*o],wind[3*o+1],wind[3*o+2]);
                    double w=gl[o]->play(0,1,i);
                    // fprintf(fp2[o],"%g %g\n",max_w,w);
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
