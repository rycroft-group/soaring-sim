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

/** This test verify the gpr prediction accuracy for random gliders in dynamic field. */

// Multi glider
// Compare known field gliders under static field and dynamic field
// Only mcts is used, don't need gpr

int main(int argc, char **argv){

    if(argc!=3){
        fputs("Syntax: ./gpr_test5 <integer> <float> \n",stderr);
        return 1;
    }

    int memory=atoi(argv[1]);
    double time_scale=atof(argv[2]);

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
    int depth=1000;
    int nmode=64;
    const int ngl=20;
    int action_scale=50;
    const double box_size=50.;
    double alpha=M_PI*M_PI/(9.*(pow(sqrt(3.)*M_PI,-2/3.)-pow(sqrt(3.)*M_PI*nmode,-2/3.)));

    bool frozen=true;
    double wind_val=0.05;
    double* r=new double[3*ngl];
    double* wind=new double[3*ngl];
    bool debug=false;

    // Compute a time step for wind forward
    const double dt_pad=0.06;
    const double gl_dt=0.01;

    // Create file to store values
    char fn[100];
    sprintf(fn,"gpr_test_result_frozen_%d/gpr_test_t%d_m%d.dat",nmode,int(time_scale),memory);
    FILE *fp=fopen(fn,"wb+");
    if(fp==NULL){
	fputs("Can't open file\n",stderr);
	exit(1);
    }

    for(int l=0;l<10;l++){
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
            gp_test[o]->init(100,50,sqrt(3.),50.);
            //gp_test[o]->init();
            gl_kf[o]=new game_glider(1.,1/15.,wind_val,wind_state,frozen,tf,*gp_test[o]);
            // Start at random locations
            double rx_,ry_,rz_,ux_,uy_,uz_,vx_,vy_,vz_,wx_,wy_,wz_,bank_;
            rx_=box_size*gsl_rng_uniform(rng);
            ry_=box_size*gsl_rng_uniform(rng);
            rz_=box_size*gsl_rng_uniform(rng);

            gl_kf[o]->init(rx_,ry_,rz_,1.,0.,-1./15);

        // Open file and store initial positions
        }
     // Compute a time step for wind forward
        double tf_dt=dt_pad*tf.est_max_timestep();
        double sint=action_scale*gl_kf[0]->dt;
        int tf_l=static_cast<int>(sint/tf_dt)+1;
        tf_dt=sint/tf_l;
            int i=1;
            while(i<6000){
                tf.transform();
            //printf("i %d\n",i);
                if(debug) printf("i %d\n",i);
                #pragma omp parallel for schedule(dynamic)

                for(int o=0;o<ngl;o++){

                    // Random actions
                   int action=gsl_rng_uniform_int(rng,action_size);
                    gl_kf[o]->bank+=gl_kf[o]->mu_f(action);
            }
            // Random tf glider play
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
                    // After the memory is full, do prediction
                    if(k==action_scale-1 && i<5000){
                        double time=0.;
                        if(!frozen) time=i*gl_kf[o]->dt;
                        // printf("time %g i %d\n",time,i);
                        gp_test[o]->add_measurement(r[3*o]/box_size,r[3*o+1]/box_size,r[3*o+2]/box_size,time,wind_val*wind[3*o],wind_val*wind[3*o+1],wind_val*wind[3*o+2]);
                    }
                    // After the memory is full, do prediction
                    if(i>5000){
			double wx_,wy_,wz_;
			double time=0.;
			if(!frozen) time=i*gl_kf[o]->dt;
                        gp_test[o]->predict(r[3*o]/box_size,r[3*o+1]/box_size,r[3*o+2]/box_size,time,wx_,wy_,wz_);
                //        printf("%g %g %g %g %g %g\n",wind_val*wind[3*o],wind_val*wind[3*o+1],wind_val*wind[3*o+2],wx_,wy_,wz_);
			fprintf(fp,"%d %g %g %g %g %g %g\n",i,wind_val*wind[3*o],wind_val*wind[3*o+1],wind_val*wind[3*o+2],wx_,wy_,wz_);
                    }
		    // Perform actions with known wind
                    gl_kf[o]->wind=3;
                    if(debug) printf("wx %g wy %g wz %g\n",wind[3*o],wind[3*o+1],wind[3*o+2]);
                    double w=gl_kf[o]->play(0,1,i);
                    // fprintf(fp2[o],"%g %g\n",max_w,w);
                    gl_kf[o]->wind=wind_state;

                }
               if(!frozen) tf.step_forward(gl_kf[0]->dt);
            }

        // wind forward
        if(!frozen){
            // for(int o=0;o<tf_l;o++) tf.step_forward(tf_dt);
            }

        }

        // After the memory is full, do prediction

        for(int o=0;o<ngl;o++){
            delete gl_kf[o];
            delete gp_test[o];
        }
        delete [] gl_kf;
        delete [] gp_test;

    }
     if(fp!=NULL) fclose(fp);
    delete [] wind;
    delete [] r;
    double t1=wtime();
    printf("running time is %g\n",t1-t0);
}
