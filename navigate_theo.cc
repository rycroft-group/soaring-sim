#include "common.hh"
#include "navigate_theo.hh"
//#include "matio.h"

#include <cstdlib>
#include <cstring>
#include <cmath>

/** Constructs the reinforcement learning navigation class, initializing
 * constants, memory for the odor concentration and reward, and the
 * random number generator.
 * \param[in] (m_,n_) the number of gridpoints in the grid.
 * \param[in] n_block*n_block subspaces for space state
 * \param[in] number of frames for dynamic concentration
 * \param[in] (sx,sy) the size of the grid.
 * \param[in] alpha the learning rate.
 * \param[in] gamma the discount factor.
 * \param[in] filename_ the output directory filename. */
navigate_theo::navigate_theo(int m_,int n_,int n_block_,int n_frame_,double sx,double sy,double alpha,double gamma,const char *filename_)
    : m(m_), n(n_),n_block(n_block_), block_size(n/n_block), n_frame(n_frame_),
    mn(m*n), mnf(m*n*n_frame), dx(sx/m), dy(sy/n), filename(filename_),
    c1(new double[mn]), r1(new double[mn]), c2(new double[mn]),
    r2(new double[mn]), c3(new double[mn]), r3(new double[mn]),
    buf(new char[strlen(filename)+128]),
    rng(gsl_rng_alloc(gsl_rng_taus)), rl(2000,5,alpha,gamma) {}

/** The class destructor frees the dynamically allocated memory and
 * the random number generator. */
navigate_theo::~navigate_theo() {
    gsl_rng_free(rng);
    delete [] buf;
    delete [] r3;
    delete [] c3;
    delete [] r2;
    delete [] c2;
    delete [] r1;
    delete [] c1;
}

/** Sets up the odor concentration and reward field (Gaussian).
 * 1,2,3: (30,40), (90,40), (60,90)
 * 4,5,6: (60,40), (90,90), (30,90)
*/
void navigate_theo::init_conc_and_reward() {
    double *cp1=c1,*rp1=r1;
    double fx,fy;
    for(int j=0;j<n;j++) {
        fy=(0.5+j)*dy-40;
        for(int i=0;i<m;i++,rp1++,cp1++) {
            fx=(0.5+i)*dx-60;
            *cp1=100*exp(-0.005*(fx*fx+fy*fy));
            *rp1=fx*fx+fy*fy<50?1:0;
        }
    }

    double *cp2=c2,*rp2=r2;
    for(int j=0;j<n;j++) {
        fy=(0.5+j)*dy-90;
        for(int i=0;i<m;i++,rp2++,cp2++) {
            fx=(0.5+i)*dx-90;
            *cp2=100*exp(-0.005*(fx*fx+fy*fy));
            *rp2=fx*fx+fy*fy<50?1:0;
        }
    }

    double *cp3=c3,*rp3=r3;
    for(int j=0;j<n;j++) {
        fy=(0.5+j)*dy-90;
        for(int i=0;i<m;i++,rp3++,cp3++) {
            fx=(0.5+i)*dx-30;
            *cp3=100*exp(-0.005*(fx*fx+fy*fy));
            *rp3=fx*fx+fy*fy<50?1:0;
        }
    }

    set_conc_extrema();
}

/** Sets up the odor concentration and reward field (Rosenbrock). */
void navigate_theo::init_conc_and_reward2() {
    double *cp=c1,*rp=r1;
    for(int j=0;j<n;j++) {
        double y=j;
        for(int i=0;i<m;i++,rp++,cp++) {
            double x=i,fx=-0.4*x+3,fy=-0.01*x*x-0.2*y+10;
            // add random noise
            // *cp=1000*exp(-0.01*fx*fx-0.3*fy*fy)+gsl_ran_gaussian(rng,0.001);
            *cp=1000*exp(-0.01*fx*fx-0.3*fy*fy);

            *rp=fx*fx+30*fy*fy<10?1:0;
        }
    }
    set_conc_extrema();
}

/** Sets up the odor concentration and reward field (zig-zag). */
void navigate_theo::init_conc_and_reward3() {
    double *cp=c1,*rp=r1;
    for(int j=0;j<n;j++) {
        double y=j;
        for(int i=0;i<m;i++,rp++,cp++) {
            double x=i,fx=0.5-x/50,fy=0.2*sin(0.3*y);
            *cp=y*exp(-10*(fx+fy)*(fx+fy));
            *rp=y*exp(-10*(fx+fy)*(fx+fy))>40?1:0;
        }
    }

    set_conc_extrema();
}

/** Sets up the odor concentration and reward field (experiment) - 3 ports. */
void navigate_theo::init_conc_and_reward4() {
    double *f1=read_cfile("nt1.odr/3d_cfinal1.mat"); //read the file at first port
    double *f2=read_cfile("nt1.odr/3d_cfinal2.mat"); //read the file at second port
    double *f3=read_cfile("nt1.odr/3d_cfinal3.mat"); //read the file at third port

    // concentration data in cfinal 1
    double *cp1=c1,*rp1=r1;
    for(int j=0;j<n;j++) {
        for(int i=0;i<m;i++) {
            for(int k=0;k<n_frame;k++,rp1++,cp1++){
                 *cp1=f1[k+i*n_frame+j*m*n_frame];
                 *rp1=f1[n_frame-1+i*n_frame+j*m*n_frame]>300?1:0;
            }
        }
    }

    // concentration data in cfinal 2
    double *cp2=c2,*rp2=r2;
    for(int j=0;j<n;j++) {
        for(int i=0;i<m;i++) {
            for(int k=0;k<n_frame;k++,rp2++,cp2++){
                 *cp2=f2[k+i*n_frame+j*m*n_frame];
                 *rp2=f2[n_frame-1+i*n_frame+j*m*n_frame]>300?1:0;
            }
        }
    }

    // concentration data in cfinal 3
    double *cp3=c3,*rp3=r3;
    for(int j=0;j<n;j++) {
        for(int i=0;i<m;i++) {
            for(int k=0;k<n_frame;k++,rp3++,cp3++){
                 *cp3=f3[k+i*n_frame+j*m*n_frame];
                 *rp3=f3[n_frame-1+i*n_frame+j*m*n_frame]>300?1:0;
            }
        }
    }

    set_conc_extrema();
    puts("Initialization done!");

    delete [] f1;
    delete [] f2;
    delete [] f3;

}

/** Read concentration field data from cfilename.
 * \return The int array containing all the information */
double* navigate_theo::read_cfile(const char* cfilename) {

    // Open the MAT file using the matio library
    /*
    mat_t *matfp;
    matvar_t *matvar;
    matfp=Mat_Open(cfilename,MAT_ACC_RDONLY);
    if(NULL==matfp) {
        fprintf(stderr,"Error opening MAT file %s\n",cfilename);
        exit(1);
    }

    // Perform several basic checks to ensure that there is a 2D array of
    // double precision numbers available in the file
    if((matvar=Mat_VarReadNext(matfp))==NULL) {
        fprintf(stderr,"Error reading data from MAT file %s\n",cfilename);
        exit(1);
    }
    if(matvar->rank!=3) {
        fputs("Data in MAT file is not a 3D array\n",stderr);
        exit(1);
    }

    puts("Start reading file");

    // Begin custom edit to convert to transpose.
    size_t t=matvar->dims[0];
    size_t x=matvar->dims[1];
    size_t y=matvar->dims[2];
    // size_t p=matvar->dims[3];

    size_t txyp=t*x*y;

    int *tmp=new int[txyp];
    // static double f[14400];
    double *f=new double[txyp];
    memcpy(tmp,matvar->data,txyp*sizeof(int));

    // print the nonzero value of matrix
    // Save the data
    for(int k=0;k<t;k++){
        for(int i=0;i<x;i++){
            for(int j=0;j<y;j++){
                f[k*y*x+i*y+j]=tmp[k*y*x+i*y+j];
            }
        }
    }

    puts("Finish reading");
    delete [] tmp;

    Mat_VarFree(matvar);
    if((matvar=Mat_VarReadNextInfo(matfp))!=NULL) {
        fputs("Warning: ignoring extra entries in MAT file\n",stderr);
        Mat_VarFree(matvar);
    }

    Mat_Close(matfp);

    return f;*/
    return NULL;
}

/** Simulates a training epsiode to update the value array.
 * \param[in] (i,j) the starting position.
 * \param[in] steps the number of steps in the epsiode.
 * \param[in] eps the probability of taking a random action.
 * \param[in] fp a file handle to save the training data to. If this is the
 *               null pointer, then no output will be saved.
 * \return The total reward for this episode. */
double navigate_theo::training_episode(int i,int j,int steps,double eps,FILE *fp,int port) {
    double rwd=0,reach_rwd=0;
    int a,new_s,s=state(i,j,0,port);
    double prob_field;
    int count_port1,count_port2,count_port3;
    int odor_level;

    if(fp!=NULL) fprintf(fp,"%d %d %d",i,j,s);

    for(int k=1;k<=steps;k++) {
        // Choose the best action, with a probability of taking a random action
        a=rl.pick_best_with_random(s,eps);

        // a=rl.pick_best_with_prob(s,eps);

        // Move the actor, and calculate the reward from taking this step
        move(i,j,a);

        // change concentration fields

        if(k%500==0){
            // prob_field=gsl_rng_uniform(rng);
            // port=prob_field<1.0/3?1:(prob_field<2.0/3?2:3);
            // port=prob_field<1.0/2?1:2;
            port=port%3+1; // start from initialize port, changes in order 1-2-3-1
        }
        // port=(k/2000)%3+1; // port changes in order 1-2-3

        if(port==1) {rwd=r1[i+j*m]; count_port1++;}
        else if(port==2) {rwd=r2[i+j*m]; count_port2++;}
        else {rwd=r3[i+j*m]; count_port3++;}

        // Record the time reach reward
        if(rwd!=0 and reach_rwd==0) reach_rwd=k;
        if(fp!=NULL) fprintf(fp," %d %g\n",a,rwd);

        // Perform the reinforcement learning update
        new_s=state(i,j,k,port);
        odor_level=strong_odor_cue(i,j,k,port);
        rl.update(s,a,new_s,rwd,k);
        if(fp!=NULL) fprintf(fp,"%d %d %d",i,j,s);

        // Update the total reward and state
        s=new_s;

    }
    if(fp!=NULL) fputs("\n\n\n",fp);
    // printf("port1=%d,port2=%d,port3=%d\n",count_port1,count_port2,count_port3);

    return reach_rwd;
}

/** Simulates a test episode to present final result.
 * \param[in] (i,j) the starting position.
 * \param[in] steps the number of steps in the epsiode.
 * \param[in] eps the probability of taking a random action.
 * \param[in] fp a file handle to save the training data to. If this is the
 *               null pointer, then no output will be saved.
 * \return The total reward for this episode. */
double navigate_theo::test_episode(int i,int j,int steps,double eps,FILE *fp,int port) {
    double rwd=0,reach_rwd=0;
    int a,new_s,s=state(i,j,0,port);
    if(fp!=NULL) fprintf(fp,"%d %d %d",i,j,s);
    for(int k=1;k<=steps;k++) {

        // Choose the best action, with a probability of taking a random action
        a=rl.pick_best_with_random(s,eps);

        // Move the actor, and calculate the reward from taking this step
        move(i,j,a);
        /*
        if(k%100==0){
            double prob_field=gsl_rng_uniform(rng);
            port=prob_field<1.0/3?1:(prob_field<0.667?2:3);
            // port=prob_field<1.0/2?1:2;
        }*/
        port=(k/500)%3+1;

        if(port==1) rwd=r1[i+j*m];
        else if(port==2) rwd=r2[i+j*m];
        else rwd=r3[i+j*m];

        // Record the time reach reward
        if(rwd!=0 and reach_rwd==0) reach_rwd=k;
        if(fp!=NULL) fprintf(fp," %d %g\n",a,rwd);
        new_s=state(i,j,k,port);
        if(fp!=NULL) fprintf(fp,"%d %d %d",i,j,s);

        // Update the total reward and state
        s=new_s;
    }
    if(fp!=NULL) fputs("\n\n\n",fp);
    return reach_rwd;
}

/** Simulates a episode doing gradient descent.
 * \param[in] (i,j) the starting position.
 * \param[in] steps the number of steps in the epsiode.
 * \param[in] eps the probability of taking a random action. */
void navigate_theo::grad_descent(int episodes,int steps,int port) {
    int i,j;
    // Output the file to sotre the gradient descent events
    FILE *fp=NULL;
    sprintf(buf,"%s/grad.dat",filename);
    fp=safe_fopen(buf,"w");

    for(int l=0;l<episodes;l++) {

        // Find a starting location that has no reward
        do {
            i=gsl_rng_uniform_int(rng,m);
            j=gsl_rng_uniform_int(rng,n);
            // i=10;
            // j=40;

        } while(r1[i+m*j]==1);

        // Perform the gradient descent episode
        double rwd=0,trwd=0;
        int a,new_s,s=state(i,j,0,port);
        if(fp!=NULL) fprintf(fp,"%d %d %d",i,j,s);
        for(int k=1;k<=steps;k++) {

            // Choose the action according to gradient
            if(discrete_conc_grad(i,j,k,port)==3){
                a=1;}
            else{
                a=discrete_conc_grad(i,j,k,port);}

            // Move the actor, and calculate the reward from taking this step
            move(i,j,a);
            if(port==1) rwd=r1[i+m*j];
            else if(port==2) rwd=r2[i+m*j];
            else rwd=r3[i+m*j];

            if(fp!=NULL) fprintf(fp," %d %g\n",a,rwd);

            // Perform the reinforcement learning update
            new_s=state(i,j,k,port);
            if(fp!=NULL) fprintf(fp,"%d %d %d",i,j,s);

            // Update the total reward and state
            trwd+=rwd;
            s=new_s;
        }
        if(fp!=NULL) fputs("\n\n\n",fp);
    }
    if(fp!=NULL) fclose(fp);
}

/** Performs the Q-learning algorithm by training.
 * \param[in] episodes the number of training episodes
 * \param[in] steps the number of steps in each episode.
 * \param[in] eps the probability of taking a random step during the training.
 * \param[in] freq the frequency with which to store the training events. */
void navigate_theo::train(int episodes,int steps,double eps,int freq) {
    int i,j,port;
    int count_mode1=0;
    int count_mode2=0;
    int count_mode3=0;
    double prob_field;

    // Output the file to sotre the training events
    FILE *fp=NULL;
    if(freq>0) {
        sprintf(buf,"%s/train.dat",filename);
        fp=safe_fopen(buf,"w");
    }

    for(int l=0;l<episodes;l++) {

        // change the concentration field
        prob_field=gsl_rng_uniform(rng);
        port=prob_field<1.0/3?1:(prob_field<2.0/3?2:3);
        // port=prob_field<1.0/2?1:2;
        // port=3;

        // printf("the mode is %d, prob is %f in train function\n",port,prob_field);

        if(port==1) {count_mode1++;}
        else if(port==2) {count_mode2++;}
        else if(port==3) {count_mode3++;}

        // Find a starting location that has no reward

        i=gsl_rng_uniform_int(rng,n_block);
        j=gsl_rng_uniform_int(rng,n_block);
        i*=block_size;
        j*=block_size;
            // i = 60;
            // j = 80;

        // Perform the training episode, and store the data if needed
        // reach_rwd[l%freq]=training_episode(i,j,steps,eps,freq>0&&l%freq==0?fp:NULL);
        training_episode(i,j,steps,eps,freq>0&&l%freq==0?fp:NULL,port);

        // save q values in each episode
        if(l%freq==0){
            sprintf(buf,"nt1.odr/Q_value/Q_value_%d.bin",l/freq);
            write_Q(buf);
            sprintf(buf,"nt1.odr/Freq_value/F_value_%d.bin",l/freq);
            write_Freq(buf);
        }

    }
    printf("mode1 = %d, mode2=%d, mode3=%d\n",count_mode1,count_mode2,count_mode3);

    // Perform test after training

    for(int l=0;l<500;l++) {

        // Find a starting location that has no reward
        do {
            // i=gsl_rng_uniform_int(rng,m);
            // j=gsl_rng_uniform_int(rng,n);
            i=gsl_rng_uniform_int(rng,n_block);
            j=gsl_rng_uniform_int(rng,n_block);
            i*=block_size;
            j*=block_size;
            // i = 60;
            // j = 60;

        } while(r1[i+m*j]==1);

        // Perform the test episode, and store the data
      test_episode(i,j,steps,0,fp,1);

    }
    if(fp!=NULL) fclose(fp);
}

/** Calculates the state associated with a given position and time.
 * \param[in] (i,j) the current position.
 * \param[in] k the number of steps (time).
 * \return The state. */
int navigate_theo::state(int i,int j,int k,int port) {
    // return discrete_conc_grad(i,j,k,port); // return odor gradient
    // return k; // return time state
    // return i/block_size+n_block*(j/block_size); // return  n_block*n_block space state
    // return i+m*j; // return space state
    // return discrete_conc_grad(i,j,k,port)+5*strong_odor_cue(i,j,k,port); // return odor and time
    return strong_odor_cue(i,j,k,port)+5*(i/block_size+n_block*(j/block_size)); //return odor and space
    // return strong_odor_cue(i,j,k,port);
    // return weak_odor_cue(i,j,k,port)+2*(i/block_size+n_block*(j/block_size));

}

/** Calculates a state from 0 to 4 based on the concentration gradient
 * \param[in] (i,j) the current position.
 * \param[in] port is the current concentration field
 * \return The integer state. */
int navigate_theo::discrete_conc_grad(int i,int j,int k,int port) {

    // Calculate the unnormalized concentration gradient
    double *cp;
    if(port==1) {cp=c1+(i+j*m); }
    else if(port==2) {cp=c2+(i+j*m);}
    else {cp=c3+(i+j*m);}
           // cy0=j-1>0?(j-1<n-1?0.5*(cp[m-1]-cp[-m-1]):cp[-1]-cp[-m-1]):cp[m-1]-cp[-1],
           // cy1=j>0?(j<n-1?0.5*(cp[m]-cp[-m]):*cp-cp[-m]):cp[m]-*cp,
           // cy2=j+1>0?(j+1<n-1?0.5*(cp[m+1]-cp[-m+1]):cp[1]-cp[-m+1]):cp[m+1]-cp[1],
           // cy=0.25*(cy0+2*cy1+cy2),
           // // cy=0.25*(cy0+2*cy1+cy2)+gsl_ran_gaussian(rng,1), //add some noise
           // cx0=i-1>0?(i-1<m-1?0.5*(cp[-m+1]-cp[-m-1]):cp[-m]-cp[-m-1]):cp[-m+1]-cp[-m],
           // cx1=i>0?(i<m-1?0.5*(cp[1]-cp[-1]):*cp-cp[-1]):cp[1]-*cp,
           // cx2=i+1>0?(i+1<m-1?0.5*(cp[m+1]-cp[m-1]):cp[m]-cp[m-1]):cp[m+1]-cp[m],
           // cx=0.25*(cx0+2*cx1+cx2);
           // // cx=0.25*(cx0+2*cx1+cx2)+gsl_ran_gaussian(rng,1); // add some noise

    // add noise
    // int rand_gauss_factor1=gsl_ran_gaussian(rng,2);
    // double   cy=j>0?(j<n-1?0.5*(cp[m]-cp[-m])+rand_gauss_factor1:*cp-cp[-m]+rand_gauss_factor1):cp[m]-*cp+rand_gauss_factor1;
    // int rand_gauss_factor2=-gsl_ran_gaussian(rng,2);
    // double cx=i>0?(i<m-1?0.5*(cp[1]-cp[-1])+rand_gauss_factor2:*cp-cp[-1]+rand_gauss_factor2):cp[1]-*cp+rand_gauss_factor2;

    // no noise
    double cy=j>0?(j<n-1?0.5*(cp[m]-cp[-m]):*cp-cp[-m]):cp[m]-*cp;
    double cx=i>0?(i<m-1?0.5*(cp[1]-cp[-1]):*cp-cp[-1]):cp[1]-*cp;

    // Determine the state from the gradient.
    // 0: right, 1: up, 2: left, 3: down, 4: zero.x
    if(cx*cx+cy*cy<navigate_zgt_sq*(dx*dx+dy*dy)) return 4;
    return fabs(cx)>fabs(cy)?(cx>0?0:2):(cy>0?1:3);
}

/** Calculates a state from 0 to 1 based on the concentration field at port
 * \param[in] (i,j) the current position.
 * \param[in] port is the current concentration field
 * \return The integer state. */
int navigate_theo::weak_odor_cue(int i,int j,int k,int port) {
    double *cp;
    if(port==1) {cp=c1+(i+j*m); }
    else if(port==2) {cp=c2+(i+j*m);}
    else {cp=c3+(i+j*m);}

    if(*cp>100) return 1;
    else return 0;
}

/** Calculates a state from 0 to 10 based on the concentration field at port
 * \param[in] (i,j) the current position.
 * \param[in] port is the current concentration field
 * \return The integer state. */
int navigate_theo::strong_odor_cue(int i,int j,int k,int port) {
    double *cp;
    int level=5;
    double gap=cmax/level;
    if(port==1) {cp=c1+(i+j*m); }
    else if(port==2) {cp=c2+(i+j*m);}
    else {cp=c3+(i+j*m);}

    if(*cp<0.1) return 0;
    else if(*cp<10) return 1;
    else if(*cp<50) return 2;
    else if(*cp<80) return 3;
    else return 4;
    /*
    for(int i=0;i<level-1;i++){
        if(gap*i<=*cp && *cp<gap*(i+1)) {
            return i;

        }
    }
    return level-1;*/
}

/** Determines if a given action will take the actor out of bounds.
 * \param[in] (i,j) the actor's current position.
 * \param[in] a the proposed action.
 * \return True if out of bounds, false otherwise. */
inline bool navigate_theo::out_of_bounds(int i,int j,int a) {
    switch(a) {
        case 0: return i==m-1;
        case 1: return j==n-1;
        case 2: return i==0;
        case 3: return j==0;
        default: return false;
    }
}

/** Moves the actor according to a given action.
 * \param[in,out] (i,j) the actor's position, to be updated.
 * \param[in] a the action. */
inline void navigate_theo::move(int &i,int &j,int a) {
    int step_size=block_size;
    switch(a) {
        case 0: if(i<m-step_size) i+=step_size;return;
        case 1: if(j<n-step_size) j+=step_size;return;
        case 2: if(i>step_size-1) i-=step_size;return;
        case 3: if(j>step_size-1) j-=step_size;return;
        case 4: return;
    }
}

/** Finds the extremal values of the concentration field. */
void navigate_theo::set_conc_extrema() {
    cmin=cmax=*c1;
    for(double *cp=c1+1;cp<c1+mn;cp++) {
        if(*cp>cmax) cmax=*cp;
        else if(*cp<cmin) cmin=*cp;
    }
}

/** Writes the concentration and reward fields. */
void navigate_theo::write_fields() {
    sprintf(buf,"%s/conc1.fld",filename);
    output_field(0);
    sprintf(buf,"%s/rwd1.fld",filename);
    output_field(1);

}

/** Outputs a field to a file in a format that can be read by Gnuplot.
 * \param[in] mode the field to output (0: concentration, 1: reward). */
void navigate_theo::output_field(const int mode) {

    // Open the output file, using the filename previously stored
    // in the temporary buffer
    FILE *outf=safe_fopen(buf,"wb");

    // Output the first line of the file
    int i,j;
    float *buf=new float[m+1],*bp=buf+1,*be=bp+m;
    *buf=m;
    for(i=0;i<m;i++) *(bp++)=(i+0.5)*dx;
    fwrite(buf,sizeof(float),m+1,outf);

    // Output the field values to the file

    double *fp=mode==0?c1:r1;
    for(j=0;j<n;j++) {
        *buf=(j+0.5)*dy;bp=buf+1;
        while(bp<be) *(bp++)=*(fp++);
        fwrite(buf,sizeof(float),m+1,outf);
    }

    // Close the file and free temporary buffer
    fclose(outf);
    delete [] buf;
}
