#include "game.hh"
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <gsl/gsl_sf.h>

game::game(int board_size_):board_size(board_size_),rest_move(board_size),board(new int[board_size]){
	reset();
}

game::~game(){
	//delete [] board;
}

game_glider::game_glider(double cL_,double cD_,double wind_val_,int wind_,bool frozen_,turb_fluid_grid &tf_,gpr &g_)
    : game(0), cL(cL_), cD(cD_),wind_val(wind_val_), frozen(frozen_),nframes(10000),  wind(wind_), duration(100),
    dt_pad(0.06), wind_field_x(NULL),wind_field_y(NULL),wind_field_z(NULL), tf(tf_), g(g_),rng(gsl_rng_alloc(gsl_rng_taus2)){
  tf.init_steady_state();
  gsl_rng_set(rng,time(NULL));
  // Create glider and initialize its position
  // init(0.,0.,0.,1.,0.,-1./15);

  // Compute a timestep
  dt=dt_pad*tf.est_max_timestep();
  double sint=duration/nframes;
  int l=static_cast<int>(sint/dt)+1;
  dt=sint/l;
  // dt=0.01;

  if(frozen) dt=0.01; // set dt a constant when tf is frozen
 //printf("dt %g\n",dt);
  tau_min=pow(sqrt(3.)*M_PI,(-2./3))/(3.094*tf.Cinv);
  tau_max=pow(sqrt(3.)*M_PI*tf.m,(-2./3))/(3.094*tf.Cinv);
}

game_glider::~game_glider(){
  delete [] wind_field_z;
  delete [] wind_field_y;
	delete [] wind_field_x;
}

inline const char* game::c(int i){
	return board[i]==0?" ":(board[i]==1?"o":"x");
}

/*************************** GLIDER **************************/

void game_glider::init(double rx_,double ry_,double rz_,double ux_,double uy_,double uz_) {
    rx=rx_;
    ry=ry_;
    rz=rz_;
    ux=ux_;
    uy=uy_;
    uz=uz_;
    vx=ux_;
    vy=uy_;
    vz=uz_;
    bank=0.;
    if(wind==2) init_tf();
}

/** Init a turbulent wind field with only wz nozero. */
void game_glider::init_tf(){
/**
  wind_field_x=new double[256*256];
  wind_field_y=new double[256*256];
  wind_field_z=new double[256*256];

  FILE *fp;
  float* temp=new float[256*256];
  char buf[50];
  for(int k=0;k<1;k++){
    // Load wx
    sprintf(buf,"tf.odr/ux.%d",k);
    fp=fopen(buf,"rb");
    fread(temp,sizeof(float),256*256,fp);
    for(int i=0;i<256*256;i++) wind_field_x[k*256*256+i]=temp[i];
    fclose(fp);

    // Load wy
    sprintf(buf,"tf.odr/uy.%d",k);
    fp=fopen(buf,"rb");
    fread(temp,sizeof(float),256*256,fp);
    for(int i=0;i<256*256;i++) wind_field_y[k*256*256+i]=temp[i];
    fclose(fp);

    // Load wz
    sprintf(buf,"tf.odr/uz.%d",k);
    fp=fopen(buf,"rb");
    fread(temp,sizeof(float),256*256,fp);
    for(int i=0;i<256*256;i++) wind_field_z[k*256*256+i]=temp[i]      fclose(fp);
  }

  delete [] temp;
*/
  }
bool game_glider::valid_move(int k){
	return true;
}

void game_glider::remove(int k){}

/** Update velocity and position given bank angle.
 * \param[in] k index to bank angle.
 * \param[in] p player id.
 * \param[in] time the current time frame. */
double game_glider::play(int k,int p,int time) {
    double mu=mu_f(k);
    double t=dt*time;

    double wx0=wx,wy0=wy,wz0=wz;
    // Get wind position at current location

    switch(wind){
      case 0:
        tf.vel(rx,ry,rz,wx,wy,wz);
        // if(time%50==0) g.add_measurement(rx/50.,ry/50.,rz/50.,0.,wx,wy,wz);
        wx*=wind_val,wy*=wind_val,wz*=wind_val;
        break;
      case 1:
        get_wind(rx,ry,rz,wx,wy,wz);
        break;
      case 2:
        get_tf(time,rx,ry,rz,wx,wy,wz);
        wx*=wind_val,wy*=wind_val,wz*=wind_val;
        break;
      case 3:
	// the wind is already updated in the main function
        // if(time%50==0) g.add_measurement(rx/50.,ry/50.,rz/50.,0.,wx,wy,wz);
        // wx*=0.05,wy*=0.05,wz*=0.05;
        break;
    }

    // wx=0,wy=0,wz=0;
    //wx=gsl_ran_gaussian(rng,1.);
    //wy=gsl_ran_gaussian(rng,1.);
    //wz=gsl_ran_gaussian(rng,1.);
    // Update velocity
    // Update position
    rx+=ux*dt;
    ry+=uy*dt;
    rz+=uz*dt;

    // Calculate force
    du(bank,vx,vy,vz);
    ux+=dt*dux;
    uy+=dt*duy;
    uz+=dt*duz;

    //update v
    vx=ux-wx;
    vy=uy-wy;
    vz=uz-wz;

    double v=sqrt(vx*vx+vy*vy+vz*vz);

    // Integrate the velocity field
//    if(!frozen) tf.step_forward(dt);

    double dwx=wx-wx0,dwy=wy-wy0,dwz=wz-wz0;
    // return energy gain
    return -cD*v*v*v+wz-vx*dwx-vy*dwy-vz*dwz;

}

void game_glider::play_multi(int ngl,int* k,int* time,double* r,double* u,double* v,double* duu){
 /**
    // Calculate the velocities at the random sample points
    double* wind=new double[3*ngl];
    tf.allocate_vel_table(ngl);
    tf.vel_multi(ngl,r,wind);

    // TODO
    //if(time%50==0) g.add_measurement(rx/50.,ry/50.,rz/50.,0.,wx,wy,wz);

    double mu,t;
    for(int i=0;i<ngl;i++){
      mu=mu_f(k[i]);
      t=dt*time[i];

      r[3*i]+=u[3*i]*dt;
      r[3*i+1]+=u[3*i+1]*dt;
      r[3*i+2]+=u[3*i+2]*dt;

      du(mu,v[3*i],v[3*i+1],v[3*i+2]);
      u[3*i]+=dt*duu[3*i];
      u[3*i+1]+=dt*duu[3*i+1];
      u[3*i+2]+=dt*duu[3*i+2];

      v[3*i]=u[3*i]-wind[3*i];
      v[3*i+1]=u[3*i+1]-wind[3*i+1];
      v[3*i+2]=u[3*i+2]-wind[3*i+2];

    }*/
}

/** Update velocity and position given bank angle with estimated wind.
 * \param[in] k index to bank angle.
 * \param[in] p player id.
 * \param[in] time the current time frame.
 * return the integral part of the energy gain. */
double game_glider::simulate_play(int k,int p,int time) {
    double mu=mu_f(k);

    // Estimate wind position at current location
    double t=0.;
    if(!frozen) t=time*dt;
    //printf("rx is %f, ry is %f, rz is %f, t is %f ux %g vx %g\n",rx,ry,rz,t,vx,vy);

    switch(wind) {
        case 0:
            g.predict(rx/tf.bx,ry/tf.by,rz/tf.bz,t,wx,wy,wz);
            // wx*=0.05,wy*=0.05,wz*=0.05;
            // Bound for predicted wind velocit
	           // wx=0.,wy=0.,wz=0.;
             //wx+=wind_val*gsl_ran_gaussian(rng,1.);
             //wy+=wind_val*gsl_ran_gaussian(rng,1.);
             //wz+=wind_val*gsl_ran_gaussian(rng,1.);
            if(wx>1.) wx=1.*wind_val;else if(wx<-1.) wx=-1.*wind_val;
            if(wy>1.) wy=1.*wind_val;else if(wy<-1.) wy=-1.*wind_val;
            if(wz>1.) wz=1.*wind_val;else if(wz<-1.) wz=-1.*wind_val;
            break;
        case 1:
            get_wind(rx,ry,rz,wx,wy,wz);
            break;
        case 2:
            get_tf(time,rx,ry,rz,wx,wy,wz);
            wx*=wind_val,wy*=wind_val,wz*=wind_val;
            break;

    }
   // printf("get predict wx is %f, wy is %f, wz is %f\n",wx,wy,wz);

    // Update position
    rx+=ux*dt;
    ry+=uy*dt;
    rz+=uz*dt;

    // Calculate force
    du(bank,vx,vy,vz);
//printf("dux %g duy %g\n",dux,duy);
    // Update velocity
    ux+=dt*dux;
    uy+=dt*duy;
    uz+=dt*duz;

    // Update air velocity
    vx=ux-wx;
    vy=uy-wy;
    vz=uz-wz;

    double v=sqrt(vx*vx+vy*vy+vz*vz);

    // Return energy gain
    double total_gain=-cD*v*v*v+dux*wx+duy*wy+duz*wz+wz;
    if(std::isinf(total_gain)||std::isnan(total_gain)) return 0.;
    else return total_gain;

}

/** Pick an action using fixed policy described in soaring papaer.
 *  Return action index.
 * \param[in] time the current time frame. */
int game_glider::pick_action_policy(int time) {

    double t=dt*time,wz0=wz;
//    double wx0=wx,wy0=wy,wz0=wz;

    // Get wind position at current location
    switch(wind) {
        case 0:
	         tf.vel(rx,ry,rz,wx,wy,wz);
            wx*=wind_val,wy*=wind_val,wz*=wind_val;
            // Bound for predicted wind velocity
            if(wx>1.) wx=1.;else if(wx<-1.) wx=-1.;
            if(wy>1.) wy=1.;else if(wy<-1.) wy=-1.;
            if(wz>1.) wz=1.;else if(wz<-1.) wz=-1.;
            break;
        case 1:
            get_wind(rx,ry,rz,wx,wy,wz);
            break;
        case 2:
            get_tf(time,rx,ry,rz,wx,wy,wz);
            wx*=wind_val,wy*=wind_val,wz*=wind_val;
            break;
    }

    double az=(wz-wz0)/dt;

    // Get left wing and right wing position
    double dx,dy,dz;
    double phi=atan2(vy,vx);
    dx=sin(bank)*sin(phi);
    dy=sin(bank)*cos(phi);
    dz=cos(bank);
    double l=1.;

    // Calculate torque
    double wxl,wyl,wzl,wxr,wyr,wzr;
    //tf.vel(rx-dx*l/2,ry-dy*l/2,rz-dz*l/2,wxl,wyl,wzl);
    //tf.vel(rx+dx*l/2,ry+dy*l/2,rz+dz*l/2,wxr,wyr,wzr);
    get_wind(rx-dx*l/2,ry-dy*l/2,rz-dz*l/2,wxl,wyl,wzl);
    get_wind(rx+dx*l/2,ry+dy*l/2,rz+dz*l/2,wxr,wyr,wzr);
    double tau=(wzl-wzr)*l;

    // 0: negative high, 1: low value, 2: positive high
    az=az>-0.05?(az<0.05?1:2):0;
    tau=tau>-1.?(tau<1.?1:2):0;

    // Pick an action
    // Use the fixed policy in soaring paper, with low ^urms assumption (Fig. 4A)
    double d5=M_PI/18.; // 5 degree
    if(tau==0){
      if(az==0){
        if(bank>=-d5) return 2;
        else if(bank>=-2*d5) return 0;
        else return 1;
      }
      else if(az==1){
        if(bank>=3*d5) return 0;
        else if(bank<=d5&&bank>=-2*d5) return 2;
        else return 0;
      }
      else{
        if(bank>=2*d5) return 0;
        else if(bank==d5) return 1;
        else return 2;
      }
    }

    if(tau==1){
      if(az==0){
        if(bank>=0) return 2;
        else if(bank>=-2*d5) return 0;
        else return 0;
      }
      else if(az==1){
        if(bank>=3*d5) return 0;
        else if(bank>=0||bank<=-3*d5) return 2;
        else return 0;
      }
      else{
        if(bank>=d5) return 0;
        else return 2;
      }

    }

    if(tau==2){
      if(az==0){
        if(bank>=3*d5) return 1;
        else if(bank==2*d5) return 2;
        else return 0;
      }
      else if(az==1){
        if(bank>=3*d5||bank==-2*d5) return 1;
        else if(bank>=-d5) return 0;
        else return 2;

      }
      else{
        if(bank>=0) return 0;
        else if(bank==-d5) return 1;
        else return 2;
      }
    }

  return 1;

}

void game_glider::get_wind(double rx_,double ry_,double rz_,double &wx_,double &wy_,double &wz_){
	if(rx_>15&&rx_<45&&ry_>15&&ry_<45){
             double fac=(rx_-30)*(rx_-30)+(ry_-30)*(ry_-30);
             wz_=0.5*exp(-0.1*sqrt(fac));
     }
     else{
             wz_=0.;
     }
     wx_=0.,wy_=0.;
}

void game_glider::get_tf(int time,double rx_,double ry_,double rz_,double &wx_,double &wy_,double &wz_){
  tf.lin_interp(rx_,ry_,rz_,wx_,wy_,wz_);
  //printf(" rx %g ry %g rz %g wx %g wy %g  wz %g\n",rx_,ry_,rz_,wx_,wy_,wz_);
  // Add gpr predict step for future time frames
  if(!frozen){
	  double ddt=abs((time%50+1)*dt);
	  //printf("ddt %g wx %g wy %g  wz %g\n",ddt,wx_,wy_,wz_);
	  //double kt=ddt/(tau_min-tau_max)*(gsl_sf_gamma_inc(-1.,ddt/tau_min)-gsl_sf_gamma_inc(-1.,ddt/tau_max));
	  gsl_sf_result gamma0;
	  gsl_sf_result gamma1;
	  gsl_sf_gamma_inc_e(-1.,ddt/tau_min,&gamma0);
	  gsl_sf_gamma_inc_e(-1.,ddt/tau_max,&gamma1);
	  double kt=ddt/(tau_min-tau_max)*(gamma0.val-gamma1.val);
      printf("t kt %g %g %g %d %g\n",rx_,ry_,rz_,time,kt);
	  wx_*=kt,wy_*=kt,wz_*=kt;
  }
}

/** Return bank angle according to the action index.
 * \param[in] k index to bank angle. */
double game_glider::mu_f(int k){

  // Check the current bank angle and impose restriction on glider

  switch(k){
    case 0:
      if(bank<=-4*M_PI/18.) return 0.;
      else return -M_PI/18.;
      break;
    case 1:
      return 0.;
      break;
    case 2:
      if(bank>=4*M_PI/18.) return 0.;
      else return M_PI/18.;
      break;
    // case 3:
    //   return M_PI/36.;
    // case 4:
    //   return -M_PI/36.;
  }
  return 0.;

}

void game_glider::du(double mu,double vx_,double vy_,double vz_){
    // Caculate related parameters
    double v2,v,gamma,phi;
    v2=vx_*vx_+vy_*vy_+vz_*vz_;
    v=sqrt(v2);
    gamma=asin(-vz/v);
    phi=atan2(vy,vx);
//   if(std::isinf(gamma)){printf("gamma %g vz %g v %g mu %g\n",gamma,vz,v,mu);}
//   if(std::isnan(phi)){printf("phi %g vy %g vx %g mu %g\n",phi,vy,vx,mu);}
    double cos_mu=cos(mu),sin_mu=sin(mu);
    double cos_gam=cos(gamma),sin_gam=sin(gamma);
    double cos_phi=cos(phi),sin_phi=sin(phi);
//  printf("v2 %g mu %g wx %g vx %g\n", v2,mu,wx,vx);
//printf("wx %g ux %g\n",wx,ux);
    dux=cL*v2*cos_mu*sin_gam*cos_phi-cL*v2*sin_mu*sin_phi-cD*v2*cos_gam*cos_phi;
    duy=cL*v2*cos_mu*sin_gam*sin_phi+cL*v2*sin_mu*cos_phi-cD*v2*cos_gam*sin_phi;
    duz=cL*v2*cos_mu*cos_gam+cD*v2*sin_gam-1.;

}

void game_glider::reset(double rx_,double ry_,double rz_,double ux_,double uy_,double uz_,double vx_,double vy_,double vz_,double bank_){
  rx=rx_,ry=ry_,rz=rz_;
  ux=ux_,uy=uy_,uz=uz_;
  vx=vx_,vy=vy_,vz=vz_;
  bank=bank_;
}

void game_glider::output(double time,int k,FILE *fp) {
    double mu=mu_f(k);
    du(bank,vx,vy,vz);
    fprintf(fp,"%g %g %g %g %g %g %g %g %g %g %g %g %g %d %g\n",time*dt,rx,ry,rz,ux,uy,uz,vx,vy,vz,dux,duy,duz,k,bank);

}
