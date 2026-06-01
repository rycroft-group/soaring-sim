#ifndef GAME_HH
#define GAME_HH

#include <cstdio>
#include "tf_grid.hh"
#include "gpr.hh"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>

class game {
    public:
        /** The size of the board. Size of state space */
        const int board_size;
        /** The number of left moves in the game. */
        int rest_move;
        /** The board of the game. More general: state space */
        int* board;
        game(int board_size_);
        ~game();
        inline const char* c(int i);
        inline void reset(){
            for(int i=0;i<board_size;i++) board[i]=0;
            rest_move=board_size;
        }
        virtual void print()=0;
        virtual bool valid_move(int k)=0;
        virtual double play(int k,int p,int time)=0;
        virtual void remove(int k)=0;
        virtual bool full()=0;
};

class game_glider: public game{
    public:
        /** The lift coefficient. */
        const double cL;
        /** The drag coefficient. */
        const double cD;
        /** The wind vel/glider vel ratio. */
        const double wind_val;
        /** The total number of output points. */
        const int nframes;
        /** The duration to integrate over. */
        const double duration;
        /** The padding factor to apply to the timestep. */
        const double dt_pad;
        /** The wind field. */
        double* wind_field_x;
        double* wind_field_y;
        double* wind_field_z;
        /** The timestep. */
        double dt;
        /** The glider x position. */
        double rx;
        /** The glider y position. */
        double ry;
        /** The glider z position. */
        double rz;
        /** The glider x velocity. */
        double ux;
        /** The glider y velocity. */
        double uy;
        /** The glider z velocity. */
        double uz;
        /** The glider x air velocity. */
        double vx;
        /** The glider y air velocity. */
        double vy;
        /** The glider z air velocity. */
        double vz;
        /** The aero force. */
        double dux;
        double duy;
        double duz;
        /** The wind velocity. */
        double wx;
        double wy;
        double wz;
	double tau_min;
	double tau_max;
        /** The current bank angle. */
        double bank;
        /** The wind field. 0: gpr wind. 1: fixed wind. 2: actual wind value. */
        int wind;
        /** Whether the wind field is static. */
        bool frozen;
        // it doesn't need a board?
        game_glider(double cL_,double cD_,double wind_val_,int wind_,bool frozen_,turb_fluid_grid &tf_,gpr &g_);
        ~game_glider();
        virtual void print() {}
        virtual bool valid_move(int k);
        virtual double play(int k,int p,int time); // equal to step forward function in glider
        virtual void remove(int k);
        virtual bool full() {return false;}
        double simulate_play(int k, int p,int time);
        int pick_action_policy(int time);
        void output(double time,int k,FILE *fp);
        void reset(double rx_,double ry_,double rz_,double ux_,double uy_,double uz_,double vx_,double vy_,double vz_,double bank_);
        double mu_f(int k);
        void get_wind(double rx_,double ry_,double rz_,double &wx_,double &wy_,double &wz_);
        void init(double rx_,double ry_,double rz_,double ux_,double uy_,double uz_);
        void init_tf();
        void play_multi(int ngl,int* k_,int* time_,double* r_,double* u_,double* v_,double* du_);
    private:
        void du(double mu,double vx_,double vy_,double vz_);
        void get_tf(int time,double rx_,double ry_,double rz_,double &wx_,double &wy_,double &wz_);
        /** A reference to a turbulent fluid simulation. */
        turb_fluid_grid &tf;
        /** A reference to gaussian progress regression. */
        gpr &g;
	gsl_rng *rng;
};

#endif
