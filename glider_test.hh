#ifndef GLIDER_TEST_HH
#define GLIDER_TEST_HH

#include "common.hh"
#include "glider.hh"
#include "turb_fluid.hh"

class glider_test {
    public:
        /** The glider model to use. */
        glider_model &gm;
        /** A pointer to the turbulent fluid class for computing the wind. If
         * this is set to zero, then no wind is used. */
        turb_fluid* const tf;
        /** The maximum timestep. */
        const double dt_max;
        /** The test duration. */
        const double duration;
        /** The interval between output frames. */
        const double fr_int;
        /** The initial bank angle. */
        const int bank0;
        /** The number of frames for output. */
        const int nframes;
        glider_test(glider_model &gm_,turb_fluid *tf_,int bank0_,double duration_,int nframes_);
        void integrate(integration_type itype,double dt_pad,double &h,glider &g,FILE *fp=NULL);
        inline void integrate(integration_type itype,double dt_pad,double &h,glider &g,const char* filename) {
            FILE *fp=safe_fopen(filename,"w");
            integrate(itype,dt_pad,h,g,fp);
            fclose(fp);
        }
    private:
        /** Chooses a timestep size that is the largest value smaller than
         * dt_pad*dt_max, such that a given interval length is a perfect
         * multiple of this timestep.
         * \param[in] dt_pad the padding factor applied to the maximum
         *                   timestep.
         * \param[out] dt the timestep size.
         * \return The number of timesteps the fit into the interval. */
        inline int timestep_select(double dt_pad,double &dt) {
            int l=static_cast<int>(fr_int/(dt_pad*dt_max))+1;
            dt=fr_int/l;
            return l;
        }
        /** The initial x velocity. */
        double vx0;
        /** The initial z velocity. */
        double vz0;
        /** Temporary space for computing the wind field. */
        double w[3];
        /** Temporary space for the change in wind field. */
        double dw[3];
};

#endif
