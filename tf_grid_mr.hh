#ifndef TF_GRID_MR_HH
#define TF_GRID_MR_HH

#include "tf_grid.hh"

class turb_fluid_grid_mr : public turb_fluid_grid {
    public:
        /** The number of segments in time for Hermite interpolation. */
        const int h_segs;
        /** The maximum time for Hermite interpolation. */
        const double Tmr;
        /** The scale factor used in the internal coordinate transformation. */
        const double rho;
        /** The inverse factor. */
        const double irho;
        /** A pointer to the gridded velocity information at the segment
         * control points. */
        double** const umr;
        turb_fluid_grid_mr(int m_,int n_,int o_,double ax_,double bx_,double ay_,double by_,double az_,double bz_,double C_,double alpha_,int h_segs_,double Tmr_,unsigned long seed=1);
        ~turb_fluid_grid_mr();
        void calc_mr_fields();
        void lin_interp_mr(double T,double x,double y,double z,double &ux,double &uy,double &uz);
        void cub_interp_mr(double T,double x,double y,double z,double &ux,double &uy,double &uz);
};

#endif
