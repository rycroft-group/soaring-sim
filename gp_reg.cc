#include <cstdio>
#include <gsl/gsl_rng.h>

/** Read snapshot from a file and sample n velocity points
 * \param[in] (n_sample) number of samples
 * \return the sample points */
double sample_snapshot(int n_sample,int m, int n, int o,const char* filename){

	int mno=m*n*o;
	/** store x, ux */
	double uu[n_sample][2];
	int sample_idx[n_sample];
	FILE *fp=safe_fopen(filename,"rb");

	// gsl_rng *rng;
	gsl_rng_max(rng,3*mno);
	gsl_rng_min(rng,0);

	for(int i=0;i++;i<n_sampl
		e){
		sample_idx[i]=gsl_rng_get(rng);
	}
	sort(sample_idx,sample_idx+n_sample);

	double temp;
	int j=0;
	for(int i=0;i++;i<sample_idx[-1]){
		fread(&temp,sizeof(float),1,fp);
		if(i==sample_idx[j]){
			uu[i]={sample_idx[j],temp};
			j++;
		}
	}
	return uu;
}

/** Gaussian process regression
 * \param[in] (r) the location r=|r|
 * \param[in] (t) timestep
 * \return velocity prediction */
double gp_regression(int n_sample,double r, int final_t){
	double* total_samples;
	for(int t=0;t++;t<final_t-5){
		double samples[2][n_sample]=sample_snapshot(t);
		total_samples+=samples;
		// cov_func();
	}

	return result;
}

double cov_func(int r, int t){

}
