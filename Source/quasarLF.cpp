/*
 *  Created on: 03/12/2013
 *       Author: F. Bellagamba
*/

#include "standard.h"
#include "source.h"

QuasarLF::QuasarLF
	(double my_red                      // redshift
	, double my_mag_limit				// magnitude limit
	, long *my_seed):
		red(my_red), mag_limit(my_mag_limit), seed(my_seed)
{
	// the one used in the paper(s)
	COSMOLOGY cosmo(0.3,0.7,0.7,-1.);

	// k-correction for i band and quasar SED in the redshift range [0,5.5]
	// We should check that the input is in a desired range (see also fitting models for LF below)
	std::ifstream file_in("kcorr_Richards06.txt");
	double red_arr[550];
	double kcorr_arr[550];
	for (int i = 0; i < 550; i++)
	{
		file_in >> red_arr[i] >> kcorr_arr[i];
		if (fabs(red_arr[i]-red)<=0.005) kcorr = kcorr_arr[i];
	}

	// QLF described as double power law

	// PLE model: good fit at redshift 0.3 - 2.2
	if (red < 2.2)
	{
		mstar = -22.85 - 2.5*(1.241*red -0.249*red*red);
		log_phi = -5.96;
		alpha = -1.16;
		beta = -3.37;
	}
	// LEDE model: good fit at redshift 2.2 - 3.5
	else
	{
		mstar = -26.57 - 0.809*(red-2.2);
		log_phi = -5.93 - 0.689*(red-2.2);
		alpha = -1.29;
		beta = -3.51;
	}
	// limits of the integration in absolute magnitudes
	double mag_min, mag_max;
	// luminosity distance
	dl = cosmo.lumDist(red);
	// conversion to absolute magnitude
	mag_max = mag_limit - 5*log10(dl*1.e+05) - kcorr;
	mag_min = 10 - 5*log10(dl*1.e+05) - kcorr;

	arr_nbin = 10000;
	mag_arr = new double[arr_nbin];
	lf_arr = new double[arr_nbin];

	// integrate the function to get normalisation for P(m)
	norm = nintegrateQLF(&QuasarLF::lf_kernel, mag_min, mag_max, 0.001);

	// saves array of M_i, P(m_i) for future random extraction
	for (int i = 0; i < arr_nbin; i++)
	{
		mag_arr[i] = mag_min + (mag_max-mag_min)/(arr_nbin-1)*i;
		lf_arr[i] = nintegrateQLF(&QuasarLF::lf_kernel, mag_min, mag_arr[i], 0.001)/norm;
	}

}

/// returns random apparent magnitude according to the luminosity function
double QuasarLF::getRandomMag()
{
	// extracts random number r between [0,1]
	// m:P(m) = r is the desired random magnitude
	double r = ran2(seed);
	bool found = false;
	int k;
	double p;

	for (int i = 0; i < arr_nbin && found == false; i++)
	{
		if (lf_arr[i] > r)
		{
			k = i-1;
			p = (r-lf_arr[i-1])/(lf_arr[i]-lf_arr[i-1]);
			found = true;
		}
	}

	// interpolates
	double mag_out = mag_arr[k] + p*(mag_arr[k+1]-mag_arr[k]);
	// converts back to apparent magnitude
	mag_out += 5*log10(dl*1.e+05) + kcorr;

	return mag_out;
}

/// returns random flux according to the luminosity function
/// must be divided by an angular area in rad^2 to be a SurfaceBrightness for the ray-tracer
double QuasarLF::getRandomFlux()
{
	double mag = getRandomMag();
	return pow(10, -0.4*(mag+48.6))*inv_hplanck;
}

double QuasarLF::nintegrateQLF(pt2MemFunc func, double a,double b,double tols) const
{
	int jmax = 34;
	int jmaxp = 35;
	int k = 6;
   double ss,dss;
   double s2[jmaxp],h2[jmaxp+1];
   int j;
   double ss2=0;

   h2[1]=1.0;
   for (j=1;j<=jmax;j++) {
     s2[j]=trapzQLFlocal(func,a,b,j,&ss2);
	if (j>=k) {
	   polintD(&h2[j-k],&s2[j-k],k,0.0,&ss,&dss);
	   if(fabs(dss) <= tols*fabs(ss)) return ss;
	}
	h2[j+1]=0.25*h2[j];
	}
   std::cout << "s2= "; for(j=1;j<=jmax;j++) std::cout << s2[j] << " ";
   std::cout << "\n";
   std::cout << "Too many steps in routine nintegrateQLF\n";
   return 0.0;
}

double QuasarLF::trapzQLFlocal(pt2MemFunc func, double a, double b, int n, double *s2) const
{
  double x,tnm,del,sum;
   int it,j;

   if (n == 1) {
	return (*s2=0.5*(b-a)*( (this->*func)(a) +(this->*func)(b) ));
   } else {

	for (it=1,j=1;j<n-1;j++) it <<= 1;
	tnm=it;
	del=(b-a)/tnm;
	x=a+0.5*del;
	for (sum=0.0,j=1;j<=it;j++,x+=del) sum += (this->*func)(x);
	*s2=0.5*(*s2+(b-a)*sum/tnm);

	return *s2;
   }
}

double QuasarLF::lf_kernel (double mag) const
{
	double den = pow(10,0.4*(alpha+1)*(mag-mstar)) + pow(10,0.4*(beta+1)*(mag-mstar));
	return 1./den;
}

