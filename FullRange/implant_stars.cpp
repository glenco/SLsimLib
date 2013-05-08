/*
 * implant_stars.c
 *
 *  Created on: Jun 5, 2010
 *      Author: R.B. Metcalf
 */

#include "slsimlib.h"
#include <iostream> // remove
#include <sstream> // remove
using namespace std;
/** \ingroup ChangeLens
 * \brief  Implants stars into the lens around the images.
 *
 * lens->Nstars and lens->stars_fstars must be set before calling
 * allocates all memory for stars
 */

void BaseAnaLens::implant_stars(Point *centers,unsigned long Nregions,long *seed, int mftype){
	PosType r,theta,NstarsPerImage;
	unsigned long i,j,m,k;
	//float maxr;

	if(stars_N < 1.0  || star_fstars <= 0) return;
	if(star_fstars > 1.0){ std::printf("fstars > 1.0\n"); exit(0); }
	if(!(stars_implanted) ){

		star_masses = new float[stars_N];
		stars = new unsigned long[stars_N];
		stars_xp = Utilities::PosTypeMatrix(stars_N,3);
		star_theta_force = 1.0e-1;

		assert(Nregions > 0);
		star_Nregions = Nregions;
		star_region = new double[Nregions];
		star_kappa = new double[Nregions];
		star_xdisk = Utilities::PosTypeMatrix(Nregions,2);

	}else{
		// free star_tree
		delete star_tree;
	}

	if(stars_N < 1  || star_fstars <= 0){
		stars_implanted = true;
		stars_N = 0;
		star_fstars = 0.0;

		for(j=0,m=0;j<Nregions;++j){
			star_region[j] = 0.0;
			star_kappa[j] = 0.0;
			star_xdisk[j][0] = centers[j].x[0];
			star_xdisk[j][1] = centers[j].x[1];

		}
		return;
	}

	for(j=0,m=0;j<Nregions;++j){

		assert( centers[j].kappa > 0.0);

		NstarsPerImage = stars_N/star_Nregions;
		//cout << "NstarsPerImage: " << NstarsPerImage << endl; // remove

		star_region[j] = 1.0/sqrt(pi*star_fstars*centers[j].kappa*Sigma_crit
				/star_massscale/(float)NstarsPerImage);


		cout << star_region[j] << " " << star_fstars <<  " " << centers[j].kappa << " " << Sigma_crit << " " << star_massscale << " " << NstarsPerImage << endl;
		//cout << "Sigma_star= " << (float)NstarsPerImage*star_massscale/(pi*pow(star_region[j],2)) << endl;
		// cutoff based on comparison of star deflection to smooth component
		//rcut = 4*sqrt(star_massscale/pi/Sigma_crit
		//		/( centers[j].kappa+sqrt(pow(centers[j].gamma[0],2)+pow(centers[j].gamma[1],2)) ) );

		star_kappa[j] = star_fstars*centers[j].kappa;

		star_xdisk[j][0] = centers[j].x[0];
		star_xdisk[j][1] = centers[j].x[1];

		//printf("kappa = %e  star_region = %e\n",star_kappa[j],star_region[j]);
		char *fname = "stars0.dat";
		if(j==0){fname = "stars0.dat";}
		if(j==1){fname = "stars1.dat";}
		if(j==2){fname = "stars2.dat";}
		if(j==3){fname = "stars3.dat";}
		ofstream fstars(fname);
		ofstream mstars("masses.dat");

		star_masses=stellar_mass_function(mftype, NstarsPerImage, seed);

		for(i=0;i<NstarsPerImage;++i,++m){
			//m=j*NstarsPerImage+i;
			r = star_region[j]*sqrt(ran2(seed));
			theta=2*pi*ran2(seed);
			stars_xp[m][0] = centers[j].x[0] + r*cos(theta);
			stars_xp[m][1] = centers[j].x[1] + r*sin(theta);
			stars_xp[m][2] = 0.0;
			cout << m << " " << star_masses[m] << endl;

			//if(maxr<r){
			//	maxr=r;
			//}
			// check to make see if star is in another centers region
			//   and remove it
			for(k=0;k<j;++k){
				if(  star_region[k] > sqrt(pow(centers[k].x[0]-stars_xp[m][0],2)
						+ pow(centers[k].x[1]-stars_xp[m][1],2)) ){
					//--NstarsPerImage;
					//--i;
					--m;
					break;
				}
			}
			fstars << scientific << stars_xp[m][0] << " " << stars_xp[m][1] << endl;
			mstars << scientific << star_masses[m] << endl;

			//cout << "max r" << maxr << " " << star_region[j] << endl;
			//printf("%e %e\n",stars_xp[m][0],stars_xp[m][1]);
		}
		fstars.close();
		mstars.close();
	}

	assert(m <= stars_N);
	stars_N = m;
	//std::printf("last star x = %e %e\n",stars_xp[stars_N-1][0],stars_xp[stars_N-1][1]);

	float dummy=0;
	//star_tree = new ForceTree(stars_xp,stars_N,star_masses,&dummy
	//		,false,false,5,2,false,star_theta_force);

	star_tree = new QuadTree(stars_xp,star_masses,&dummy,stars_N
			,false,false,0,4,star_theta_force);

	// visit every branch to find center of mass and cutoff scale */
	stars_implanted = true;

	return ;
}

// This allows the stars to be turned off after they have been implanted.
/*void AnaLens::toggleStars(bool implanted){
	stars_implanted = implanted;
}
*/

/// subtracts the mass in stars from the smooth model to compensate
/// for the mass of the stars the lensing quantities are all updated not replaced
void BaseAnaLens::substract_stars_disks(double *ray,double *alpha
		,KappaType *kappa,KappaType *gamma){

	if(!(stars_implanted)) return;

	double xcm,ycm,r;
	float tmp;
	double mass;
	int i;

	for(i=0;i<star_Nregions;++i){
		xcm = star_xdisk[i][0] - ray[0];
		ycm = star_xdisk[i][1] - ray[1];
		r=sqrt(xcm*xcm + ycm*ycm);

		if(r < star_region[i]){
			alpha[0] += star_kappa[i]*xcm;
			alpha[1] += star_kappa[i]*ycm;
			*kappa -= star_kappa[i];
		}else{

			mass = star_kappa[i]*pow(star_region[i],2)/r/r;
			alpha[0] += mass*xcm;
			alpha[1] += mass*ycm;

			tmp = 2*mass/r/r;
			gamma[0] += 0.5*(xcm*xcm-ycm*ycm)*tmp;
			gamma[1] += xcm*ycm*tmp;
		}
	}

	return;
}

// random stellar masses according to IMF of choice
/* mtype defines the stellar mass function
 * 0 - always the same stellar mass (e.g. 1Msol)
 * 1 - salpeter imf, i.e. slope = -2.35
 * 2 - broken power law, requires lower mass end slope (powerlo), high mass slope (powerhi), bending point (bendmass)
 * 3 - further IMF models may follow
 */
float* BaseAnaLens::stellar_mass_function(int mftype, unsigned long Nstars, long *seed, double minmass, double maxmass, double bendmass, double powerlo, double powerhi){
	//if(!(stars_implanted)) return;
	unsigned long i;
	double powerp1,n0,n1;
	float *star_masses = new float[Nstars];

	if(mftype==0){
		for(i = 0; i < Nstars; i++){
			star_masses[i]=1.0;
		}
		//std::fill_n(star_masses, Nstars, 1.0);

	}

    if(mftype==1){
    	powerp1 = powerhi+1.0;
    	n0 = (pow(maxmass,powerp1)) /powerp1;
    	n1 =  n0 - (pow(minmass,powerp1)) / powerp1;
    	for(i = 0; i < Nstars; i++){
    		star_masses[i] = pow( (-powerp1*(n1*ran2(seed)-n0)),(1.0/powerp1) );
    	//cout << powerp1 << " " << n0 << " " << n1 << endl;
    	}
	}

    /*

    bendmass=1.
    powerl=0.
    powerlp1=powerl+1.
    dummax=maxmass/bendmass
    dummin=minmass/bendmass

    N1=(1./powerlp1)-(dummin**(powerlp1))/(powerlp1)
    N2=(dummax**powerp1)/powerp1-(1./powerp1)
    N=N1+N2

    if(dum1.lt.(N1/N))then
      RM=((N*dum1)*(powerlp1)+dummin**powerlp1)**(1./powerlp1)
    else
      RM=((N*dum1-N1)*powerp1+1.)**(1./powerp1)
    endif
    massf = RM*bendmass
    */

    return star_masses;
}
