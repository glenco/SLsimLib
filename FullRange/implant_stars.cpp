/*
 * implant_stars.c
 *
 *  Created on: Jun 5, 2010
 *      Author: R.B. Metcalf
 */

#include "slsimlib.h"

/** \ingroup ChangeLens
 * \brief  Implants stars into the lens around the images.
 *
 * lens->Nstars and lens->stars_fstars must be set before calling
 * allocates all memory for stars
 */

void implant_stars(AnaLens *lens,Point *centers,unsigned long Nregions,long *seed){
	PosType r,theta,NstarsPerImage;
	unsigned long i,j,m,k;

	if(lens->stars_N < 1.0  || lens->star_fstars <= 0) return;
	if(lens->star_fstars > 1.0){ std::printf("fstars > 1.0\n"); exit(0); }
	if(!(lens->stars_implanted) ){

		lens->star_masses = new float[lens->stars_N];
		lens->stars = new unsigned long[lens->stars_N];
		lens->stars_xp = PosTypeMatrix(lens->stars_N,3);
		lens->star_theta_force = 1.0e-1;

		assert(Nregions > 0);
		lens->star_Nregions = Nregions;
		lens->star_region = new double[Nregions];
		lens->star_kappa = new double[Nregions];
		lens->star_xdisk = PosTypeMatrix(Nregions,2);

	}else{
		// free star_tree
		delete lens->star_tree;
	}

	if(lens->stars_N < 1  || lens->star_fstars <= 0){
		lens->stars_implanted = true;
		lens->stars_N = 0;
		lens->star_fstars = 0.0;

		for(j=0,m=0;j<Nregions;++j){
			lens->star_region[j] = 0.0;
			lens->star_kappa[j] = 0.0;
			lens->star_xdisk[j][0] = centers[j].x[0];
			lens->star_xdisk[j][1] = centers[j].x[1];

		}
		return;
	}

	for(j=0;j<Nregions;++j){
		m = 0;
		assert( centers[j].kappa > 0.0);

		NstarsPerImage = lens->stars_N/lens->star_Nregions;

		lens->star_region[j] = 1.0/sqrt(pi*lens->star_fstars*centers[j].kappa*lens->Sigma_crit
				/lens->star_massscale/(float)NstarsPerImage);

		// cutoff based on comparison of star deflection to smooth component
		//rcut = 4*sqrt(star_massscale/pi/Sigma_crit
		//		/( centers[j].kappa+sqrt(pow(centers[j].gamma[0],2)+pow(centers[j].gamma[1],2)) ) );

		lens->star_kappa[j] = lens->star_fstars*centers[j].kappa;
		lens->star_xdisk[j][0] = centers[j].x[0];
		lens->star_xdisk[j][1] = centers[j].x[1];

		//printf("kappa = %e  star_region = %e\n",star_kappa[j],star_region[j]);

		for(i=0;i<NstarsPerImage;++i,++m){
			//m=j*NstarsPerImage+i;
			r = lens->star_region[j]*sqrt(ran2(seed));
			theta=2*pi*ran2(seed);
			lens->stars_xp[m][0] = centers[j].x[0] + r*cos(theta);
			lens->stars_xp[m][1] = centers[j].x[1] + r*sin(theta);
			lens->stars_xp[m][2] = 0.0;
			lens->star_masses[m] = 1.0;

			// check to make see if star is in another centers region
			//   and remove it
			for(k=0;k<j;++k){
				if(  lens->star_region[k] > sqrt(pow(centers[k].x[0]-lens->stars_xp[m][0],2)
						+ pow(centers[k].x[1]-lens->stars_xp[m][1],2)) ){
					--NstarsPerImage;
					--i;
					--m;
				}
			}
			//printf("%e %e\n",stars_xp[m][0],stars_xp[m][1]);
		}
	}


	assert(m <= lens->stars_N);
	lens->stars_N = m;

	//std::printf("last star x = %e %e\n",stars_xp[stars_N-1][0],stars_xp[stars_N-1][1]);

	float dummy=0;
	//lens->star_tree = new ForceTree(lens->stars_xp,lens->stars_N,lens->star_masses,&dummy
	//		,false,false,5,2,false,lens->star_theta_force);

	lens->star_tree = new QuadTree(lens->stars_xp,lens->star_masses,&dummy,lens->stars_N
			,false,false,0,4,lens->star_theta_force);

	std::cout<<"projected with 2D tree"<<std::endl;

	// visit every branch to find center of mass and cutoff scale */
	lens->stars_implanted = true;

	return ;
}

void setStars(AnaLens *lens, bool implanted){
	lens->stars_implanted = implanted;
}

/// subtracts the mass in stars from the smooth model to compensate
/// for the mass of the stars the lensing quantities are all updated not replaced
void AnaLens::substract_stars_disks(double *ray,double *alpha
		,float *kappa,float *gamma){

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
