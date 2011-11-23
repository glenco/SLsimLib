/*
 * implant_stars.c
 *
 *  Created on: Jun 5, 2010
 *      Author: R.B. Metcalf
 */

#include <slsimlib.h>

/** \ingroup ChangeLens
 * \brief  Implants stars into the lens around the images.
 *
 * lens->Nstars and lens->stars_fstars must be set before calling
 * allocates all memory for stars
 */

void implant_stars(AnaLens *lens,Point *centers,unsigned long Nregions,long *seed){
	PosType r,theta,NstarsPerImage;
	static PosType **coord;
	unsigned long i,j,m,k;

	if(lens->stars_N < 1.0  || lens->star_fstars <= 0) return;
	if(lens->star_fstars > 1.0){ printf("fstars > 1.0\n"); exit(0); }
	if(!(lens->stars_implanted) ){

		coord = PosTypeMatrix(0,2,0,2);
		coord[0][0] = coord[1][1] = coord[2][2] = 1.0;
		coord[0][1] = coord[0][2] = coord[1][0] = coord[1][2] = 0.0;
		coord[2][0] = coord[2][1] = 0.0;

		lens->star_masses = (float *) calloc(lens->stars_N,sizeof(float));
		lens->stars = (unsigned long *) calloc(lens->stars_N,sizeof(unsigned long));
		lens->stars_xp = PosTypeMatrix(0,lens->stars_N-1,0,2);
		lens->star_theta_force = 1.0e-1;

		assert(Nregions > 0);
		lens->star_Nregions = Nregions;
		lens->star_region = (double *)calloc(Nregions,sizeof(double));
		lens->star_kappa = (double *)calloc(Nregions,sizeof(double));
		lens->star_xdisk = dmatrix(0,Nregions-1,0,1);

	}else{
		// free star_tree
		delete lens->star_tree;
	}

//#pragma omp parallel for private(j, NstarsPerImage, r, theta, m) firstprivate(i, k)
	for(j=0;j<Nregions;++j){
		m = 0;
		assert( centers[j].kappa > 0.0);

		NstarsPerImage = lens->stars_N/lens->star_Nregions;

		lens->star_region[j] = 1.0/sqrt(pi*lens->star_fstars*centers[j].kappa*lens->Sigma_crit
				/lens->star_massscale/(float)NstarsPerImage);

		// cutoff based on comparison of star deflection to smooth component
		//rcut = 4*sqrt(lens->star_massscale/pi/lens->Sigma_crit
		//		/( centers[j].kappa+sqrt(pow(centers[j].gamma[0],2)+pow(centers[j].gamma[1],2)) ) );

		lens->star_kappa[j] = lens->star_fstars*centers[j].kappa;
		lens->star_xdisk[j][0] = centers[j].x[0];
		lens->star_xdisk[j][1] = centers[j].x[1];

		//printf("kappa = %e  star_region = %e\n",lens->star_kappa[j],lens->star_region[j]);

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
			//printf("%e %e\n",lens->stars_xp[m][0],lens->stars_xp[m][1]);
		}
	}

	assert(m <= lens->stars_N);
	lens->stars_N = m;

	//printf("last star x = %e %e\n",lens->stars_xp[lens->stars_N-1][0],lens->stars_xp[lens->stars_N-1][1]);

	float dummy=0;
	lens->star_tree = new ForceTree(lens->stars_xp,lens->stars_N,lens->star_masses,&dummy
			,false,true,5,2,false,lens->star_theta_force);

	//printf("projected with 2D tree\n");

	// visit every branch to find center of mass and cutoff scale */
	lens->stars_implanted = true;

	return ;
}

/// subtracts the mass in stars from the smooth model to compensate
/// for the mass of the stars the lensing quantities are all updated not replaced
void substract_stars_disks(AnaLens *lens,PosType *ray,PosType *alpha
		,PosType *kappa,PosType *gamma){

	if(!(lens->stars_implanted)) return;

	PosType xcm,ycm,r,tmp;
	double mass;
	int i;

	for(i=0;i<lens->star_Nregions;++i){
		xcm = lens->star_xdisk[i][0] - ray[0];
		ycm = lens->star_xdisk[i][1] - ray[1];
		r=sqrt(xcm*xcm + ycm*ycm);

		if(r < lens->star_region[i]){
			alpha[0] += lens->star_kappa[i]*xcm;
			alpha[1] += lens->star_kappa[i]*ycm;
			*kappa -= lens->star_kappa[i];
		}else{

			mass = lens->star_kappa[i]*pow(lens->star_region[i],2)/r/r;
			alpha[0] += mass*xcm;
			alpha[1] += mass*ycm;

			tmp = 2*mass/r/r;
			gamma[0] += 0.5*(xcm*xcm-ycm*ycm)*tmp;
			gamma[1] += xcm*ycm*tmp;
		}
	}

	return;
}
