/*
 * implant_stars.c
 *
 *  Created on: Jun 5, 2010
 *      Author: R.B. Metcalf
 */

#include "slsimlib.h"

using namespace std;



/** 
 * \brief  Implant or randomize stars into the lens around the images.
 *
 * The first time this is called the memory for the stars is allocated.
 * On subsequent calls the star positions and masses are randomized within thier 
 * regions, but thier number and average density in stars and total
 * is not changed.  This can become out of sync with the host if its 
 * parameters are changed.
 *
 * lens->Nstars and lens->stars_fstars must be set before calling
 * allocates all memory for stars
 *
 * 
 */
void LensHalo::implant_stars(
      PosType **centers      /// Warning:: positions need to be in physical Mpc on the lens plane
      ,int Nregions          /// number of regions where stars should be implanted
      ,long *seed
      ,IMFtype type
                             ){
	PosType r,theta,NstarsPerImage;
	unsigned long i,j,m,k;
  
	if(stars_N < 1.0  || star_fstars <= 0) return;
	if(star_fstars > 1.0){ std::printf("fstars > 1.0\n"); exit(0); }
	if(!(stars_implanted) ){
    
		//star_masses = new float[stars_N];
		stars_index.resize(stars_N);
    stars_xp.resize(stars_N);
		star_theta_force = 1.0e-1;
		assert(Nregions > 0);
		star_Nregions = Nregions;
		star_region.resize(Nregions);
		star_Sigma.resize(Nregions);
		star_xdisk.resize(Nregions);
	}
  
  PosType mean_mstar[Nregions];
	if(stars_N < 1  || star_fstars <= 0){
		stars_implanted = true;
		stars_N = 0;
		star_fstars = 0.0;
    
		for(j=0,m=0;j<Nregions;++j){
			star_region[j] = 0.0;
			star_Sigma[j] = 0.0;
			mean_mstar[j] = 0.0;
			star_xdisk[j][0] = centers[j][0];
			star_xdisk[j][1] = centers[j][1];
          
		}
		return;
	}
  
	//params.get("main_sub_mass_max",sub_Mmax)
	NstarsPerImage = stars_N/star_Nregions;
  stars_xp.resize(stars_N);
  {
    std::vector<float> star_masses=stellar_mass_function(type, stars_N, seed, main_stars_min_mass, main_stars_max_mass, bend_mstar,lo_mass_slope,hi_mass_slope);
  
    for(int i=0; i<stars_N; ++i){
      stars_xp[i].Mass = star_masses[i];
    }
  }
  
	if (type==One){
		for(j=0;j<Nregions;++j){
			mean_mstar[j]=1.0;
		}
	}
	else if(type==Mono){
		for(j=0;j<Nregions;++j){
			mean_mstar[j]=main_stars_min_mass;
		}
	}
	else{
		for(j=0,m=0;j<Nregions;++j){
			mean_mstar[j] = 0.0;
			for(i=0;i<NstarsPerImage;++i,++m){
				mean_mstar[j]+=stars_xp[m].Mass;
			}
			mean_mstar[j]/=NstarsPerImage;
			//cout << "mean_mstar: " << j << " " << mean_mstar[j] << endl;
		}
	}
  
	for(j=0,m=0;j<Nregions;++j){
		PosType alpha[2];
		KappaType Sigma = 0, gamma[3];
		alpha[0]=alpha[1]=gamma[0]=gamma[1]=gamma[2]=0.;
		PosType xcm[2];
		xcm[0] = centers[j][0];
		xcm[1] = centers[j][1];
    
    if(!stars_implanted){
      KappaType phi;
      force_halo(alpha,&Sigma,gamma,&phi,xcm,false);
      star_Sigma[j] = star_fstars*Sigma;
      star_region[j] = 1.0/sqrt(PI*star_Sigma[j]/(mean_mstar[j]*(float)NstarsPerImage));
      
      star_xdisk[j][0] = centers[j][0];
      star_xdisk[j][1] = centers[j][1];
      
      //std::cout << "star disk centers " << star_xdisk[j][0] << "  " << star_xdisk[j][1] << std::endl;
      
    }
    
		//printf("kappa = %e  star_region = %e\n",star_Sigma[j],star_region[j]);
		/*char *fname = "stars0.dat";
		if(j==0){fname = "stars0.dat";}
		if(j==1){fname = "stars1.dat";}
		if(j==2){fname = "stars2.dat";}
		if(j==3){fname = "stars3.dat";}
		ofstream fstars(fname);
		if(j==0){fname = "masses0.dat";}
		if(j==1){fname = "masses1.dat";}
		if(j==2){fname = "masses2.dat";}
		if(j==3){fname = "masses3.dat";}
		ofstream mstars(fname);*/
    
		for(i=0;i<NstarsPerImage;++i,++m){
			//m=j*NstarsPerImage+i;
			r = star_region[j]*sqrt(ran2(seed));
			theta=2*PI*ran2(seed);
			stars_xp[m][0] = star_xdisk[j][0] + r*cos(theta);
			stars_xp[m][1] = star_xdisk[j][1] + r*sin(theta);
			stars_xp[m][2] = 0.0;
			//cout << m << " " << star_masses[m] << endl;
      
			//if(maxr<r){
			//	maxr=r;
			//}
			// check to make see if star is in another centers region
			//   and remove it
			//cout << "bla: " << j << " " <<  m << " " << i << endl;
			for(k=0;k<j;++k){
				if(  star_region[k] > sqrt(pow(centers[k][0]-stars_xp[m][0],2)
                                   + pow(centers[k][1]-stars_xp[m][1],2)) ){
					//--NstarsPerImage;
					//--i;
					--m;
					break;
				}
			}
			//fstars << scientific << stars_xp[m][0] << " " << stars_xp[m][1] << endl;
      //std::cout << stars_xp[m][0] << " " << stars_xp[m][1] << endl;
			//mstars << scientific << star_masses[m] << endl;
      
			//cout << "max r" << maxr << " " << star_region[j] << endl;
			//printf("%e %e\n",stars_xp[m][0],stars_xp[m][1]);
		}
		//fstars.close();
		//mstars.close();
    
	}
  
	assert(m <= stars_N);
	stars_N = m;
  
	//std::printf("last star x = %e %e\n",stars_xp[stars_N-1][0],stars_xp[stars_N-1][1]);
	//star_tree = new TreeForce(stars_xp,stars_N,star_masses,&dummy
	//		,false,false,5,2,false,star_theta_force);
  
  if(stars_implanted) delete star_tree;
	star_tree = new TreeQuadParticles<StarType>(stars_xp.data(),stars_N
                           ,false,false,0,4,star_theta_force);
  
	// visit every branch to find center of mass and cutoff scale */
	stars_implanted = true;
  
	return ;
}

/// For implanting one patch of stars. See other implant_stars.
void LensHalo::implant_stars(
                             PosType *center      /// Warning:: positions need to be in physical Mpc on the lens plane
                             ,long *seed
                             ,IMFtype type
                             ){
  PosType *tmp_p;// = new double *;
  tmp_p = center;
  
  implant_stars(&tmp_p,1,seed,type);
}
/// Un-implant stars.  Remove stars an any information about the number and size of star regions.
void LensHalo::remove_stars(){  
  
	if(stars_implanted){
    
    delete star_tree;
    stars_index.clear();;
		star_region.clear();
		star_Sigma.clear();
    star_xdisk.clear();
    star_Nregions = 0;
    
    stars_implanted = false;
	}
  
	return ;
}

/** \brief subtracts the mass in stars from the smooth model to compensate
* for the mass of the stars the lensing quantities are all updated not replaced
 */
void LensHalo::substract_stars_disks(PosType const *ray,PosType *alpha
		,KappaType *kappa,KappaType *gamma){

	if(!(stars_implanted)) return;

	PosType xcm,ycm,r;
	float tmp;
	PosType tmp_mass;
	int i;

  //std::cout <<  std::endl;
  //std::cout << "ray = " << ray[0] << " " << ray[1] << std::endl;
 	for(i=0;i<star_Nregions;++i){
		xcm = ray[0] - star_xdisk[i][0];
		ycm = ray[1] - star_xdisk[i][1];
		r=sqrt(xcm*xcm + ycm*ycm);

    //std::cout << "r/star_region[" << i << "] = " << r/star_region[i] << std::endl;

		if(r < star_region[i]){
			alpha[0] += star_Sigma[i]*xcm;
			alpha[1] += star_Sigma[i]*ycm;
			*kappa -= star_Sigma[i];
		}else{
      
			tmp_mass = star_Sigma[i]*pow(star_region[i],2)/r/r;
			alpha[0] += tmp_mass*xcm;
			alpha[1] += tmp_mass*ycm;

			tmp = 2*tmp_mass/r/r;
			gamma[0] += 0.5*(xcm*xcm-ycm*ycm)*tmp;
			gamma[1] += xcm*ycm*tmp;
		}
	}

	return;
}

/** \brief random stellar masses according to IMF of choice
 *
 * mtype defines the stellar mass function
 * 0 - always the same stellar mass (e.g. 1Msol)
 * 1 - salpeter imf, i.e. slope = -2.35
 * 2 - broken power law, requires lower mass end slope (powerlo), high mass slope (powerhi), bending point (bendmass)
 * 3 - further IMF models may follow
 */
std::vector<float> LensHalo::stellar_mass_function(IMFtype type, unsigned long Nstars, long *seed
		,PosType minmass, PosType maxmass, PosType bendmass, PosType powerlo, PosType powerhi){

	//if(!(stars_implanted)) return;
	unsigned long i;
	PosType powerp1,powerp2,powerp3,powerlp1,bendmass1,bendmass2,shiftmax,shiftmin,n0,n1,n2,n3,rndnr,tmp0,tmp1,tmp2;

  std::vector<float> stellar_masses(Nstars);
  
	if(type==One){
		for(i = 0; i < Nstars; i++){
				stellar_masses[i]=1.0;
		}
	}

	if(type==Mono){
		if((minmass!=maxmass)){
		    cerr << "For IMF type Mono main_stars_min_mass and main_stars_max_mass must be defined in parameter file and they must be equal" << endl;
		    exit(1);
		}
		for(i = 0; i < Nstars; i++){
			stellar_masses[i]=minmass;
		}
	}

    if(type==Salpeter){
    	if(minmass==maxmass){
    			    cerr << "For IMF type Salpeter main_stars_min_mass and main_stars_max_mass must be defined in parameter file" << endl;
    			    exit(1);
    	}
    	powerp1 = -1.35;
    	n0 = (pow(maxmass,powerp1)) /powerp1;
    	n1 =  n0 - (pow(minmass,powerp1)) / powerp1;
    	for(i = 0; i < Nstars; i++){
    		stellar_masses[i] = pow( (-powerp1*(n1*ran2(seed)-n0)),(1.0/powerp1) );
    	}
	}

    if(type==SinglePowerLaw){
       	if((minmass==maxmass) || (powerlo!=powerhi) || ((powerlo==0)&(powerhi==0))){
       			    cerr << "For IMF type SinglePowerLaw main_stars_min_mass, main_stars_max_mass and main_stars_lo_mass_slope must be defined in parameter file. main_stars_lo_mass_slope must be equal to main_stars_hi_mass_slope, main_stars_min_mass must be different from main_stars_max_mass!" << endl;
       			    exit(1);
       	}

       	powerp1 = powerlo+1.0;
       	n0 = (pow(maxmass,powerp1)) /powerp1;
       	n1 =  n0 - (pow(minmass,powerp1)) / powerp1;
       	for(i = 0; i < Nstars; i++){
       		stellar_masses[i] = pow( (-powerp1*(n1*ran2(seed)-n0)),(1.0/powerp1) );
       	}
   	}


    if(type==Chabrier){
		if(minmass==maxmass){
			cerr << "For IMF type Chabrier main_stars_min_mass and main_stars_max_mass must be defined in parameter file!" << endl;
			exit(1);
		}
		PosType chab_param[]={0.086,0.22,0.57};
		powerp1=-1.3;
		tmp0=2.0*chab_param[2]*chab_param[2];
		tmp1=chab_param[0]/(log(10.0)*log(10.0))*exp((-(log10(chab_param[1]))*(log10(chab_param[1])))/tmp0);
		tmp2=-0.5*chab_param[0]/log(10)*sqrt(PI)*sqrt(tmp0);
		n1=tmp2*(erff(log10(chab_param[1])/sqrt(tmp0))-erff((log10(chab_param[1])-log10(minmass))/sqrt(tmp0)));
		n2=tmp1/powerp1*(pow(maxmass,powerp1)-1.0);
		n0=n1+n2;
		for(i = 0; i < Nstars; i++){
			rndnr=ran2(seed);
			if(rndnr<(n1/n0)){
				stellar_masses[i]=pow(10.0,(log10(chab_param[1])-sqrt(2.*chab_param[2]*chab_param[2])*erfinv((n0*rndnr)/tmp2+erff((log10(chab_param[1])-log10(minmass))/(sqrt(2.*chab_param[2]*chab_param[2])) ) )));
			}
			else{
				stellar_masses[i]=pow( ((n0*rndnr-n1)*powerp1/tmp1 +1.0),(1./powerp1) );
			}
		}
	}

    if(type==BrokenPowerLaw){
    	if(powerlo==powerhi){
    		cerr << "For IMF type BrokenPowerLaw inner slope (main_stars_lo_mass_slope) and outer slope (main_stars_hi_mass_slope) must be defined in parameter file" << endl;
    		exit(1);
    	}
    	if((powerlo==-1) || (powerhi==-1) ){
    	    cerr << "For IMF type BrokenPowerLaw no slope of -1 is allowed." << endl;
    	    exit(1);
    	}

    	else{
    		powerlp1=powerlo+1.0;
    		powerp1 = powerhi+1.0;
    		shiftmax=maxmass/bendmass;
    		shiftmin=minmass/bendmass;
    		n1=(1./powerlp1)-(pow(shiftmin, powerlp1))/powerlp1;
    		n2=( pow(shiftmax,powerp1) )/powerp1-(1./powerp1);
    		n0=n1+n2;
    		for(i = 0; i < Nstars; i++){
    			rndnr=ran2(seed);
    			if(rndnr<(n1/n0)){
    				stellar_masses[i]=(pow( ((n0*rndnr)*(powerlp1)+ pow(shiftmin,powerlp1)),(1.0/powerlp1)))*bendmass;
    				//cout << stellar_masses[i] << endl;
    			}
    			else{
    				stellar_masses[i]=(pow( ((n0*rndnr-n1)*powerp1+1.0),(1.0/powerp1)))*bendmass;
    				//cout << stellar_masses[i] << endl;
    			}
    		}
    	}
    }

    if(type==Kroupa){
    	bendmass1=0.08;
    	bendmass2=0.5;
		if((minmass>bendmass1)||(maxmass<bendmass2)){
			cerr << "For IMF type Kroupa main_stars_min_mass<0.08 and main_stars_max_mass>0.5 must be defined in parameter file!" << endl;
			exit(1);
		}
		else{
			PosType kroupa_param[]={-0.3,-1.3,-2.3};
			powerp1=kroupa_param[0]+1.0;
			powerp2=kroupa_param[1]+1.0;
			powerp3=kroupa_param[2]+1.0;
			tmp1=pow(bendmass1,(kroupa_param[0]-kroupa_param[1]));
			tmp2=tmp1*pow(bendmass2,(kroupa_param[1]-kroupa_param[2]));
			n1=(pow(bendmass1,powerp1)-pow(minmass, powerp1))/powerp1;
			n2=tmp1*((pow(bendmass2,powerp2)-pow(bendmass1, powerp2))/powerp2);
			n3=tmp2*((pow(maxmass,powerp3)-pow(bendmass2, powerp3))/powerp3);
			n0=n1+n2+n3;

			for(i = 0; i < Nstars; i++){
				rndnr=ran2(seed);
				if(rndnr<(n1/n0)){
					stellar_masses[i]=(pow( ((n0*rndnr)*(powerp1)+ pow(minmass,powerp1)),(1.0/powerp1)));

				}
				else if(rndnr<((n1+n2)/n0)) {
					stellar_masses[i]=(pow( ((n0*rndnr-n1)*(powerp2/tmp1)+pow(bendmass1,powerp2)),(1.0/powerp2)));

				}
				else{
					tmp0= ((n0*rndnr-n1-n2)*(powerp3/tmp2)+pow(bendmass2,powerp3));
					stellar_masses[i]=(pow(tmp0,(1.0/powerp3)));
				}
			}
		}
    }
    return stellar_masses;
}

