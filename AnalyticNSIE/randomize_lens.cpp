/*
 * randomize_lens.c
 *
 *  Created on: Mar 23, 2010
 *      Author: R.B. Metcalf
 */

#include "slsimlib.h"

const float sheartol = 1.0e-3;

using namespace std;

/*
 * Einstein_routines for making random, close to elliptical
 *     lenses
 */


void AnaNSIELensHalo::RandomizeSigma(long *seed,bool tables){
	double *sigmaTable;
	int n,NsigmaTable;
	ifstream file;
	char *filename;

	filename = "GalaxyData/slacs_sigma.dat";
	file.open(filename);

	if(!file){
		std::cout << "Can't open file " << filename << std::endl;
		exit(1);
	}

	file >> NsigmaTable;
	sigmaTable=new double[NsigmaTable];

	for(n=0;n<NsigmaTable;++n) file >> sigmaTable[n];

	file.close();

	sigma = Utilities::RandomFromTable(sigmaTable,NsigmaTable,seed);

	delete[] sigmaTable;
}

void MultiLens::RandomizeSigma(long *seed,bool tables){
	if(flag_input_lens)
		analens->RandomizeSigma(seed,tables);
}

void AnaNSIELensHalo::RandomizeHost(long *seed,bool tables){
	double fo,*axisTable;
	int n;
	int NaxisTable;
	ifstream file;
	string filename;

	filename = "GalaxyData/slacs_f.dat";
	file.open(filename.c_str());

	if(!file){
		std::cout << "Can't open file " << filename << std::endl;
		exit(1);
	}

	file >> NaxisTable;

	axisTable=new double[NaxisTable];

	for(n=0;n<NaxisTable;++n) file >> axisTable[n];
	file.close();

	fo=fratio;

	// random host ellipticity
	do{
		if(tables) fratio = Utilities::RandomFromTable(axisTable,NaxisTable,seed);
		else fratio = fo + 0.1*gasdev(seed);
		//std::cout << "f=%e\n",axis_ratio);
	}while(fratio < 0.0 || fratio > 1.0);

	// unaligned shear modes
	if(perturb_Nmodes > 0) RandomlyDistortLens(seed,3);
	// aligned hexopole and octopole
	if(perturb_Nmodes > 0){
		for(n=3;n<6;++n) AlignedRandomlyDistortLens(seed
			,pa+3*pi*gasdev(seed)/180,n);
	}

	delete[] axisTable;

	return ;
}

void MultiLens::RandomizeHost(long *seed,bool tables){
	if(flag_input_lens)
		analens->RandomizeHost(seed,tables);
}

/** \ingEinstein_roup ChangeLens
 *
 *
 *	*  make a random realization of the perturbations to the lens
	 *    the normalization is such that the rms maximum surface
	 *    density of the k-th mode is rms[i] times the surface density
	 *    of the host at the Einstein radius if it were a SIS, if beta=1
	 *    this is true at all radii.
	 *    each mode is randomly oriented
	 */

void AnaNSIELensHalo::RandomlyDistortLens(long *seed, int Nmodes){

	int i,k;
	double tmp,theta;

	if(Nmodes > 0){
		// lognormal shear and kappa distribution
		perturb_modes[0]=0.015*pow(10,gasdev(seed)*perturb_rms[0]);
		tmp=0.015*pow(10,gasdev(seed)*perturb_rms[1]);
		theta=2*pi*ran2(seed);
		perturb_modes[1] = tmp*cos(theta);
		perturb_modes[2] = tmp*sin(theta);

		/*
		modes[0] = rms_perturb[0]*gasdev(seed);
		if(Nmodes > 2){
			modes[1] = rms_perturb[1]*gasdev(seed)/sqrt(2.);
			modes[2] = rms_perturb[1]*gasdev(seed)/sqrt(2.);
		}
		modes[0]= sqrt(pow(modes[1],2) + pow(modes[2],2));
		 */

		if(Nmodes > 2) perturb_modes[3] = perturb_rms[2]*pow(Einstein_ro,2-perturb_beta)
					*gasdev(seed)/ perturb_beta/ perturb_beta;
		for(i=4;i<((perturb_Nmodes < Nmodes) ? perturb_Nmodes : Nmodes) ;i+=2){
			k=i/2;
			perturb_modes[i] =   perturb_rms[k+1]*pow(Einstein_ro,2- perturb_beta)
				*gasdev(seed)/( perturb_beta* perturb_beta-k*k)/sqrt(2);
			perturb_modes[i+1] = perturb_rms[k+1]*pow(Einstein_ro,2- perturb_beta)
				*gasdev(seed)/( perturb_beta* perturb_beta-k*k)/sqrt(2);

			//std::cout << "k=%i i=%i  modes = %e %e rms_purturb=%e\n",k,i,modes[i],modes[i+1],
			//		rms_perturb[k+1]);
		}
	}
	//PrintAnaLens(lens,false,false);

	return ;
}
/** \ingEinstein_roup ChangeLens
 *
 *  make a random realization of the perturbations to the lens
 *    the normalization is such that the rms maximum surface
 *    density of the k-th mode is rms[i] times the surface density
 *    of the host at the Einstein radius if it were a SIS, if beta=1
 *    this is true at all radii.
 *    This is for randomly distorting a mode but maintaining its alignment
 *    with the ellipticity for example.
 */
void AnaNSIELensHalo::AlignedRandomlyDistortLens(long *seed,double theta,int Npole){

	double tmp;
	int i,k;

	for(i=4;i<perturb_Nmodes;i+=2){
		k=i/2;

		if( k+1 == Npole){
			tmp=gasdev(seed);
			perturb_modes[i] =   perturb_rms[k+1]*pow(Einstein_ro,2- perturb_beta)
						*tmp*cos(theta)/( perturb_beta* perturb_beta-k*k);
			perturb_modes[i+1] = perturb_rms[k+1]*pow(Einstein_ro,2- perturb_beta)
						*tmp*sin(theta)/( perturb_beta* perturb_beta-k*k);
		}
		//std::cout << "k=%i i=%i  modes = %e %e rms_purturb=%e\n",k,i,modes[i],modes[i+1],
		//		rms_perturb[k+1]);
	}

	return ;
}
/** \ingEinstein_roup ChangeLens
 *
 */
void AnaNSIELensHalo::RandomizeSubstructure2(double rangeInRei,long *seed){
	long i,k;
	double r,theta,rmax,rav[2],r2av=0,area_av=0;
	static unsigned long NsubMax;
	static double Rm = 0.0,scale = 0.0,ndensity = 0.0,Einstein_ro_save = 0.0;

	rav[0]=rav[1]=0.0;


	// scaling of substructure masses with host sigma
	//  exponent of 3 is for  M_host pEinstein_ropto sigma^3 model
	if(sigma > 0.0){
		scale = pow(sigma/sub_sigmaScale,3);

		// keep fraction of surface density at r = f R(sigma) constant
		//   R(sigma) pEinstein_ropto sigma to make all hosts have the same average density
		ndensity = sub_Ndensity*pow(sigma/sub_sigmaScale,1)/scale;
	}

	if(Einstein_ro > 0.0){
		Rm = Einstein_ro*rangeInRei + sub_Rmax*pow(scale,1./3.)
	          + pow(2*sub_Mmax*scale*Einstein_ro/pi/Sigma_crit/sheartol,1./3.);
		Einstein_ro_save = Einstein_ro;
	}

	if(!(substruct_implanted) && ndensity > 0){
		NsubMax=(unsigned long)(ndensity*pi*Rm*Rm + 5*sqrt(ndensity*pi*Rm*Rm) );
		if(NsubMax > 0){
			sub_x=Utilities::PosTypeMatrix(NsubMax,2);
			switch(sub_type){
			case pointmass:
				subs=new LensHalo[NsubMax];
				break;
			case powerlaw:
				subs=new PowerLawLensHalo[NsubMax];
				break;
			case nfw:
				subs=new NFWLensHalo[NsubMax];
				break;
			default:
				subs=new LensHalo[NsubMax];
				break;
			}
			sub_substructures = new IndexType[NsubMax];
		}
		substruct_implanted=true;
	}
	//std::cout << "Rm/re = %e\n",Rm/Einstein_ro);
	//for(i=0;i<12;++i) std::cout << "%f %f\n",poidev(ndensity*pi*Rm*Rm,seed),ndensity*pi*Rm*Rm);

	unsigned int Nsub=(int)(poidev(ndensity*pi*Rm*Rm,seed));
	Nsub = (NsubMax > Nsub) ? Nsub : NsubMax ;

	//std::cout << "scale = %e\n",scale);

	for(i=0,k=0; i < Nsub;++i){
		//for(i=0;i<NSubstruct;++i){

		r=Rm*sqrt(ran2(seed));

		do{
			if(sub_alpha == -1.0){
				float ratio = sub_Mmax/sub_Mmin;
				subs[k].set_mass(sub_Mmin*scale*pow(ratio,ran2(seed)));
			}else{
				subs[k].set_mass(sub_Mmin*scale
						*pow(ran2(seed)*(pow(sub_Mmax/sub_Mmin,sub_alpha+1)-1)+1.0
								,1.0/(sub_alpha+1)));
			}
		}while(subs[k].get_mass() < sub_Mmin);  // not sure why this is necessary

		// average density of a substructure does not scale with host
		subs[k].set_Rmax(sub_Rmax*pow(scale,1./3.)
				*pow(subs[k].get_mass()/sub_Mmax/scale,1/3.));

		subs[k].set_slope(sub_beta);

		subs[k].set_rscale(0.1*subs[k].get_Rmax());

		// maximum radius for a substructure of this mass
		rmax = (Einstein_ro_save*rangeInRei + subs[k].get_Rmax()
		     + pow(2*subs[k].get_mass()*Einstein_ro_save/pi/Sigma_crit/sheartol,1./3.) );

		//std::cout << "RcutSubstruct[%i] = %e\n",k,RcutSubstruct[k]);
		//std::cout << "%e %e %e Rm=%e\n",r/rmax,r,rmax,Rm);
		if( r < rmax){
			theta=2*pi*ran2(seed);
			sub_x[k][0]=r*cos(theta);
			sub_x[k][1]=r*sin(theta);
			assert(k<NsubMax);

			rav[0] += sub_x[k][0];
			rav[1] += sub_x[k][1];
			r2av += r*r;
			area_av += pow(subs[k].get_Rmax(),2);
			++k;
		}
	}

	sub_N = k;

/*
	if(sub_N > 2){
		// decide if it is worth doing substructure force calculation by tree or
		//   by direct summation
		r2av /= sub_N;
		r2av = r2av - (rav[0]*rav[0]+rav[1]*rav[1])/(sub_N)/(sub_N);
		area_av /= sub_N;

		//std::cout << "relative area of average particle %e \n",area_av/r2av);

		assert(r2av > 0);
		assert(area_av >= 0.0);

		if(sub_type == pointmass){  // always use tree for point mass substructures
			sub_theta_force=1.0e-1;
			sub_tree=BuildTreeNB(sub_x,sub_Rcut,sub_mass,
				false,true,sub_N,sub_substructures,2,sub_theta_force);
		}else if( area_av/r2av < 1.0e-2 && sub_N > 300 ){
			// build tree for doing substructure force calculation
			sub_theta_force=1.0e-1;
			sub_tree=BuildTreeNB(sub_x,sub_Rcut,sub_mass,
				true,true,sub_N,sub_substructures,2,sub_theta_force);
		}else sub_tree = 0;
	}else sub_tree = 0;
*/
	return;
}
/** \ingEinstein_roup ChangeLens
	 *
	 *  get a random population of substructures
	 *  This version does not scale to Einstein_ro and
	 *  sigma
	 */
void AnaNSIELensHalo::RandomizeSubstructure3(double rangeInRei,long *seed){
	long i,k;
	double r,theta,rmax,rav[2],r2av=0,area_av=0;
	static unsigned long NsubMax;
	static double Rm = 0.0,Einstein_ro_save = 0.0;

	rav[0]=rav[1]=0.0;

	if(Einstein_ro > 0.0){
		Einstein_ro_save = Einstein_ro;
		//std::cout << "Einstein_ro = %e\n",Einstein_ro);
	}else if( perturb_Nmodes >= 4 && perturb_modes[3] > 0 ){  // this is in case the the host has already been fit to image positions
		Einstein_ro_save = perturb_modes[3]/2;
	}

	assert(Einstein_ro_save > 0.0);

	Rm = Einstein_ro_save*rangeInRei;
	Rm +=  sub_Rmax
          + pow(2*sub_Mmax*Einstein_ro_save/pi/Sigma_crit/sheartol,1./3.);

	assert(Rm > 0.0);

	if(!substruct_implanted){
		NsubMax=(unsigned long)(sub_Ndensity*pi*Rm*Rm*(1+5/sqrt(sub_Ndensity*pi*Rm*Rm)) );
		sub_x=Utilities::PosTypeMatrix(NsubMax,2);
		switch(sub_type){
		case pointmass:
			subs=new LensHalo[NsubMax];
			break;
		case powerlaw:
			subs=new PowerLawLensHalo[NsubMax];
			break;
		case nfw:
			subs=new NFWLensHalo[NsubMax];
			break;
		default:
			subs=new LensHalo[NsubMax];
			break;
		}
		substruct_implanted=true;
		sub_substructures = new IndexType[NsubMax];
	}
	//std::cout << "Rm/re = %e\n",Rm/Einstein_ro);
	//for(i=0;i<12;++i) std::cout << "%f %f\n",poidev(ndensity*pi*Rm*Rm,seed),ndensity*pi*Rm*Rm);

	unsigned int Nsub=(int)(poidev(sub_Ndensity*pi*Rm*Rm,seed));

	assert(Nsub < NsubMax);

	Nsub = (NsubMax > Nsub) ? Nsub : NsubMax ;

	//std::cout << "scale = %e\n",scale);
	for(i=0,k=0; i < Nsub;++i){
		//for(i=0;i<NSubstruct;++i){

		r=Rm*sqrt(ran2(seed));

		do{
			if(sub_alpha == -1.0){
				float ratio = sub_Mmax/sub_Mmin;
				subs[k].set_mass(sub_Mmin*pow(ratio,ran2(seed)));
			}else{
				subs[k].set_mass(sub_Mmin
						*pow(ran2(seed)*(pow(sub_Mmax/sub_Mmin,sub_alpha+1)-1)+1.0
								,1.0/(sub_alpha+1)));
			}
		}while(subs[k].get_mass() < sub_Mmin);  // not sure why this is necessary

		// average density of a substructure does not scale with host
		subs[k].set_Rmax(sub_Rmax
				*pow(subs[k].get_mass()/sub_Mmax,1/3.));

		// maximum radius for a substructure of this mass
		rmax = (Einstein_ro_save*rangeInRei + subs[k].get_Rmax()
		     + pow(2*subs[k].get_mass()*Einstein_ro_save/pi/Sigma_crit/sheartol,1./3.) );

		//std::cout << "RcutSubstruct[%i] = %e\n",k,RcutSubstruct[k]);
		//std::cout << "%e %e %e Rm=%e\n",r/rmax,r,rmax,Rm);
		if( r < rmax){
			theta=2*pi*ran2(seed);
			sub_x[k][0]=r*cos(theta);
			sub_x[k][1]=r*sin(theta);
			assert(k<NsubMax);

			rav[0] += sub_x[k][0];
			rav[1] += sub_x[k][1];
			r2av += r*r;
			area_av += pow(subs[k].get_Rmax(),2);
			++k;
		}
	}

	sub_N = k;

	//std::cout << "sub_N = %li Nsub = %li ndensity =%e\n",sub_N,Nsub,sub_Ndensity);
/*
	if(sub_N > 2){
		// decide if it is worth doing substructure force calculation by tree or
		//   by direct summation
		r2av /= sub_N;
		r2av = r2av - (rav[0]*rav[0]+rav[1]*rav[1])/(sub_N)/(sub_N);
		area_av /= sub_N;

		//std::cout << "relative area of average particle %e \n",area_av/r2av);

		assert(r2av > 0);
		assert(area_av >= 0.0);

		if(sub_type == pointmass){  // always use tree for point mass substructures
			sub_theta_force=1.0e-1;
			sub_tree=BuildTreeNB(sub_x,sub_Rcut,sub_mass,
				false,true,sub_N,sub_substructures,2,sub_theta_force);
		}else if( area_av/r2av < 1.0e-2 && sub_N > 300 ){
			// build tree for doing substructure force calculation
			sub_theta_force=1.0e-1;
			sub_tree=BuildTreeNB(sub_x,sub_Rcut,sub_mass,
				true,true,sub_N,sub_substructures,2,sub_theta_force);
		}else sub_tree = 0;
	}else sub_tree = 0;
*/
	return;
}

namespace Utilities{
/** \ingEinstein_roup Utill
 * \brief Generates a random deviates drawn fEinstein_rom appEinstein_roximately the same as the values of table
 *
 */
	double RandomFromTable(double *table,unsigned long Ntable,long *seed){
		double y;
		unsigned long j;

		y=ran2(seed)*(Ntable-1);
		j=(int)(y);
		return (table[j+1]-table[j])*(y-j) + table[j];
	}
}

double AnaNSIELensHalo::FractionWithinRe(double rangeInRei){
	double B;

	B = (sub_Rmax/pow(sub_Mmax,1./3.)
			+ pow(2*Einstein_ro/pi/Sigma_crit/1.e-3,1./3.) );

	return 1+(1+sub_alpha)*(
		2*rangeInRei*Einstein_ro*B*(
				pow(sub_Mmax,sub_alpha+4./3.)-pow(sub_Mmin,sub_alpha+4./3.)
				)/(sub_alpha+4./3.)
		+ B*B*(
				pow(sub_Mmax,sub_alpha+5./3.)-pow(sub_Mmin,sub_alpha+5./3.)
				)/(sub_alpha+5./3.)
		)/(pow(rangeInRei*Einstein_ro,2)*(pow(sub_Mmax,sub_alpha+1.0)-pow(sub_Mmin,sub_alpha+1.0))
				);
}

double BaseNSIELensHalo::averageSubMass(){
	// average mass of substructures
	return sub_Mmax*(sub_alpha+1)
				  /(sub_alpha+2)*(1-pow(sub_Mmin/sub_Mmax,sub_alpha+2))/
				  (1-pow(sub_Mmin/sub_Mmax,sub_alpha+1));
}
