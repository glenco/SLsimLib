/*
 * randomize_lens.c
 *
 *  Created on: Mar 23, 2010
 *      Author: R.B. Metcalf
 */

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <nr.h>
#include <nrD.h>
#include <nrutil.h>
#include <analytic_lens.h>
#define sheartol 1.0e-3

/*
 * routines for making random, close to elliptical
 *     lenses
 */

//extern COSMOLOGY cosmo;

void RandomizeHost(AnaLens *lens,double r_source_phys,long *seed,Boolean tables
		,CosmoHndl cosmo){
	static double fo=0.0,*axisTable,*sigmaTable,**zTable;
	int n,i;
	static int init=0,NaxisTable,NreTable,NsigmaTable,NzTable;
	FILE *file;
	char filename[50];
	//double re_onsource;

	if(init==0 && tables){

		printf("> reading lens distribution tables");
		++init;
		// read in Einstein radius projected onto the source plane

		sprintf(filename,"GalaxyData/z_table.txt");
		printf(">reading from %s\n",filename);
		file=fopen(filename,"r");
		fscanf(file,"%i",&NzTable);
		//printf(">%i\n",NzTable);
		zTable=(double **) dmatrix(0,NzTable-1,0,1);
		for(n=0;n<NzTable;++n) fscanf(file,"%le  %le",&zTable[n][0],&zTable[n][1]);
		//for(n=0;n<NzTable;++n) printf("%le  %le\n",zTable[n][0],zTable[n][1]);
		fclose(file);

		sprintf(filename,"GalaxyData/slacs_sigma.dat");
		printf(">reading from %s\n",filename);
		file=fopen(filename,"r");
		fscanf(file,"%i",&NsigmaTable);
		sigmaTable=(double *)calloc(NsigmaTable,sizeof(double));
		//printf("sigma\n");
		for(n=0;n<NsigmaTable;++n){
			fscanf(file,"%le",&sigmaTable[n]);
			//printf("%e\n",sigmaTable[n]);
		}
		fclose(file);

		sprintf(filename,"GalaxyData/slacs_f.dat");
		printf(">reading from %s\n",filename);
		file=fopen(filename,"r");
		fscanf(file,"%i",&NaxisTable);
		//printf("%i\n",NaxisTable);
		axisTable=(double *)calloc(NaxisTable,sizeof(double));
		//printf("a/b\n");
		for(n=0;n<NaxisTable;++n){
			fscanf(file,"%le",&axisTable[n]);
			//printf("%le\n",axisTable[n]);
		}
		fclose(file);
	}


	//if(tables){
	//	re_onsource=RandomFromTable(reTable,NreTable,seed);
	//}else{
	//re_onsource=1.2e-2 + 5.5e-3*gasdev(seed);
	//}
	//do{ re_onsource=12.0e-3 + 5.5e-3*gasdev(seed);
	//}while(re_onsource < 1.0e-3 || re_onsource > 25.0e-3);

	//lens->r_source=r_source_phys*lens->ro/re_onsource;

	// random host sigma
	if(tables){
		// choose random set of redshifts
		i=(int)(NzTable*ran2(seed));
		lens->zsource = zTable[i][0];
		lens->zlens = zTable[i][1];

		//lens->zlens = RandomFromTable(zlTable,NzlTable,seed);
		//lens->zsource = RandomFromTable(zsTable,NzsTable,seed);

		lens->host_sigma = RandomFromTable(sigmaTable,NsigmaTable,seed);
		lens->host_ro = 4*pi*pow(lens->host_sigma/2.99792e5,2)
		        *angDist(lens->zlens,lens->zsource,cosmo)*angDist(0,lens->zlens,cosmo)
				/angDist(0,lens->zsource,cosmo);
		lens->Sigma_crit = pow(lens->host_sigma/2.99792e5,2)/Grav/lens->host_ro;
		lens->source_r = r_source_phys*angDist(0,lens->zlens,cosmo)
		                /angDist(0,lens->zsource,cosmo);
		lens->MpcToAsec=60*60*180*(1+lens->zsource)/pi/angDist(0,lens->zlens,cosmo);
	}

	if(fo==0.0) fo=lens->host_axis_ratio;

	// random host ellipticity
	do{
		if(tables) lens->host_axis_ratio = RandomFromTable(axisTable,NaxisTable,seed);
		else lens->host_axis_ratio = fo + 0.1*gasdev(seed);
		//printf("f=%e\n",lens->axis_ratio);
	} while(lens->host_axis_ratio < 0.0 || lens->host_axis_ratio > 1.0);

	// unaligned shear modes
	if(lens->perturb_Nmodes > 0) RandomlyDistortLens(lens,seed,3);
	// aligned hexopole and octopole
	if(lens->perturb_Nmodes > 0){
		for(n=3;n<6;++n) AlignedRandomlyDistortLens(lens,seed
			,lens->host_pos_angle+3*pi*gasdev(seed)/180,n);
//		for(n=3;n<6;++n) AlignedRandomlyDistortLens(lens,seed
//			,lens->theta+2*pi*ran2(seed),n);
	}

	return ;
}

void RandomlyDistortLens(AnaLens *lens,long *seed, int Nmodes){
	/*  make a random realization of the perturbations to the lens
	 *    the normalization is such that the rms maximum surface
	 *    density of the k-th mode is rms[i] times the surface density
	 *    of the host at the Einstein radius if it were a SIS, if beta=1
	 *    this is true at all radii
	 *    each mode is randomly oriented
	 */

	int i,k;
	double tmp,theta;

	if(Nmodes > 0){
		// lognormal shear and kappa distribution
		lens->perturb_modes[0]=0.015*pow(10,gasdev(seed)*lens->perturb_rms[0]);
		tmp=0.015*pow(10,gasdev(seed)*lens->perturb_rms[1]);
		theta=2*pi*ran2(seed);
		lens->perturb_modes[1] = tmp*cos(theta);
		lens->perturb_modes[2] = tmp*sin(theta);

		/*
		lens->modes[0] = lens->rms_perturb[0]*gasdev(seed);
		if(Nmodes > 2){
			lens->modes[1] = lens->rms_perturb[1]*gasdev(seed)/sqrt(2.);
			lens->modes[2] = lens->rms_perturb[1]*gasdev(seed)/sqrt(2.);
		}
		lens->modes[0]= sqrt(pow(lens->modes[1],2) + pow(lens->modes[2],2));
		 */

		if(Nmodes > 2) lens->perturb_modes[3] = lens->perturb_rms[2]*pow(lens->host_ro,2-lens->perturb_beta)
					*gasdev(seed)/ lens->perturb_beta/ lens->perturb_beta;
		for(i=4;i<((lens->perturb_Nmodes < Nmodes) ? lens->perturb_Nmodes : Nmodes) ;i+=2){
			k=i/2;
			lens->perturb_modes[i] =   lens->perturb_rms[k+1]*pow(lens->host_ro,2- lens->perturb_beta)
				*gasdev(seed)/( lens->perturb_beta* lens->perturb_beta-k*k)/sqrt(2);
			lens->perturb_modes[i+1] = lens->perturb_rms[k+1]*pow(lens->host_ro,2- lens->perturb_beta)
				*gasdev(seed)/( lens->perturb_beta* lens->perturb_beta-k*k)/sqrt(2);

			//printf("k=%i i=%i  modes = %e %e rms_purturb=%e\n",k,i,lens->modes[i],lens->modes[i+1],
			//		lens->rms_perturb[k+1]);
		}
	}
	//PrintAnaLens(lens,False,False);

	return ;
}
void AlignedRandomlyDistortLens(AnaLens *lens,long *seed,double theta,int Npole){
	/*  make a random realization of the perturbations to the lens
	 *    the normalization is such that the rms maximum surface
	 *    density of the k-th mode is rms[i] times the surface density
	 *    of the host at the Einstein radius if it were a SIS, if beta=1
	 *    this is true at all radii
	 *    each mode is randomly oriented
	 */

	double tmp;
	int i,k;

	for(i=4;i<lens->perturb_Nmodes;i+=2){
		k=i/2;

		if( k+1 == Npole){
			tmp=gasdev(seed);
			lens->perturb_modes[i] =   lens->perturb_rms[k+1]*pow(lens->host_ro,2- lens->perturb_beta)
						*tmp*cos(theta)/( lens->perturb_beta* lens->perturb_beta-k*k);
			lens->perturb_modes[i+1] = lens->perturb_rms[k+1]*pow(lens->host_ro,2- lens->perturb_beta)
						*tmp*sin(theta)/( lens->perturb_beta* lens->perturb_beta-k*k);
		}
		//printf("k=%i i=%i  modes = %e %e rms_purturb=%e\n",k,i,lens->modes[i],lens->modes[i+1],
		//		lens->rms_perturb[k+1]);
	}

	//PrintAnaLens(lens,False,False);

	return ;
}

void RandomizeSubstructure2(AnaLens *lens,double rangeInRei,long *seed){
	long i,k;
	double r,theta,rmax,rav[2],r2av=0,area_av=0;
	static unsigned long NsubMax;
	static double Rmax = 0.0,scale = 0.0,ndensity = 0.0,host_ro_save = 0.0;

	rav[0]=rav[1]=0.0;


	// scaling of substructure masses with host sigma
	//  exponent of 3 is for  M_host propto sigma^3 model
	if(lens->host_sigma > 0.0){
		scale = pow(lens->host_sigma/lens->sub_sigmaScale,3);

		// keep fraction of surface density at r = f R(sigma) constant
		//   R(sigma) propto sigma to make all hosts have the same average density
		ndensity = lens->sub_Ndensity*pow(lens->host_sigma/lens->sub_sigmaScale,1)/scale;
	}

	if(lens->host_ro > 0.0){
		Rmax = lens->host_ro*rangeInRei + lens->sub_Rmax*pow(scale,1./3.)
	          + pow(2*lens->sub_Mmax*scale*lens->host_ro/pi/lens->Sigma_crit/sheartol,1./3.);
		host_ro_save = lens->host_ro;
	}

	if(!(lens->substruct_implanted) && ndensity > 0){
		NsubMax=(unsigned long)(ndensity*pi*Rmax*Rmax + 5*sqrt(ndensity*pi*Rmax*Rmax) );
		if(NsubMax > 0){
			lens->sub_x=dmatrix(0,NsubMax-1,0,1);
			lens->sub_Rcut=(float *)calloc(NsubMax,sizeof(float));
			lens->sub_mass=(float *)calloc(NsubMax,sizeof(float));
			lens->sub_substructures = (IndexType *)calloc(NsubMax,sizeof(IndexType));
		}
		lens->substruct_implanted=True;
	}
	//printf("Rmax/re = %e\n",Rmax/lens->ro);
	//for(i=0;i<12;++i) printf("%f %f\n",poidev(ndensity*pi*Rmax*Rmax,seed),ndensity*pi*Rmax*Rmax);

	unsigned int Nsub=(int)(poidev(ndensity*pi*Rmax*Rmax,seed));
	Nsub = (NsubMax > Nsub) ? Nsub : NsubMax ;

	//printf("scale = %e\n",scale);

	for(i=0,k=0; i < Nsub;++i){
		//for(i=0;i<lens->NSubstruct;++i){

		r=Rmax*sqrt(ran2(seed));

		do{
			if(lens->sub_alpha == -1.0){
				lens->sub_mass[k] = lens->sub_Mmin*scale
						*pow(lens->sub_Mmax/lens->sub_Mmin,ran2(seed));
			}else{
				lens->sub_mass[k] = lens->sub_Mmin*scale
						*pow(ran2(seed)*(pow(lens->sub_Mmax/lens->sub_Mmin,lens->sub_alpha+1)-1)+1.0
								,1.0/(lens->sub_alpha+1));
			}
		}while(lens->sub_mass[k] < lens->sub_Mmin);  // not sure why this is necessary

		// average density of a substructure does not scale with host
		if(lens->sub_type != pointmass) lens->sub_Rcut[k]=lens->sub_Rmax*pow(scale,1./3.)
				*pow(lens->sub_mass[k]/lens->sub_Mmax/scale,1/3.);

		// maximum radius for a substructure of this mass
		rmax = (host_ro_save*rangeInRei + lens->sub_Rcut[k]
		     + pow(2*lens->sub_mass[k]*host_ro_save/pi/lens->Sigma_crit/sheartol,1./3.) );

		//printf("lens->RcutSubstruct[%i] = %e\n",k,lens->RcutSubstruct[k]);
		//printf("%e %e %e Rmax=%e\n",r/rmax,r,rmax,Rmax);
		if( r < rmax){
			theta=2*pi*ran2(seed);
			lens->sub_x[k][0]=r*cos(theta);
			lens->sub_x[k][1]=r*sin(theta);
			assert(k<NsubMax);

			rav[0] += lens->sub_x[k][0];
			rav[1] += lens->sub_x[k][1];
			r2av += r*r;
			area_av += pow(lens->sub_Rcut[k],2);
			++k;
		}
	}

	lens->sub_N = k;

/*
	if(lens->sub_N > 2){
		// decide if it is worth doing substructure force calculation by tree or
		//   by direct summation
		r2av /= lens->sub_N;
		r2av = r2av - (rav[0]*rav[0]+rav[1]*rav[1])/(lens->sub_N)/(lens->sub_N);
		area_av /= lens->sub_N;

		//printf("relative area of average particle %e \n",area_av/r2av);

		assert(r2av > 0);
		assert(area_av >= 0.0);

		if(lens->sub_type == pointmass){  // always use tree for point mass substructures
			lens->sub_theta_force=1.0e-1;
			lens->sub_tree=BuildTreeNB(lens->sub_x,lens->sub_Rcut,lens->sub_mass,
				False,True,lens->sub_N,lens->sub_substructures,2,lens->sub_theta_force);
		}else if( area_av/r2av < 1.0e-2 && lens->sub_N > 300 ){
			// build tree for doing substructure force calculation
			lens->sub_theta_force=1.0e-1;
			lens->sub_tree=BuildTreeNB(lens->sub_x,lens->sub_Rcut,lens->sub_mass,
				True,True,lens->sub_N,lens->sub_substructures,2,lens->sub_theta_force);
		}else lens->sub_tree = 0;
	}else lens->sub_tree = 0;
*/
	return;
}
void RandomizeSubstructure3(AnaLens *lens,double rangeInRei,long *seed){
	/*
	 *  get a random population of substructures
	 *  This version does not scale to lens->host_ro and
	 *  lens->host_sigma
	 */
	long i,k;
	double r,theta,rmax,rav[2],r2av=0,area_av=0;
	static unsigned long NsubMax;
	static double Rmax = 0.0,host_ro_save = 0.0;

	rav[0]=rav[1]=0.0;

	if(lens->host_ro > 0.0){
		host_ro_save = lens->host_ro;
		//printf("lens->host_ro = %e\n",lens->host_ro);
	}else if( lens->perturb_Nmodes >= 4 && lens->perturb_modes[3] > 0 ){  // this is in case the the host has already been fit to image positions
		host_ro_save = lens->perturb_modes[3]/2;
	}

	assert(host_ro_save > 0.0);

	Rmax = host_ro_save*rangeInRei;
	Rmax +=  lens->sub_Rmax
          + pow(2*lens->sub_Mmax*host_ro_save/pi/lens->Sigma_crit/sheartol,1./3.);

	assert(Rmax > 0.0);

	if(!lens->substruct_implanted){
		NsubMax=(unsigned long)(lens->sub_Ndensity*pi*Rmax*Rmax*(1+5/sqrt(lens->sub_Ndensity*pi*Rmax*Rmax)) );
		lens->sub_x=dmatrix(0,NsubMax-1,0,1);
		lens->sub_Rcut=(float *)calloc(NsubMax,sizeof(float));
		lens->sub_mass=(float *)calloc(NsubMax,sizeof(float));
		lens->substruct_implanted=True;
		lens->sub_substructures = (IndexType *)calloc(NsubMax,sizeof(IndexType));
	}
	//printf("Rmax/re = %e\n",Rmax/lens->ro);
	//for(i=0;i<12;++i) printf("%f %f\n",poidev(ndensity*pi*Rmax*Rmax,seed),ndensity*pi*Rmax*Rmax);

	unsigned int Nsub=(int)(poidev(lens->sub_Ndensity*pi*Rmax*Rmax,seed));

	assert(Nsub < NsubMax);

	Nsub = (NsubMax > Nsub) ? Nsub : NsubMax ;

	//printf("scale = %e\n",scale);

	for(i=0,k=0; i < Nsub;++i){
		//for(i=0;i<lens->NSubstruct;++i){

		r=Rmax*sqrt(ran2(seed));

		do{
			if(lens->sub_alpha == -1.0){
				lens->sub_mass[k] = lens->sub_Mmin
						*pow(lens->sub_Mmax/lens->sub_Mmin,ran2(seed));
			}else{
				lens->sub_mass[k] = lens->sub_Mmin
						*pow(ran2(seed)*(pow(lens->sub_Mmax/lens->sub_Mmin,lens->sub_alpha+1)-1)+1.0
								,1.0/(lens->sub_alpha+1));
			}
		}while(lens->sub_mass[k] < lens->sub_Mmin);  // not sure why this is necessary

		// average density of a substructure does not scale with host
		if(lens->sub_type != pointmass) lens->sub_Rcut[k] = lens->sub_Rmax
				*pow(lens->sub_mass[k]/lens->sub_Mmax,1/3.);

		// maximum radius for a substructure of this mass
		rmax = (host_ro_save*rangeInRei + lens->sub_Rcut[k]
		     + pow(2*lens->sub_mass[k]*host_ro_save/pi/lens->Sigma_crit/sheartol,1./3.) );

		//printf("lens->RcutSubstruct[%i] = %e\n",k,lens->RcutSubstruct[k]);
		//printf("%e %e %e Rmax=%e\n",r/rmax,r,rmax,Rmax);
		if( r < rmax){
			theta=2*pi*ran2(seed);
			lens->sub_x[k][0]=r*cos(theta);
			lens->sub_x[k][1]=r*sin(theta);
			assert(k<NsubMax);

			rav[0] += lens->sub_x[k][0];
			rav[1] += lens->sub_x[k][1];
			r2av += r*r;
			area_av += pow(lens->sub_Rcut[k],2);
			++k;
		}
	}

	lens->sub_N = k;

	//printf("sub_N = %li Nsub = %li ndensity =%e\n",lens->sub_N,Nsub,lens->sub_Ndensity);
/*
	if(lens->sub_N > 2){
		// decide if it is worth doing substructure force calculation by tree or
		//   by direct summation
		r2av /= lens->sub_N;
		r2av = r2av - (rav[0]*rav[0]+rav[1]*rav[1])/(lens->sub_N)/(lens->sub_N);
		area_av /= lens->sub_N;

		//printf("relative area of average particle %e \n",area_av/r2av);

		assert(r2av > 0);
		assert(area_av >= 0.0);

		if(lens->sub_type == pointmass){  // always use tree for point mass substructures
			lens->sub_theta_force=1.0e-1;
			lens->sub_tree=BuildTreeNB(lens->sub_x,lens->sub_Rcut,lens->sub_mass,
				False,True,lens->sub_N,lens->sub_substructures,2,lens->sub_theta_force);
		}else if( area_av/r2av < 1.0e-2 && lens->sub_N > 300 ){
			// build tree for doing substructure force calculation
			lens->sub_theta_force=1.0e-1;
			lens->sub_tree=BuildTreeNB(lens->sub_x,lens->sub_Rcut,lens->sub_mass,
				True,True,lens->sub_N,lens->sub_substructures,2,lens->sub_theta_force);
		}else lens->sub_tree = 0;
	}else lens->sub_tree = 0;
*/
	return;
}


/// for generating a random deviates drawn from approximately the same as the
//      values of table
double RandomFromTable(double *table,unsigned long Ntable,long *seed){
	double y;
	unsigned long j;

	y=ran2(seed)*(Ntable-1);
	j=(int)(y);
	return (table[j+1]-table[j])*(y-j) + table[j];
}

double FractionWithinRe(AnaLens *lens,double rangeInRei){
	double B;

	B = (lens->sub_Rmax/pow(lens->sub_Mmax,1./3.)
			+ pow(2*lens->host_ro/pi/lens->Sigma_crit/1.e-3,1./3.) );

	return 1+(1+lens->sub_alpha)*(
		2*rangeInRei*lens->host_ro*B*(
				pow(lens->sub_Mmax,lens->sub_alpha+4./3.)-pow(lens->sub_Mmin,lens->sub_alpha+4./3.)
				)/(lens->sub_alpha+4./3.)
		+ B*B*(
				pow(lens->sub_Mmax,lens->sub_alpha+5./3.)-pow(lens->sub_Mmin,lens->sub_alpha+5./3.)
				)/(lens->sub_alpha+5./3.)
		)/(pow(rangeInRei*lens->host_ro,2)*(pow(lens->sub_Mmax,lens->sub_alpha+1.0)-pow(lens->sub_Mmin,lens->sub_alpha+1.0))
				);
}

double averageSubMass(AnaLens *lens){
	// average mass of substructures
	return lens->sub_Mmax*(lens->sub_alpha+1)
				  /(lens->sub_alpha+2)*(1-pow(lens->sub_Mmin/lens->sub_Mmax,lens->sub_alpha+2))/
				  (1-pow(lens->sub_Mmin/lens->sub_Mmax,lens->sub_alpha+1));
}
