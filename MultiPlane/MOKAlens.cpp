/*
 * MOKAlens.cpp
 *
 *  Created on: Jun 8, 2012
 *      Author: mpetkova
 */

#ifdef WITH_MOKA

#include <MOKAlens.h>
#include <string>
#include <sstream>
#include <fstream>
#include <utilities.h>

using namespace std;

/*
 * \ingroup Constructor
 * \brief allocates and reads the MOKA map in
 */
MOKALens::MOKALens(std::string paramfile) : Lens(){
	map = new MOKAmap;

	readParamfile(paramfile);

	getDims(MOKA_input_file,&(map->nx),&(map->ny));

	map->convergence.resize(map->nx*map->ny);
	map->alpha1.resize(map->nx*map->ny);
	map->alpha2.resize(map->nx*map->ny);
	map->gamma1.resize(map->nx*map->ny);
	map->gamma2.resize(map->nx*map->ny);
	map->gamma3.resize(map->nx*map->ny);
	readImage(MOKA_input_file
			,&map->convergence
			,&map->alpha1
			,&map->alpha2
			,&map->gamma1
			,&map->gamma2
    		        ,&map->LH);
		  /*
			,&(map->boxl)
			,&(map->boxlMpc)
			,&(map->zlens)
			,&(map->zsource)
			,&(map->omegam)
			,&(map->omegal)
			,&(map->h)
			,&(map->DL));
		  */

	map->boxl = map->LH.boxlarcsec;
	map->zlens = map->LH.zl;
	map->zsource = map->LH.zs;
	map->omegam = map->LH.omegam;
	map->omegal = map->LH.omegal;
	map->DL = map->LH.DL;
	map->center[0] = map->center[1] = 0.0;

	map->boxlMpc = map->LH.boxlMpc;
	map->h = map->LH.h;

	map->boxlMpc /= map->h;

	/// to radians
	map->boxl *= pi/180/3600.;
}

MOKALens::~MOKALens(){
	map->convergence.resize(0);
	map->alpha1.resize(0);
	map->alpha2.resize(0);
	map->gamma1.resize(0);
	map->gamma2.resize(0);
	map->gamma3.resize(0);
	delete map;
}

/*
 * sets the cosmology and the lens and the source according to the MOKA map parameters
 */
void MOKALens::setInternalParams(CosmoHndl cosmo, SourceHndl source){
	cosmo->setOmega_matter(map->LH.omegam,true);
	cosmo->sethubble(map->LH.h);
	setZlens(map->LH.zl);
	source->zsource = map->LH.zs;

	double Ds = cosmo->angDist(0,map->LH.zs);
	double Dl = cosmo->angDist(0,map->LH.zl);
	double Dls = cosmo->angDist(map->LH.zl,map->LH.zs);
	double fac = Ds/Dls/Dl;

	/// converts to the code units
	int i, j;
	for(i=0;i<map->nx;i++)
		for(j=0;j<map->ny;j++){
			int index = i+map->ny*j;
			map->convergence[index] *= fac;
			map->gamma1[index] *= fac;
			map->gamma2[index] *= fac;
		}
}


/** \ingroup ImageFinding
 * \brief Reads in a parameter file and sets up a MOKA lens map.
 *
 * Sets many parameters within the MOKA lens model
 */

void MOKALens::readParamfile(std::string filename){
  const int MAXPARAM = 50;
  string label[MAXPARAM], rlabel, rvalue;
  void *addr[MAXPARAM];
  int id[MAXPARAM];
  std::stringstream ss;
  int i, n;
  int myint;
  double mydouble;
  string mystring;
  string escape = "#";
  char dummy[300];
  int flag;

  n = 0;

  // id[] = 0 double, 1 int, 2 string

  addr[n] = &zlens;
  id[n] = 0;
  label[n++] = "z_lens";

  addr[n] = &MOKA_input_file;
  id[n] = 2;
  label[n++] = "MOKA_input_file";

  std::ifstream file_in(filename.c_str());
  if(!file_in){
    cout << "Can't open file " << filename << endl;
    exit(1);
  }

  // output file
  while(!file_in.eof()){
	  file_in >> rlabel >> rvalue;
	  file_in.getline(dummy,100);

	  if(rlabel[0] == escape[0])
		  continue;

	  flag = 0;

	  for(i = 0; i < n; i++){
		  if(rlabel == label[i]){

			  flag = 1;
			  ss << rvalue;

			  switch(id[i]){
			  case 0:
				  ss >> mydouble;
				  *((double *)addr[i]) = mydouble;
				  break;
			  case 1:
				  ss >> myint;
				  *((int *)addr[i]) = myint;
				  break;
			  case 2:
				  ss >> mystring;
				  *((string *)addr[i]) = mystring;
				  break;
			  }

			  ss.clear();
			  ss.str(std::string());

			  id[i] = -1;
		  }
	  }
  }

  for(i = 0; i < n; i++){
	  if(id[i] > 0){
		  ERROR_MESSAGE();
		  cout << "parameter " << label[i] << " needs to be set!" << endl;
		  exit(0);
	  }
  }

  file_in.close();

  set = true;

}


double MOKALens::getZlens(){
	return zlens;
}

void MOKALens::setZlens(double z){
	zlens = z;
}

/*
 * saves the image, by rading off the calues from the image tree
 * and then saving to a fits file and computing the radial profile
 * of the convergence
 */
void MOKALens::saveImage(GridHndl grid,bool saveprofiles){
	std::stringstream f;
	std::string filename;

	f << "noisy_" << MOKA_input_file;
	filename = f.str();

	MoveToTopList(grid->i_tree->pointlist);

	do{
		long index = IndexFromPosition(grid->i_tree->pointlist->current->x,map->nx,map->boxl,map->center);
		if(index > -1){

			map->convergence[index] = grid->i_tree->pointlist->current->kappa;
			map->gamma1[index] = grid->i_tree->pointlist->current->gamma[0];
			map->gamma2[index] = grid->i_tree->pointlist->current->gamma[1];
			map->gamma3[index] = grid->i_tree->pointlist->current->gamma[2];
		}
	}while(MoveDownList(grid->i_tree->pointlist)==true);

	map->boxl *= 180/pi*3600;

	writeImage(filename
			,map->convergence
			,map->gamma1
			,map->gamma2
			,map->gamma3
			,map->nx
			,map->ny
   		        ,map->LH);
	/*
			,map->boxl
			,map->boxlMpc
			,map->zlens
			,map->zsource
			,map->omegam
			,map->omegal
			,map->h
			,map->DL);
	*/

	if(saveprofiles == true){
	  std:: cout << " saving profile " << std:: endl;
	            saveProfiles();
		 // saveKappaProfile();
		 // saveGammaProfile();
	}
}

/*
 * computing and saving the radial profile of the convergence
 */
void MOKALens::saveKappaProfile(){

	/* measuring the differential and cumulative profile*/
	double xmin = -map->boxlMpc*0.5*map->h;
	double xmax =  map->boxlMpc*0.5*map->h;
	double drpix = map->boxlMpc/map->nx*map->h;
	std::valarray<float> pxdist(map->nx*map->ny);
	for(int i=0; i<map->nx; i++ ) for(int j=0; j<map->ny; j++ ){
		pxdist[i+map->ny*j]= sqrt(pow((xmin+(drpix*0.5)+i*drpix),2) +
				pow((xmin+(drpix*0.5)+j*drpix),2));
	}
	double dr0 = 8.*(0.5*map->boxlMpc*map->h)/(map->nx/2.);
	int nbin = int(0.5*map->boxlMpc*map->h/dr0);
	//
	std:: cout << "   " << std:: endl;
	std:: cout << " nbins = " << nbin << "  dr0 = " << dr0 << std:: endl;
	std:: cout << " ______________________________________________________ " << std:: endl;
	std:: cout << " computing profiles assuming spherical symmetry";
	// - - - - - - - - - - - - - - - - -
	double *kprofr = estprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax);
	double *sigmakprof = estsigmaprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax,kprofr);
	double *ckprofr = estcprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax);
	double *sigmackprof = estsigmacprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax,kprofr);
	std::ostringstream fprof;
	fprof << "MAP_radial_prof_kappa_" << MOKA_input_file << ".dat";
	std:: ofstream filoutprof;
	std:: string filenameprof = fprof.str();
	filoutprof.open(filenameprof.c_str());
	filoutprof <<"# r      kappa     sig_k     ckappa     sig_ck" << std:: endl;
	int l;
	for(l=0;l<nbin;l++){
		filoutprof << dr0*l + dr0/2. << "  "
				<< kprofr[l] << "  " << sigmakprof[l] << "  "
				<< ckprofr[l] << "  " << sigmackprof[l] <<
				std:: endl;
	}
	filoutprof.close();
}

void MOKALens::saveGammaProfile(){

	/* measuring the differential and cumulative profile*/
	double xmin = -map->boxlMpc*0.5*map->h;
	double xmax =  map->boxlMpc*0.5*map->h;
	double drpix = map->boxlMpc/map->nx*map->h;
	std::valarray<float> pxdist(map->nx*map->ny);
	for(int i=0; i<map->nx; i++ ) for(int j=0; j<map->ny; j++ ){
		pxdist[i+map->ny*j]= sqrt(pow((xmin+(drpix*0.5)+i*drpix),2) +
				pow((xmin+(drpix*0.5)+j*drpix),2));
	}
	double dr0 = 8.*(0.5*map->boxlMpc*map->h)/(map->nx/2.);
	int nbin = int(0.5*map->boxlMpc*map->h/dr0);
	//
	std:: cout << "   " << std:: endl;
	std:: cout << " nbins = " << nbin << "  dr0 = " << dr0 << std:: endl;
	std:: cout << " ______________________________________________________ " << std:: endl;
	std:: cout << " computing profiles assuming spherical symmetry";
	// - - - - - - - - - - - - - - - - -
	double *kprofr = estprof(map->gamma2,map->nx,map->ny,pxdist,dr0,xmax);
	double *sigmakprof = estsigmaprof(map->gamma2,map->nx,map->ny,pxdist,dr0,xmax,kprofr);
	double *ckprofr = estcprof(map->gamma2,map->nx,map->ny,pxdist,dr0,xmax);
	double *sigmackprof = estsigmacprof(map->gamma2,map->nx,map->ny,pxdist,dr0,xmax,kprofr);
	std::ostringstream fprof;
	fprof << "MAP_radial_prof_gamma_" << MOKA_input_file << ".dat";
	std:: ofstream filoutprof;
	std:: string filenameprof = fprof.str();
	filoutprof.open(filenameprof.c_str());
	filoutprof <<"# r      kappa     sig_k     ckappa     sig_ck" << std:: endl;
	int l;
	for(l=0;l<nbin;l++){
		filoutprof << dr0*l + dr0/2. << "  "
				<< kprofr[l] << "  " << sigmakprof[l] << "  "
				<< ckprofr[l] << "  " << sigmackprof[l] <<
				std:: endl;
	}
	filoutprof.close();
}

void MOKALens::saveProfiles(){
	/* measuring the differential and cumulative profile*/
	double xmin = -map->boxlMpc*0.5*map->h;
	double xmax =  map->boxlMpc*0.5*map->h;
	double drpix = map->boxlMpc/map->nx*map->h;
	std:: vector<double> x(map->nx);	
	fill_linear (x,map->nx,xmin,xmax); // physical
	std::valarray<float> pxdist(map->nx*map->ny);
	std::valarray<float> red_sgE(map->nx*map->ny),red_sgB(map->nx*map->ny),sgm(map->nx*map->ny); 
	for(int i=0; i<map->nx; i++ ) for(int j=0; j<map->ny; j++ ){
		pxdist[i+map->ny*j]= sqrt(pow((xmin+(drpix*0.5)+i*drpix),2) +
				pow((xmin+(drpix*0.5)+j*drpix),2));
		// reduced shear E a B
		double dx=x[i];
		double dy=x[j];
		double p=atan2( dy, dx ); // check gamma 1 and gamma 2 definition
		red_sgE[i+map->ny*j] = (-map->gamma1[i+map->ny*j]*cos(2*p)-map->gamma2[i+map->ny*j]*sin(2*p))/(1.-map->convergence[i+map->ny*j]);
		red_sgB[i+map->ny*j] = (map->gamma1[i+map->ny*j]*sin(2*p)-map->gamma2[i+map->ny*j]*cos(2*p))/(1.-map->convergence[i+map->ny*j]);
		sgm[i+map->ny*j] = sqrt(pow(map->gamma1[i+map->ny*j],2) + pow(map->gamma2[i+map->ny*j],2));
	}

	double dr0 = 8.*(0.5*map->boxlMpc*map->h)/(map->nx/2.);
	int nbin = int(xmax/dr0);                           
	double inarcsec  = 10800./M_PI/map->LH.DL*60.;                             
	//                                                                                           
	std:: cout << "   " << std:: endl;                                                           
	std:: cout << " nbins = " << nbin << "  dr0 = " << dr0 << std:: endl;                        
	std:: cout << " ______________________________________________________ " << std:: endl;      
	std:: cout << " computing profiles assuming spherical symmetry";                                
	// - - - - - - - - - - - - - - - - -                                                         
	double *kprofr = estprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax);                                          
	double *sigmakprof = estsigmaprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax,kprofr);                          
	double *ckprofr = estcprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax);                                          
	double *sigmackprof = estsigmacprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax,kprofr);                          
	double *gamma1profr = estprof(red_sgE,map->nx,map->ny,pxdist,dr0,xmax); // reduced shear
	double *sigmagamma1prof = estsigmaprof(red_sgE,map->nx,map->ny,pxdist,dr0,xmax,gamma1profr);              
	double *gamma0profr = estprof(red_sgB,map->nx,map->ny,pxdist,dr0,xmax);  
	double *sigmagamma0prof = estsigmaprof(red_sgB,map->nx,map->ny,pxdist,dr0,xmax,gamma0profr);              
	double *gamma2profr = estprof(sgm,map->nx,map->ny,pxdist,dr0,xmax);  
	double *sigmagamma2prof = estsigmaprof(sgm,map->nx,map->ny,pxdist,dr0,xmax,gamma2profr);              
	std::ostringstream fprof;
	fprof << "MAP_radial_prof_" << MOKA_input_file << ".dat";
	std:: ofstream filoutprof;
	std:: string filenameprof = fprof.str();
	filoutprof.open(filenameprof.c_str());
	filoutprof <<"# r      kappa     sig_k     ckappa     sig_ck    redgE   sig_redgE   redgB   sig_redgE   g   sig_g   theta   Anulus_area" << std:: endl;
	for(int l=0;l<nbin;l++){
	  double Aanulus = M_PI*((dr0*l+dr0)*(dr0*l+dr0)-(dr0*l)*(dr0*l)); 
	  filoutprof << dr0*l + dr0/2. << "  " 
		     << kprofr[l] << "  " << sigmakprof[l] << "  " 
		     << ckprofr[l] << "  " << sigmackprof[l] << "   " 
		     << gamma1profr[l] << "  " << sigmagamma1prof[l] << "  " 
		     << gamma0profr[l] << "  " << sigmagamma0prof[l] << "  " 
	    	     << gamma2profr[l] << "  " << sigmagamma2prof[l] << "   "
		     << (dr0*l + dr0/2.)*inarcsec << "   "  << Aanulus*inarcsec*inarcsec << "   " <<  
	    std:: endl;
	}
	filoutprof.close();
	/*
	delete[] kprofr;
	delete[] sigmakprof;
	delete[] ckprofr;
	delete[] sigmackprof;
	delete[] gamma0profr;
	delete[] sigmagamma0prof;
	delete[] gamma1profr;
	delete[] sigmagamma1prof;
	delete[] gamma2profr;
	delete[] sigmagamma2prof;
	*/
}

/** \ingroup DeflectionL2
   *
   * \brief Routine for obtaining the deflection and other lensing quantities for
   * a MOKA map (MOKALens), for just one ray!!
   *
*/
void MOKALens::rayshooterInternal(double *xx, double *alpha, double *gamma, double *kappa, bool kappa_off){

	long index = IndexFromPosition(xx,map->nx,map->boxlMpc,map->center);

	if(index > -1){
		alpha[0] = map->alpha1[index];
		alpha[1] = map->alpha2[index];
		gamma[0] = map->gamma1[index];
		gamma[1] = map->gamma2[index];
		gamma[2] = 0.0;
		*kappa = map->convergence[index];
	}
	else{
		alpha[0] = alpha[1] = 0.0;
		gamma[0] = gamma[1] = gamma[2] = 0.0;
		*kappa = 0.0;
	}


	return;
}
#endif
