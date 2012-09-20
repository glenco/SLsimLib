/*
 * MOKAlens.cpp
 *
 *  Created on: Jun 8, 2012
 *      Author: mpetkova
 */


#include <MOKAlens.h>
#include <string>
#include <sstream>
#include <fstream>
#include <utilities.h>

using namespace std;

/**
 * \brief allocates and reads the MOKA map in
 */
MOKALens::MOKALens(std::string paramfile) : Lens(){
#ifndef ENABLE_FITS
	std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
	exit(1);
#endif

	map = new MOKAmap;
	LH = new LensHalo;

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
			,LH);


	int i,j;
	if(flag_background_field==1){
	  for(i=0;i<map->nx;i++) for(j=0;j<map->ny;j++){
	      map->convergence[i+map->nx*j] = 0;
	      map->alpha1[i+map->nx*j] = 0;
	      map->alpha2[i+map->nx*j] = 0;
	      map->gamma1[i+map->nx*j] = 0;
	      map->gamma2[i+map->nx*j] = 0;
	    }
	}
	initMap();
}

/**
 * \ingroup Constructor
 * \brief allocates and reads the MOKA map in
 */
MOKALens::MOKALens(std::string paramfile,LensHalo *halo) : Lens(){
	map = new MOKAmap;
	LH = halo;

	readParamfile(paramfile);

	map->nx = map->ny = LH->npix;

	map->convergence.resize(map->nx*map->ny);
	map->alpha1.resize(map->nx*map->ny);
	map->alpha2.resize(map->nx*map->ny);
	map->gamma1.resize(map->nx*map->ny);
	map->gamma2.resize(map->nx*map->ny);
	map->gamma3.resize(map->nx*map->ny);

	initMap();
}

MOKALens::~MOKALens(){
	map->convergence.resize(0);
	map->alpha1.resize(0);
	map->alpha2.resize(0);
	map->gamma1.resize(0);
	map->gamma2.resize(0);
	map->gamma3.resize(0);
	delete map;
	delete LH;
}

void MOKALens::initMap(){
	map->boxl = LH->boxlarcsec;
	map->zlens = LH->zl;
	map->zsource = LH->zs;
	map->omegam = LH->omegam;
	map->omegal = LH->omegal;
	map->DL =LH->DL;
	map->center[0] = map->center[1] = 0.0;

	map->boxlMpc = LH->boxlMpc;
	map->h = LH->h;

	map->boxlMpc /= map->h;

	/// to radians
	map->boxl *= pi/180/3600.;

	double xmin = -map->boxlMpc*0.5*map->h;
	double xmax =  map->boxlMpc*0.5*map->h;
	fill_linear (map->x,map->nx,xmin,xmax); // physical
	map->inarcsec  = 10800./M_PI/LH->DL*60.;
}

/** \brief sets the cosmology and the lens and the source according to the MOKA map parameters
 */
void MOKALens::setInternalParams(CosmoHndl cosmo, SourceHndl source){
	cosmo->setOmega_matter(map->omegam,true);
	cosmo->sethubble(map->h);
	setZlens(map->zlens);
	source->setZ(map->zsource);

	double Ds = cosmo->angDist(0,map->zsource);
	double Dl = cosmo->angDist(0,map->zlens);
	double Dls = cosmo->angDist(map->zlens,map->zsource);
	double fac = Ds/Dls/Dl;

	/// converts to the code units
	if(flag_MOKA_analyze == 0){
	  int i, j;
	  for(i=0;i<map->nx;i++)
	    for(j=0;j<map->ny;j++){
	      int index = i+map->ny*j;
	      map->convergence[index] *= fac;
	      map->gamma1[index] *= fac;
	      map->gamma2[index] *= fac;
	    }
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

  addr[n] = &flag_MOKA_analyze;
  id[n] = 1;
  label[n++] = "MOKA_analyze";

  addr[n] = &flag_background_field;
  id[n] = 1;
  label[n++] = "background_field";

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
	  if(id[i] > 0 && addr[i] != &flag_MOKA_analyze){
		  ERROR_MESSAGE();
		  cout << "parameter " << label[i] << " needs to be set!" << endl;
		  exit(0);
	  }

	  if(id[i] >= 0 && addr[i] == &flag_MOKA_analyze){
		  flag_MOKA_analyze = 0; //false, no analyzis, prepare for ray-shooting
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

/**
 * saves the image, by rading off the calues from the image tree
 * and then saving to a fits file and computing the radial profile
 * of the convergence
 */
void MOKALens::saveImage(GridHndl grid,bool saveprofiles){
	std::stringstream f;
	std::string filename;

	if(flag_background_field==1) f << MOKA_input_file << "_only_noise.fits";
	else f << MOKA_input_file << "_noisy.fits";
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
		    ,LH);

	if(saveprofiles == true){

	  std:: cout << " saving profile " << std:: endl;
	  double RE3;
	  saveProfiles(RE3);
	  estSignLambdas();
	  double RE1,RE2;
	  EinsteinRadii(RE1,RE2);
	  std::ostringstream fEinr;
	  if(flag_background_field==1) fEinr << MOKA_input_file << "_only_noise_Einstein.radii.dat";
	  else fEinr << MOKA_input_file << "_noisy_Einstein.radii.dat";
	  std:: ofstream filoutEinr;
	  std:: string filenameEinr = fEinr.str();
	  filoutEinr.open(filenameEinr.c_str());
	  filoutEinr << "# effective        median      from_profiles" << std:: endl;
	  filoutEinr << RE1 << "   " << RE2 << "    " << RE3 << std:: endl;
	  filoutEinr.close();
	}
}


/**
 * computing and saving the radial profile of the convergence, reduced tangential and parallel shear and of the shear
 *  */
void MOKALens::saveProfiles(double &RE3){
	/* measuring the differential and cumulative profile*/
	double xmin = -map->boxlMpc*0.5*map->h;
	double xmax =  map->boxlMpc*0.5*map->h;
	double drpix = map->boxlMpc/map->nx*map->h;

	std::valarray<float> pxdist(map->nx*map->ny);
	std::valarray<float> red_sgE(map->nx*map->ny),red_sgB(map->nx*map->ny),sgm(map->nx*map->ny); 
	int i, j;
	for(i=0; i<map->nx; i++ ) for(j=0; j<map->ny; j++ ){
		pxdist[i+map->ny*j]= sqrt(pow((xmin+(drpix*0.5)+i*drpix),2) +
				pow((xmin+(drpix*0.5)+j*drpix),2));
		// reduced shear E a B
		double dx=map->x[i];
		double dy=map->x[j];
		double p=atan2( dy, dx ); // check gamma 1 and gamma 2 definition
		red_sgE[i+map->ny*j] = (-map->gamma1[i+map->ny*j]*cos(2*p)-map->gamma2[i+map->ny*j]*sin(2*p))/(1.-map->convergence[i+map->ny*j]);
		red_sgB[i+map->ny*j] = (map->gamma1[i+map->ny*j]*sin(2*p)-map->gamma2[i+map->ny*j]*cos(2*p))/(1.-map->convergence[i+map->ny*j]);
		sgm[i+map->ny*j] = sqrt(pow(map->gamma1[i+map->ny*j],2) + pow(map->gamma2[i+map->ny*j],2));
	}

	double dr0 = 8.*(0.5*map->boxlMpc*map->h)/(map->nx/2.);
	int nbin = int(xmax/dr0);                           
	int nbggal=0;
	int ih;
	// check if backgroundgalXarcmin2.d file exist
	std:: string filebgXarcmin2 = "backgroundgalXarcmin2.d";
	std:: ifstream infilebgXarcmin2;
	infilebgXarcmin2.open(filebgXarcmin2.c_str());
	// ... if so it read it!
	if(infilebgXarcmin2.is_open()){
	  std:: cout << " I will read  backgroundgalXarcmin2.d file" << std:: endl;
	  infilebgXarcmin2 >> nbggal;
	  infilebgXarcmin2 >> ih;
	  std:: cout << "  and consider " << nbggal <<  " per arcmin2 as background points " << std:: endl;
	  std:: cout << " in building the cluster profile" << std:: endl;
	  infilebgXarcmin2.close();
	}    
	double dntbggal=double(nbggal)*(map->boxl*map->boxl)/3600.;
	int ntbggal=int(dntbggal+0.5);
	std:: vector<int> runi,runj;
	int ibut;
	int lrun;
	if(ntbggal>0){
	  std:: cout << " ~ ~ ~ " << ntbggal << " number of bg gal in the field" << std:: endl;
	  // fill the vector with random coordinates of background galaxies
	  srand(ih+94108);
	  for(lrun=0;lrun<ntbggal;lrun++){
	    ibut = rand()%map->nx;
	    if(ibut>map->nx || ibut<0) std:: cout << "1. ibut out of npix range = " << ibut << std:: endl;
	    runi.push_back(ibut);
	    ibut = rand()%map->nx;
	    if(ibut>map->nx || ibut<0) std:: cout << "2. ibut out of npix range = " << ibut << std:: endl;
	    runj.push_back(ibut);
	  }
	  if(runi.size()>ntbggal || runj.size()>ntbggal){
	    std:: cout << " random points in the field out of range " << std:: endl;
	    std:: cout << " runi.size() = " << runi.size() << std:: endl;
	    std:: cout << " runj.size() = " << runj.size() << std:: endl;	      
	  }
	  dr0 = nXbin*(0.5*map->boxlMpc*map->h)/(0.5*ntbggal);
	  nbin = int(xmax/dr0);                                                        
	}
	//                                                                                           
	std:: cout << "   " << std:: endl;                                                           
	std:: cout << " nbins = " << nbin << "  dr0 = " << dr0 << std:: endl;                        
	std:: cout << " ______________________________________________________ " << std:: endl;      
	std:: cout << " computing profiles assuming spherical symmetry";                                
	// - - - - - - - - - - - - - - - - -                                                         
	double *kprofr = estprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal);                                          
	double *sigmakprof = estsigmaprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal,kprofr);                          
	double *ckprofr = estcprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal);                                          
	double *sigmackprof = estsigmacprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal,kprofr);                          
	double *gamma1profr = estprof(red_sgE,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal); // reduced shear
	double *sigmagamma1prof = estsigmaprof(red_sgE,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal,gamma1profr);              
	double *gamma0profr = estprof(red_sgB,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal);  
	double *sigmagamma0prof = estsigmaprof(red_sgB,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal,gamma0profr);              
	double *gamma2profr = estprof(sgm,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal);  
	double *sigmagamma2prof = estsigmaprof(sgm,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal,gamma2profr);              
	std::ostringstream fprof;

	if(flag_background_field==1) fprof << MOKA_input_file << "_only_noise_MAP_radial_prof.dat";
	else fprof << MOKA_input_file << "_noisy_MAP_radial_prof.dat";
	std:: ofstream filoutprof;
	std:: string filenameprof = fprof.str();
	filoutprof.open(filenameprof.c_str());
	filoutprof <<"# r      kappa     sig_k     ckappa     sig_ck    redgE   sig_redgE   redgB   sig_redgE   g   sig_g   theta   Anulus_area" << std:: endl;
	int l;
	std:: vector<double> lrOFr(nbin),lprofFORre(nbin);
	for(l=0;l<nbin;l++){
	  double Aanulus = M_PI*((dr0*l+dr0)*(dr0*l+dr0)-(dr0*l)*(dr0*l)); 
	  filoutprof << dr0*l + dr0/2. << "  " 
		     << kprofr[l] << "  " << sigmakprof[l] << "  " 
		     << ckprofr[l] << "  " << sigmackprof[l] << "   " 
		     << gamma1profr[l] << "  " << sigmagamma1prof[l] << "  " 
		     << gamma0profr[l] << "  " << sigmagamma0prof[l] << "  " 
	    	     << gamma2profr[l] << "  " << sigmagamma2prof[l] << "   "
		     << (dr0*l + dr0/2.)*map->inarcsec << "   "  << Aanulus*map->inarcsec*map->inarcsec << "   " <<  
	    std:: endl;
	  lprofFORre[l] = log10(kprofr[l] + gamma2profr[l]);
	  lrOFr[l] = log10(dr0*l + dr0/2.);
	}
	filoutprof.close();
	RE3 = InterpolateYvec(lprofFORre,lrOFr,0.);  
	RE3 = pow(10.,RE3)*map->inarcsec;
}

/** \ingroup DeflectionL2
   *
   * \brief Routine for obtaining the deflection and other lensing quantities for
   * a MOKA map (MOKALens), for just one ray!!
   *
*/
void MOKALens::rayshooterInternal(double *xx, double *alpha, float *gamma, float *kappa, bool kappa_off){

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

/**
 * compute the signal of \lambda_r and \lambda_t
 */
void MOKALens::estSignLambdas(){
  map->Signlambdar.resize(map->nx*map->ny);
  map->Signlambdat.resize(map->nx*map->ny);
  double gamma,lambdar,lambdat;
  int i, j;
  for(i=0;i<map->nx;i++)
    for(j=0;j<map->ny;j++){
      gamma = sqrt(pow(map->gamma1[i+map->ny*j],2) + 
			  pow(map->gamma2[i+map->ny*j],2));
      lambdat=1-map->convergence[i+map->ny*j]-gamma;
      lambdar=1-map->convergence[i+map->ny*j]+gamma;
      
      if(lambdar>=0) map->Signlambdar[i+map->ny*j]=1;
      else map->Signlambdar[i+map->ny*j]=-1;
      
      if(lambdat>=0) map->Signlambdat[i+map->ny*j]=1;
      else map->Signlambdat[i+map->ny*j]=-1;     
    }
}

/**
 * measure the effective and the median Einstein radii of the connected critical 
 * points present at the halo center
 */
void MOKALens::EinsteinRadii(double &RE1, double &RE2){
  double signV;
  //  std:: vector<double> xci1,yci1;
  std:: vector<double> xci2,yci2;
  // open file readable by ds9
  std::ostringstream fcrit;
  if(flag_background_field==1) fcrit << MOKA_input_file << "_only_noise_Criticals.reg";
  else fcrit << MOKA_input_file << "_noisy_Criticals.reg";
  std:: ofstream filoutcrit;
  std:: string filenamecrit = fcrit.str();
  filoutcrit.open(filenamecrit.c_str());
  // define the critical points in the map
  int i, j;
  for(i=1;i<map->nx-1;i++)
    for(j=1;j<map->ny-1;j++){
      signV=map->Signlambdar[i-1+map->ny*j]+map->Signlambdar[i+map->ny*(j-1)]+
	map->Signlambdar[i+1+map->ny*j]+map->Signlambdar[i+map->ny*(j+1)];      
      if(fabs(signV)<4.){
	// xci1.push_back(map->x[i]);
	// yci1.push_back(map->x[j]);
	filoutcrit << "circle(" << i << "," << j << ",0.5)" << std:: endl;
      }
      signV=map->Signlambdat[i-1+map->ny*j]+map->Signlambdat[i+map->ny*(j-1)]+
	map->Signlambdat[i+1+map->ny*j]+map->Signlambdat[i+map->ny*(j+1)];      
      if(fabs(signV)<4.){
	xci2.push_back(map->x[i]);
	yci2.push_back(map->x[j]);
	filoutcrit << "circle(" << i << "," << j << ",0.5)" << std:: endl;
      }
    }
  filoutcrit.close();
  double pixDinL = map->boxlMpc*map->h/double(map->nx);
  /* measure the Einstein radius */
  std:: vector<double> xci,yci;	
  //for(int ii=0;ii<xci1.size();ii++){
  //  xci.push_back(xci1[ii]);
  //  yci.push_back(yci1[ii]);
  //}
  for(int ii=0;ii<xci2.size();ii++){
    xci.push_back(xci2[ii]);
    yci.push_back(yci2[ii]);
  }
  // xci1.clear();
  // yci1.clear();
  xci2.clear();
  yci2.clear();
  int nc = xci.size();
  std:: vector<int> groupid(nc);
  int largestgroupid = fof(pixDinL,xci,yci,groupid);
  std:: vector<double> xcpoints,ycpoints;
  double xercm=0;
  double yercm=0;
  for(int ii=0;ii<nc;ii++){
    if(groupid[ii] == largestgroupid){
      xcpoints.push_back(xci[ii]);
      ycpoints.push_back(yci[ii]);
      xercm+=xci[ii];
      yercm+=yci[ii];
    }
  }
  nc = xcpoints.size();
  xercm=xercm/double(nc);
  yercm=yercm/double(nc);
  if(nc>0){
    std:: vector<double>::iterator maxit, minit; 
    // find the min and max elements in the vector
    maxit = max_element(xcpoints.begin(), xcpoints.end());
    minit = min_element(xcpoints.begin(), xcpoints.end());
    double xmincpoints,xmaxcpoints;
    xmaxcpoints = *maxit;
    xmincpoints = *minit;
    int imin = locate(map->x,xmincpoints);
    imin = ((imin > 0) ? imin:0);
    imin = ((imin < map->nx-1) ? imin:map->nx-1);
    // imin=GSL_MIN( GSL_MAX( imin, 0 ), map->nx-1 );
    int imax = locate(map->x,xmaxcpoints);
    imax = ((imax > 0) ? imax:0);
    imax = ((imax < map->nx-1) ? imax:map->nx-1);
    // imax=GSL_MIN( GSL_MAX( imax, 0 ), map->nx-1 );
    std:: vector<double> ysup,yinf,xsup,xinf;
    for(int ii=imin;ii<=imax;ii++){
      std:: vector<double>ybut;
      int condition=0;
      for(int ji=0;ji<nc;ji++){
	if(fabs(xcpoints[ji]-map->x[ii])<pixDinL/2){
	  if(condition==0){
	    xsup.push_back(xcpoints[ji]);
	    xinf.push_back(xcpoints[ji]);
	    condition=1;
	  }
	  ybut.push_back(ycpoints[ji]);
	}
      }
      if(ybut.size()>0){
	std:: vector<double>::iterator ymax, ymin; 
	// Find the min and max elements in the vector
	ymax = max_element(ybut.begin(), ybut.end());
	ymin = min_element(ybut.begin(), ybut.end());
	double ymincpoints,ymaxcpoints;
	ymaxcpoints = *ymax;
	ymincpoints = *ymin;  
	ysup.push_back(ymaxcpoints);
	yinf.push_back(ymincpoints);
      }  
      if(ybut.size()==1){
	double ymincpoints,ymaxcpoints;
	ymaxcpoints = ybut[0];
	ymincpoints = ybut[0];
	ysup.push_back(ymaxcpoints);
	yinf.push_back(ymincpoints);
      }  
    }
    nc = yinf.size();
    int npixIN=0;
    std:: vector<double> RE;
    for(int ii=0;ii<nc;ii++){
      RE.push_back(sqrt(pow(xinf[ii]-xercm,2.) + pow(yinf[ii]-yercm,2)));
      std:: vector<double> ycounts;
      for(int ji=0;ji<map->nx;ji++){
	if(map->x[ji]>=yinf[ii] && map->x[ji]<=ysup[ii]){
	  ycounts.push_back(map->x[ji]);
	}
      }
      int ncounts = ycounts.size();
      npixIN=npixIN+ncounts;
    }
    for(int ii=nc-1;ii>=0;ii--){
      RE.push_back(sqrt(pow(xsup[ii]-xercm,2) + pow(ysup[ii]-yercm,2)));
    }
    RE1=map->inarcsec*sqrt(pixDinL*pixDinL*npixIN/M_PI);
    RE2=map->inarcsec*median(RE);
    if(RE2!=RE2) RE2=0.;
  }
  else{
    RE1=0.;
    RE2=0.;
  }
}
/**
 * saves MAP properties, computing the radial profile
 * of the convergence and shear
 */
void MOKALens::saveImage(bool saveprofiles){
        std::stringstream f;
        std::string filename;  
	if(flag_background_field==1) f << MOKA_input_file << "_only_noise.fits";
	else f << MOKA_input_file << "_noisy.fits";
	filename = f.str();
	/*
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
	*/
	map->boxl *= 180/pi*3600;

	writeImage(filename
		   ,map->convergence
		   ,map->gamma1
		   ,map->gamma2
		   ,map->gamma3
		   ,map->nx
		   ,map->ny
		   ,LH);
	
	if(saveprofiles == true){
	  std:: cout << " saving profile " << std:: endl;
                    double RE3;
	            saveProfiles(RE3);
		    estSignLambdas(); 
		    double RE1,RE2;
		    EinsteinRadii(RE1,RE2);
		    std::ostringstream fEinr;
		    if(flag_background_field==1) fEinr << MOKA_input_file << "_only_noise_Einstein.radii.dat";
		    else fEinr << MOKA_input_file << "_noisy_Einstein.radii.dat";
		    std:: ofstream filoutEinr;
		    std:: string filenameEinr = fEinr.str();
		    filoutEinr.open(filenameEinr.c_str());
		    filoutEinr << "# effective        median      from_profles" << std:: endl;
		    filoutEinr << RE1 << "   " << RE2 << "    " << RE3 << std:: endl;
		    filoutEinr.close();
		    // saveKappaProfile();
		    // saveGammaProfile();
	}
}

