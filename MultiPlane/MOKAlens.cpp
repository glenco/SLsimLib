/*
 * MOKAlens.cpp
 *
 *  Created on: Jun 8, 2012
 *      Author: mpetkova
 */


#include "slsimlib.h"

using namespace std;

/*
 * \brief center of mass of a map using the moving centre
 */
void cmass(int n, std::valarray<float> map, std:: vector<double> x, double &xcm, double &ycm){ 
  double shrink = 1.02;
  // 0.05% of the total number of pixels
  int nstop = int(0.05*double(n)*double(n)/100.);
  if(nstop<36) nstop = 36;
  std:: cout << "nstop = " << nstop << std:: endl;
  xcm = 0.;
  ycm = 0.;
  double tot = 0;
  //
  double rsel = x[n-1] - x[0];
  double xmin = x[0];
  double drpix = rsel/n;  
  rsel/=2.;
  // I will set in the centre with t width equal to 128 pixels
  rsel = drpix*64.;
  //
  int nsel = n*n;
  int i,j;
  double dist;
  while(nsel>nstop){
    double nxcm = 0;
    double nycm = 0;
    double xi,yi;
    tot = 0;
    nsel = 0;
    for(i=0;i<n;i++) for(j=0;j<n;j++){
	xi = (xmin+(drpix*0.5)+i*drpix);
	yi = (xmin+(drpix*0.5)+j*drpix);
	dist = sqrt(pow(xi-xcm,2)+ pow(yi-ycm,2));
	if(dist<=rsel){
	  nxcm+=map[i+n*j]*xi;
	  nycm+=map[i+n*j]*yi;
	  tot+=map[i+n*j];
	  nsel++;
	}
      }
    nxcm/=tot;
    nycm/=tot;
    rsel = rsel / shrink;
    // center of mass
    xcm = nxcm;
    ycm = nycm;
  }
}

/**
 * \brief allocates and reads the MOKA map in
 */
//MOKALens::MOKALens(std::string paramfile) : Lens(){
MOKALens::MOKALens(InputParams& params) : Lens(){
#ifndef ENABLE_FITS
	std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
	exit(1);
#endif

	map = new MOKAmap;
	LH = new LensHalo;

	assignParams(params);

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
//MOKALens::MOKALens(std::string paramfile,LensHalo *halo) : Lens(){
MOKALens::MOKALens(InputParams& params,LensHalo *halo) : Lens(){
	map = new MOKAmap;
	LH = halo;

	assignParams(params);

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
	map->boxlarcsec = LH->boxlarcsec;
	map->zlens = LH->zl;
	map->zsource = LH->zs;
	map->omegam = LH->omegam;
	map->omegal = LH->omegal;
	map->DL =LH->DL;
	map->center[0] = map->center[1] = 0.0;
	map->boxlMpc = LH->boxlMpc;
	map->h = LH->h;
	map->boxlrad = map->boxlarcsec*pi/180/3600.;

	double xmin = -map->boxlMpc*0.5;
	double xmax =  map->boxlMpc*0.5;
	fill_linear (map->x,map->nx,xmin,xmax); // physical Mpc/h
	map->inarcsec  = 10800./M_PI/LH->DL*60.; // Mpc/h to arcsec
}

/** \brief sets the cosmology and the lens and the source according to the MOKA map parameters
 */
void MOKALens::setInternalParams(CosmoHndl cosmo, SourceHndl source){
	cosmo->setOmega_matter(map->omegam,true);
	cosmo->sethubble(map->h);
	setZlens(cosmo,map->zlens,source->getZ());
	source->setZ(map->zsource);

	double fac = LH->DS/LH->DLS/LH->DL*LH->h;

	/// converts to the code units
	if(flag_MOKA_analyze == 0 || flag_MOKA_analyze == 2){
	  int i, j;
	  std::cout << "converting the units of the MOKA map" << std::endl;
	  for(i=0;i<map->nx;i++)
	    for(j=0;j<map->ny;j++){
	      int index = i*map->ny+j;
	      map->convergence[index] *= fac;
	      map->gamma1[index] *= fac;
	      map->gamma2[index] *= fac;
	    }
	}
}

/**
 * Sets many parameters within the MOKA lens model
 */

void MOKALens::assignParams(InputParams& params){


  if(!params.get("z_lens",zlens)){
	  ERROR_MESSAGE();
	  std::cout << "parameter z_lens needs to be set in parameter file " << params.filename() << std::endl;
	  exit(0);
  }

  if(!params.get("MOKA_input_file",MOKA_input_file)){
	  ERROR_MESSAGE();
	  std::cout << "parameter MOKA_input_file needs to be set in parameter file " << params.filename() << std::endl;
	  exit(0);
  }

  if(!params.get("background_field",flag_background_field)){
	  ERROR_MESSAGE();
	  std::cout << "parameter background_field needs to be set in parameter file " << params.filename() << std::endl;
	  exit(0);
  }

  if(!params.get("MOKA_analyze",flag_MOKA_analyze)) flag_MOKA_analyze = 0;

  set = true;
}

double MOKALens::getZlens(){
	return zlens;
}

void MOKALens::setZlens(CosmoHndl cosmo,double z,double dummy){
	zlens = z;
}

/**
 * saves the image, by reading off the calues from the image tree
 * and then saving to a fits file and computing the radial profile
 * of the convergence
 */
void MOKALens::saveImage(GridHndl grid,bool saveprofiles){
	std::stringstream f;
	std::string filename;

	if(flag_background_field==1) f << MOKA_input_file << "_only_noise.fits";
	else{
	  if(flag_MOKA_analyze == 0) f << MOKA_input_file << "_noisy.fits";
	  else f << MOKA_input_file << "_no_noise.fits";
	}
	filename = f.str();

	MoveToTopList(grid->i_tree->pointlist);

	do{
		long index = IndexFromPosition(grid->i_tree->pointlist->current->x,map->nx,map->boxlrad,map->center);
		if(index > -1){
			map->convergence[index] = grid->i_tree->pointlist->current->kappa;
			map->gamma1[index] = grid->i_tree->pointlist->current->gamma[0];
			map->gamma2[index] = grid->i_tree->pointlist->current->gamma[1];
			map->gamma3[index] = grid->i_tree->pointlist->current->gamma[2];
		}
	}while(MoveDownList(grid->i_tree->pointlist)==true);

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
	  double RE3,xxc,yyc;
          saveProfiles(RE3,xxc,yyc);
	  estSignLambdas();
	  double RE1,RE2;
	  EinsteinRadii(RE1,RE2,xxc,yyc);
	  std::ostringstream fEinr;
	  if(flag_background_field==1) fEinr << MOKA_input_file << "_only_noise_Einstein.radii.dat";
	  else{
	    if(flag_MOKA_analyze == 0) fEinr << MOKA_input_file << "_noisy_Einstein.radii.dat";
	    else fEinr << MOKA_input_file << "_no_noise_Einstein.radii.dat";
	  }
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
void MOKALens::saveProfiles(double &RE3,double &xxc,double &yyc){
	/* measuring the differential and cumulative profile*/
	double xmin = -map->boxlMpc*0.5;
	double xmax =  map->boxlMpc*0.5;
	double drpix = map->boxlMpc/map->nx;

	int galaxiesPerBin = 64;

	std::valarray<float> pxdist(map->nx*map->ny);
	std::valarray<float> red_sgE(map->nx*map->ny),red_sgB(map->nx*map->ny),sgm(map->nx*map->ny); 
	int i, j;
	/*
	  measure the center of mass
	  double xcm=0,ycm=0,tot=0;
	  for(i=0; i<map->nx; i++ ) for(j=0; j<map->ny; j++ ){
	  xcm+=map->convergence[i+map->ny*j]*(xmin+(drpix*0.5)+i*drpix);
	  ycm+=map->convergence[i+map->ny*j]*(xmin+(drpix*0.5)+j*drpix);
	  tot+=map->convergence[i+map->ny*j];
	  }	
	  xcm/=tot;
	  ycm/=tot;
	  xxc = xcm;
	  yyc = ycm;
	*/
	// moving center
        double xcm,ycm;
	if(flag_background_field==1){
	  xcm = 0.;
	  ycm = 0.;
	}
	else{
	  cmass(map->ny,map->convergence,map->x,xcm,ycm);
	}
	xxc = xcm;
	yyc = ycm;
	int ai = locate(map->x,xxc);
	ai = ((ai > 0) ? ai:0);
	ai = ((ai < map->nx-1) ? ai:map->nx-1);
	int bi = locate(map->x,yyc);
	bi = ((bi > 0) ? bi:0);
	bi = ((bi < map->nx-1) ? bi:map->nx-1);
	std:: cout << "  ------------ center ------------- " << std:: endl;
	std:: cout << "    " << xxc << "  " << yyc << std:: endl;
	std:: cout << "    " << ai << "  " << bi << std:: endl;
	std:: cout << "  ------------------r ------------- " << std:: endl;
	for(i=0; i<map->nx; i++ ) for(j=0; j<map->ny; j++ ){
	    pxdist[i+map->ny*j]= sqrt(pow((xmin+(drpix*0.5)+i*drpix-xcm),2) +
				pow((xmin+(drpix*0.5)+j*drpix-ycm),2));
		// reduced shear E a B
		double dx=map->x[i];
		double dy=map->x[j];
		double p=atan2( dy, dx ); // check gamma 1 and gamma 2 definition
		red_sgE[i+map->ny*j] = (-map->gamma1[i+map->ny*j]*cos(2*p)-map->gamma2[i+map->ny*j]*sin(2*p))/(1.-map->convergence[i+map->ny*j]);
		red_sgB[i+map->ny*j] = (map->gamma1[i+map->ny*j]*sin(2*p)-map->gamma2[i+map->ny*j]*cos(2*p))/(1.-map->convergence[i+map->ny*j]);
		sgm[i+map->ny*j] = sqrt(pow(map->gamma1[i+map->ny*j],2) + pow(map->gamma2[i+map->ny*j],2));
	}

	double dr0 = 8.*(0.5*map->boxlMpc)/(map->nx/2.);
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
    // number of galaxies per arcminsquare
	double dntbggal=double(nbggal)*(map->boxlarcsec*map->boxlarcsec)/3600.;
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

	  dr0 = galaxiesPerBin*(0.5*map->boxlMpc)/(0.5*ntbggal);

	  nbin = int(xmax/dr0);                                                        
	}
	//                                                                                           
	std:: cout << "   " << std:: endl;                                                           
	std:: cout << " nbins = " << nbin << "  dr0 = " << dr0 << std:: endl;                        
	std:: cout << " ______________________________________________________ " << std:: endl;      
	std:: cout << " computing profiles assuming spherical symmetry";                                
	// - - - - - - - - - - - - - - - - -

	// TODO Carlo:  These are all memory leaks!  They are never deleted!
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
	else{
	  if(flag_MOKA_analyze == 0) fprof << MOKA_input_file << "_noisy_MAP_radial_prof.dat";
	  else fprof << MOKA_input_file << "_no_noise_MAP_radial_prof.dat";
	}
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
	if((RE3-lrOFr[0])<1e-4 || (RE3-lrOFr[nbin-1])<1e-4) RE3=0;
	else RE3 = pow(10.,RE3)*map->inarcsec;
}

/** \ingroup DeflectionL2
   *
   * \brief Routine for obtaining the deflection and other lensing quantities for
   * a MOKA map (MOKALens), for just one ray!!
   *
*/
void MOKALens::rayshooterInternal(double *xx, double *alpha, KappaType *gamma, KappaType *kappa, bool kappa_off){
    
  long index = IndexFromPosition(xx,map->nx,map->boxlMpc/map->h,map->center);

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
void MOKALens::EinsteinRadii(double &RE1, double &RE2, double &xxc, double &yyc){
  double signV;
  //  std:: vector<double> xci1,yci1;
  std:: vector<double> xci2,yci2;
  // open file readable by ds9
  std::ostringstream fcrit;
  if(flag_background_field==1) fcrit << MOKA_input_file << "_only_noise_Criticals.reg";
  else{
    if(flag_MOKA_analyze == 0) fcrit << MOKA_input_file << "_noisy_Criticals.reg";
    else fcrit << MOKA_input_file << "_no_noise_Criticals.reg";
  }
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
  double pixDinL = map->boxlMpc/double(map->nx);
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
  int nearest0groupid = fof(pixDinL,xci,yci,groupid);
  if(nearest0groupid>0){
    std:: vector<double> xcpoints,ycpoints;
    double xercm=0;
    double yercm=0;
    for(int ii=0;ii<nc;ii++){
      if(groupid[ii] == nearest0groupid){
	xcpoints.push_back(xci[ii]);
	ycpoints.push_back(yci[ii]);
	xercm+=xci[ii];
	yercm+=yci[ii];
      }
    }
    nc = xcpoints.size();
    xercm=xercm/double(nc);
    yercm=yercm/double(nc);
    double distcentre=sqrt((xercm-xxc)*(xercm-xxc)+(yercm-yyc)*(yercm-yyc))*map->inarcsec;
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
    // if is not in the centre
    std:: cout << "distance " <<  distcentre << std:: endl;
    if(distcentre>1.5*RE2){
      RE1=-RE1;
      RE2=-RE2;
    }
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
	else{
	  if(flag_MOKA_analyze == 0) f << MOKA_input_file << "_noisy.fits";
	  else f << MOKA_input_file << "_no_noise.fits";
	}
	filename = f.str();

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
	            double RE3,xxc,yyc;
	            saveProfiles(RE3,xxc,yyc);
		    estSignLambdas(); 
		    double RE1,RE2;
		    EinsteinRadii(RE1,RE2,xxc,yyc);
		    std::ostringstream fEinr;
		    if(flag_background_field==1) fEinr << MOKA_input_file << "_only_noise_Einstein.radii.dat";
		    else{
		      if(flag_MOKA_analyze == 0) fEinr << MOKA_input_file << "_noisy_Einstein.radii.dat";
		      else fEinr << MOKA_input_file << "_no_noise_Einstein.radii.dat";
		    }
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

