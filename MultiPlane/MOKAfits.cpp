/*
 * fits.cpp
 *
 *  Created on: Jun 19, 2012
 *      Author: mpetkova
 */

#include "MOKAlens.h"
#include <fstream>

#ifdef ENABLE_FITS
#include <CCfits/CCfits>
using namespace CCfits;

#endif

#ifdef ENABLE_FFTW
#include "fftw3.h"
#endif

void LensHaloMOKA::getDims(){
#ifdef ENABLE_FITS
	try{
	  std::auto_ptr<FITS> ff(new FITS (MOKA_input_file, Read));
	  
	  PHDU *h0=&ff->pHDU();
	  
	  map->nx=h0->axis(0);
	  map->ny=h0->axis(1);
	  std:: cout << "nx           ny " << std:: endl;
	  std:: cout << map->nx << "   " << map->ny << std:: endl;
	}
	catch(FITS::CantOpen){
	  std::cout << "can not open " << MOKA_input_file << std::endl;
	  exit(1);
	}
#else
	std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
	exit(1);
#endif
}

/**
 * \brief reads in the fits file for the MOKA map and saves it in the structure map
 */
void LensHaloMOKA::readImage(){
#ifdef ENABLE_FITS
	std:: cout << " reading MOKA file: " << MOKA_input_file << std:: endl;
	std::auto_ptr<FITS> ff(new FITS (MOKA_input_file, Read));

	PHDU *h0=&ff->pHDU();

        h0->readAllKeys();
 
	// try to read MVIR, if it exists is a MOKA map
	bool moka;
  try
  {
    h0->readKey ("MVIR",map->m);
    moka=true;
  }
  catch(CCfits::HDU::NoSuchKeyword)
  {
    moka=false;
  }
  
	// principal HDU is read
	h0->read(map->convergence);
  
	if(moka)
  {
    // check if they exist if not it has to compute them
	  int nhdu = h0->axis(2);
	  
    std::cout << "nhdu = " << nhdu << std::endl ;
    
	  if(nhdu==1)
    {
	    std:: cout << "  preProcessing Map MOKA fits file has only KAPPA map" << std:: endl;
	    PreProcessFFTWMap();
	  }
	  else
    {
	    ExtHDU &h1=ff->extension(1);
	    h1.read(map->alpha1);
	    ExtHDU &h2=ff->extension(2);
	    h2.read(map->alpha2);
	    ExtHDU &h3=ff->extension(3);
	    h3.read(map->gamma1);
	    ExtHDU &h4=ff->extension(4);
	    h4.read(map->gamma2);
      ExtHDU &h5=ff->extension(5);
	    h5.read(map->phi);
      
	    std::cout << *h0 << h1 << h2 << h3  << h4 << h5 << std::endl;
	  }
    
	  /* these are always present in each fits file created by MOKA */
	  h0->readKey ("SIDEL",map->boxlarcsec);
	  h0->readKey ("SIDEL2",map->boxlMpc);  // recall you that MOKA Mpc/h
	  h0->readKey ("ZLENS",map->zlens);
	  h0->readKey ("ZSOURCE",map->zsource);
	  h0->readKey ("OMEGA",map->omegam);
	  h0->readKey ("LAMBDA",map->omegal);
	  h0->readKey ("H",map->h);
	  h0->readKey ("W",map->wq);
	  h0->readKey ("MSTAR",map->mstar);
	  // h0->readKey ("MVIR",map->m);
	  h0->readKey ("CONCENTRATION",map->c);
	  h0->readKey ("DL",map->Dlens);
	  h0->readKey ("DLS",map->DLS);
	  h0->readKey ("DS",map->DS);
    
	}
  else
  {
    
	  int npixels = map->nx;
	  if(map->ny<map->nx) npixels = map->ny;
	  std:: valarray<double> mapbut(map->nx*map->ny);
	  for(int i=0;i<map->nx;i++) for(int j=0;j<map->ny;j++)
    {
	      mapbut[i+map->nx*j] = map->convergence[i+map->nx*j];
    }
	  map->convergence.resize(npixels*npixels);
    
	  // FITS file must contain the following keywords:
	  /*
            ZLENS or REDSHIFT                                 double
            PHYSICALSIZE                                      double
	  */
	  
	  try
    {
            h0->readKey("ZLENS",map->zlens);
	  }
	  catch(CCfits::HDU::NoSuchKeyword)
    {
            try
            {
              h0->readKey("REDSHIFT",map->zlens);
            }
            catch(CCfits::HDU::NoSuchKeyword)
            {
              std::cout << "unable to read map zlmap" << std::endl;
              std::cout << "I will STOP here!!!" << std::endl;
              exit(1);
            }
	  }
	  map->Dlens = cosmo.angDist(0.,map->zlens);	  
	  double inarcsec  = 180./M_PI/map->Dlens*60.*60.;
	  double pixLMpc,pixelunit;
	  try
    {
	    // physical size in degrees
            h0->readKey("PHYSICALSIZE",map->boxlarcsec);
            map->boxlarcsec*=60.*60;
            pixLMpc = map->boxlarcsec/npixels/inarcsec;  
            map->boxlMpc = pixLMpc*npixels;
	    try {
	      h0->readKey("PIXELUNIT",pixelunit);
	      pixelunit/=pixLMpc/pixLMpc;	      
	    }
	    catch(CCfits::HDU::NoSuchKeyword)
      {
	      std::cout << " unable to read map pixelunit" << std::endl;
	      std::cout << " check this out in MOKAfits.cpp " << std::endl;
	      std:: cout << " I will STOP here !!! " << std:: endl; 
	      exit(1);
	    }
	  }
	  catch(CCfits::HDU::NoSuchKeyword)
    {
            std::cout << "unable to read map physical size and pixelunit" << std::endl;
            std::cout << "assuming is the PixelizMap file" << std::endl;
	    
            map->boxlarcsec = 4*60.*60.;    // the square is 4x4 by hand
            pixLMpc = map->boxlarcsec/npixels/inarcsec;  
            map->boxlMpc = pixLMpc*npixels;
            pixelunit = 1.e+10/cosmo.gethubble()/pixLMpc/pixLMpc; // by hand
	  }
	  
	  double avkappa = 0;
	  
	  // made square // need to be
	  for(int i=0;i<npixels;i++) for(int j=0;j<npixels;j++)
    {
	      map->convergence[i+npixels*j] = mapbut[i+map->nx*j]*pixelunit;
	      avkappa += map->convergence[i+npixels*j];
    }
	  avkappa /= (npixels*npixels);
	  
	  //std:: cout << "  <<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>> " << std:: endl;
	  //std:: cout << "     average kappa = " << avkappa << std:: endl;
	  //std:: cout << "  <<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>> " << std:: endl;
	  for(int i=0;i<npixels;i++) for(int j=0;j<npixels;j++)
    {
	      map->convergence[i+npixels*j] = (map->convergence[i+npixels*j] - avkappa); 
    }
	  
	  // kappa is not divided by the critical surface density
	  // they don't need to be preprocessed by fact
	  // create alpha and gamma arrays by FFT
	  map->nx = map->ny = npixels;

#ifdef ENABLE_FFTW
	  std:: cout << "  preProcessing Map " << std:: endl;
	  PreProcessFFTWMap();
	  
#else
	  std::cout << "Please enable the preprocessor flag ENABLE_FFTW !" << std::endl;
	  exit(1);
#endif
	}
	
	std:: cout << map->boxlMpc << "  " << map->boxlarcsec << std:: endl;
	
	for(size_t i=0;i<map->convergence.size();++i)
  {
	  assert(map->convergence[i] == map->convergence[i]);
	}
#else
	std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
	exit(1);
#endif
}


/**
 * \brief write the fits file of the new MOKA map from the structure map
 */
void LensHaloMOKA::writeImage(std::string filename){
#ifdef ENABLE_FITS
	long naxis=2;
	long naxes[2]={map->nx,map->ny};

	std::auto_ptr<FITS> fout(0);

	try{
		fout.reset(new FITS(filename,FLOAT_IMG,naxis,naxes));
	}
	catch(FITS::CantCreate){
		std::cout << "Unable to open fits file " << filename << std::endl;
		ERROR_MESSAGE();
		exit(1);
	}

	std::vector<long> naxex(2);
	naxex[0]=map->nx;
	naxex[1]=map->ny;

	PHDU *phout=&fout->pHDU();

	phout->write( 1,map->nx*map->ny,map->convergence );

	phout->addKey ("SIDEL",map->boxlarcsec,"arcsec");
	phout->addKey ("SIDEL2",map->boxlMpc,"Mpc/h");
	phout->addKey ("ZLENS",map->zlens,"lens redshift");
	phout->addKey ("ZSOURCE",map->zsource, "source redshift");
	phout->addKey ("OMEGA",map->omegam,"omega matter");
	phout->addKey ("LAMBDA",map->omegal,"omega lamda");
	phout->addKey ("H",map->h,"hubble/100");
	phout->addKey ("W",map->wq,"dark energy equation of state parameter");
	phout->addKey ("MSTAR",map->mstar,"stellar mass of the BCG in Msun/h");
	phout->addKey ("MVIR",map->m,"virial mass of the halo in Msun/h");
	phout->addKey ("CONCENTRATION",map->c,"NFW concentration");
	phout->addKey ("DL",map->Dlens,"Mpc/h");
	phout->addKey ("DLS",map->DLS,"Mpc/h");
	phout->addKey ("DS",map->DS,"Mpc/h");

	ExtHDU *eh1=fout->addImage("gamma1", FLOAT_IMG, naxex);
	eh1->write(1,map->nx*map->ny,map->gamma1);
	ExtHDU *eh2=fout->addImage("gamma2", FLOAT_IMG, naxex);
	eh2->write(1,map->nx*map->ny,map->gamma2);
	ExtHDU *eh3=fout->addImage("gamma3", FLOAT_IMG, naxex);
	eh3->write(1,map->nx*map->ny,map->gamma3);

  ExtHDU *eh4=fout->addImage("phi", FLOAT_IMG, naxex);
	eh4->write(1,map->nx*map->ny,map->phi);
  
	std::cout << *phout << std::endl;
#else
	std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
	exit(1);
#endif
}

/**
 * routine used by fof to link nearby grid cell points
 */
void make_friendship(int ii,int ji,int np,std:: vector<int> &friends, std:: vector<double> &pointdist){
  for(int jj=0;jj<np;jj++){
    if(friends[ji+np*jj]!=0){
      if(friends[ji+np*jj]<0){
    	  friends[ii+np*jj]=-(ii+1);
      }
      else{
    	  friends[ii+np*jj]=(ii+1);
      }
      friends[ji+np*jj]=0;
    }
  }  
  friends[ii+np*ji]=-(ii+1);
}

/*
 * given a set of grid points xci and yci and an interpixeled distance l return the id of the
 * fof group nearest to the centre
 */

int fof(double l,std:: vector<double> xci, std:: vector<double> yci, std:: vector<int> &groupid){
  int np = xci.size();
  std:: vector<int> friends(np*np);
  std:: vector<double> pointdist(np*np);
  for(int ii = 0;ii<np; ii++) for(int ji = 0;ji<np; ji++){
      pointdist[ii+np*ji] = sqrt( pow(xci[ii] - xci[ji],2) + pow(yci[ii] - yci[ji],2));
      groupid[ii] = 0;
      friends[ii+np*ji]=0;
    }
  for(int ii=0;ii<np;ii++) for(int ji = 0;ji<np; ji++){
      if(pointdist[ii+np*ji]<=1.5*l) friends[ii+np*ji] = ii+1;
    }
  for(int ii=0;ii<np;ii++){
    int r = 0;
    while(r==0){
      r=1;
      for(int ji=0;ji<np;ji++){
	if(friends[ii+np*ji]>0){
	  if(ii!=ji){
	    make_friendship(ii,ji,np,friends,pointdist);
	    r=0;
	  }
	}
      }
    }
  }
  for(int ii=0;ii<np;ii++){
    int p=0;
    for(int ji=0;ji<np;ji++){
      if(friends[ji+np*ii]!=0) p++;
      if(p==2){
	std:: cout << ji << "  " << ii << ":  "  << friends[ji+np*ii] << "  " << friends[ii+np*ji] << std:: endl;
	exit(1);
      }
    }
  }
  // count the particles in each group
  int kt = 0;
  int ng= 0;
  std:: vector<double> distcentre;
  std:: vector<int> idgroup;
  for(int ii=0;ii<np;ii++){
    int k = 0;
    double xcm=0;
    double ycm=0;
    for(int ji=0;ji<np;ji++){
      if(friends[ii+np*ji]!=0){
	k++;
	groupid[ji]=ii+1;
	xcm += xci[ji];
	ycm += yci[ji];
      }
    }
    if(k>4){
      ng++;
      xcm/=k;
      ycm/=k;
      distcentre.push_back(sqrt(xcm*xcm+ycm*ycm));
      idgroup.push_back(ii+1);
      // std:: cout << "  " << ii+1 << "  " << k << "   " << sqrt(xcm*xcm+ycm*ycm) << std:: endl;
    }
    kt = kt + k;
  }
  if(kt != np){
    std:: cout << " number of screaned particles : " << kt << std:: endl;
    std:: cout << " differes from the number of particles : " << np << std:: endl;
    std:: cout << " number of group found : " << ng << std:: endl;
    std:: cout << "     " << std:: endl;
    std:: cout << " I will STOP here!!! " << std:: endl;
    exit(1);
  }
  if(idgroup.size()>0){
    std:: vector<double>::iterator it = min_element(distcentre.begin(), distcentre.end());
    int minpos = idgroup[distance(distcentre.begin(), it)];
    // std:: cout << " nearest to the centre " << minpos << std:: endl;
    return minpos;
  }
  else return 0;
  /* Make a histogram of the data */
  /*
  std::vector< int > histogram(np,0);
  std::vector< int >::iterator it = groupid.begin();
  while(it != groupid.end()) histogram[*it++]++;
  // Print out the frequencies of the values in v
  std::copy(histogram.begin(),histogram.end(),std::ostream_iterator< int >(std::cout, " "));
  std::cout << std::endl;
  // Find the mode
  int mode = std::max_element(histogram.begin(),histogram.end()) - histogram.begin();
  return mode;
  */
}

/**
 * \brief pre-process surface mass density map computing deflection angles and shear in FFT
 */

void LensHaloMOKA::PreProcessFFTWMap(){
    #ifdef ENABLE_FFTW
    // initialize the quantities and does the zeropadding
    int zerosize = 4;      // zero padding region
    int npix_filter = 0;   // filter the map if you want on a given number of pixels

    int npixels = map->nx;
    double boxl = map->boxlMpc;
  
    // size of the new map
    int Nnpixels = int(zerosize*npixels);
    double Nboxl = boxl*zerosize;
    
    // filter the map in the fourier space to avoid spikes if you want
    double sigmag = 0;
    if(npix_filter>0) sigmag = (2.*M_PI/boxl)*npixels/npix_filter;
    
    std:: valarray<double> Nmap( Nnpixels*Nnpixels );
        
    // assume locate in a squared map and build up the new map with the zero padding region
    for( int i=0; i<Nnpixels; i++ ) for( int j=0; j<Nnpixels; j++ ){
            Nmap[i+Nnpixels*j]=0;
            if(i>=int(Nnpixels/2-npixels/2) && i<int(Nnpixels/2+npixels/2) &&
               j>=int(Nnpixels/2-npixels/2) && j<int(Nnpixels/2+npixels/2)){
                int ii = i-int(Nnpixels/2-npixels/2);
                int jj = j-int(Nnpixels/2-npixels/2);
                
                if(ii>=npixels || jj>=npixels){
                    std::cout << " 1 error mapping " << ii << "  " << jj << std::endl; 
                    exit(1); 
                } 
                if(ii<0 || jj<0){ 
                    std::cout << " 2 error mapping " << ii << "  " << jj << std::endl; 
                    exit(1); 
                } 
                Nmap[i+Nnpixels*j]=map->convergence[ii+npixels*jj];
            }
      }
    
    double *dNmap=new double[Nnpixels*Nnpixels];
    double *input=new double[Nnpixels*Nnpixels];
    fftw_complex *fNmap=new fftw_complex[Nnpixels*(Nnpixels/2+1)];
    fftw_complex *output=new fftw_complex[Nnpixels*(Nnpixels/2+1)];
    for(int k=0;k<Nnpixels*Nnpixels;k++) dNmap[k] = double(Nmap[k]);
    fftw_plan p;
    p=fftw_plan_dft_r2c_2d(Nnpixels,Nnpixels,input,output,FFTW_ESTIMATE);
    for (int i=0; i<Nnpixels*Nnpixels; i++) input[i] = dNmap[i];
    fftw_execute( p );
    for(int i=0; i<Nnpixels*(Nnpixels/2+1);i++){
        fNmap[i][0] = output[i][0];
        fNmap[i][1] = output[i][1];
    }
    delete[] input;
    delete[] output;
    delete[] dNmap;

    // fourier space
    // std:: cout << " allocating fourier space maps " << std:: endl;
    fftw_complex *fphi   = new fftw_complex[Nnpixels*(Nnpixels/2+1)];
    fftw_complex *falpha1= new fftw_complex[Nnpixels*(Nnpixels/2+1)];
    fftw_complex *falpha2= new fftw_complex[Nnpixels*(Nnpixels/2+1)];
    fftw_complex *fgamma1= new fftw_complex[Nnpixels*(Nnpixels/2+1)];
    fftw_complex *fgamma2= new fftw_complex[Nnpixels*(Nnpixels/2+1)];
    
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nnpixels; i++ ){
        // kx = i if i<n/2 else i-n
        double kx=(i<Nnpixels/2)?double(i):double(i-Nnpixels);
        kx=kx*2.*M_PI/Nboxl;
        for( int j=0; j<Nnpixels/2+1; j++ ){
            double ky=double(j);
            ky=ky*2.*M_PI/Nboxl;
            // rescale respect to the box size
            double k2 = (kx*kx + ky*ky);
            // smooth if you want
            if(npix_filter>0){
                fNmap[j+(Nnpixels/2+1)*i][0] = fNmap[j+(Nnpixels/2+1)*i][0]*exp(-k2/sigmag/sigmag/2.);
                fNmap[j+(Nnpixels/2+1)*i][1] = fNmap[j+(Nnpixels/2+1)*i][1]*exp(-k2/sigmag/sigmag/2.);
            }
            // fphi
            fphi[j+(Nnpixels/2+1)*i][0]= -2.*fNmap[j+(Nnpixels/2+1)*i][0]/k2;
            fphi[j+(Nnpixels/2+1)*i][1]= -2.*fNmap[j+(Nnpixels/2+1)*i][1]/k2;
            // null for k2 = 0 no divergence
            if(k2 == 0){
                fphi[j+(Nnpixels/2+1)*i][0] = 0.;
                fphi[j+(Nnpixels/2+1)*i][1] = 0.;
            }
            // gamma
            fgamma1[j+(Nnpixels/2+1)*i][0] = 0.5*(kx*kx-ky*ky)*fphi[j+(Nnpixels/2+1)*i][0];
            fgamma1[j+(Nnpixels/2+1)*i][1] = 0.5*(kx*kx-ky*ky)*fphi[j+(Nnpixels/2+1)*i][1];
            fgamma2[j+(Nnpixels/2+1)*i][0] = kx*ky*fphi[j+(Nnpixels/2+1)*i][0];
            fgamma2[j+(Nnpixels/2+1)*i][1] = kx*ky*fphi[j+(Nnpixels/2+1)*i][1];
            // alpha
            falpha1[j+(Nnpixels/2+1)*i][0] = -kx*fphi[j+(Nnpixels/2+1)*i][1];
            falpha1[j+(Nnpixels/2+1)*i][1] =  kx*fphi[j+(Nnpixels/2+1)*i][0];
            falpha2[j+(Nnpixels/2+1)*i][0] = -ky*fphi[j+(Nnpixels/2+1)*i][1];
            falpha2[j+(Nnpixels/2+1)*i][1] =  ky*fphi[j+(Nnpixels/2+1)*i][0];
        } 
    }
    
    fftw_destroy_plan(p);  
    delete[] fNmap;
    
    double *phi    = new double[Nnpixels*Nnpixels];
    double *alpha1 = new double[Nnpixels*Nnpixels];
    double *alpha2 = new double[Nnpixels*Nnpixels];
    double *gamma1 = new double[Nnpixels*Nnpixels];
    double *gamma2 = new double[Nnpixels*Nnpixels];

    fftw_plan pp;
    pp=fftw_plan_dft_c2r_2d(Nnpixels,Nnpixels,fphi,phi,FFTW_ESTIMATE);
    fftw_execute( pp );
    fftw_destroy_plan(pp);

    pp=fftw_plan_dft_c2r_2d(Nnpixels,Nnpixels,falpha1,alpha1,FFTW_ESTIMATE);
    fftw_execute( pp );
    fftw_destroy_plan(pp);

    pp=fftw_plan_dft_c2r_2d(Nnpixels,Nnpixels,falpha2,alpha2,FFTW_ESTIMATE);
    fftw_execute( pp );
    fftw_destroy_plan(pp);

    pp=fftw_plan_dft_c2r_2d(Nnpixels,Nnpixels,fgamma1,gamma1,FFTW_ESTIMATE);
    fftw_execute( pp );
    fftw_destroy_plan(pp);

    pp=fftw_plan_dft_c2r_2d(Nnpixels,Nnpixels,fgamma2,gamma2,FFTW_ESTIMATE);
    fftw_execute( pp );
    fftw_destroy_plan(pp);

    // std:: cout << " remapping the map in the original size " << std:: endl;

    map->phi.resize(npixels*npixels);
    map->gamma1.resize(npixels*npixels);
    map->gamma2.resize(npixels*npixels);
    map->alpha1.resize(npixels*npixels);
    map->alpha2.resize(npixels*npixels);

    for( int i=Nnpixels/2-npixels/2; i<Nnpixels/2+npixels/2; i++ ){
            for( int j=Nnpixels/2-npixels/2; j<Nnpixels/2+npixels/2; j++ ){
                int ii = i-int(Nnpixels/2-npixels/2);
                int jj = j-int(Nnpixels/2-npixels/2);

                map->phi[ii+npixels*jj] = float( phi[i+Nnpixels*j]/Nnpixels/Nnpixels);
              
                map->gamma1[ii+npixels*jj] = float( gamma1[i+Nnpixels*j]/Nnpixels/Nnpixels);
                map->gamma2[ii+npixels*jj] = float(-gamma2[i+Nnpixels*j]/Nnpixels/Nnpixels);
        
                map->alpha1[ii+npixels*jj] = float(alpha1[i+Nnpixels*j]/Nnpixels/Nnpixels);
                map->alpha2[ii+npixels*jj] = float(alpha2[i+Nnpixels*j]/Nnpixels/Nnpixels);
	    }
    }
    delete[] fphi;
    delete[] falpha1;
    delete[] falpha2;
    delete[] fgamma1;
    delete[] fgamma2;
    delete[] phi;
    delete[] alpha1;
    delete[] alpha2;
    delete[] gamma1;
    delete[] gamma2;
    #endif
}


