/*
 * fits.cpp
 *
 *  Created on: Jun 19, 2012
 *      Author: mpetkova
 */

#include <MOKAfits.h>
#include <fstream>

#ifdef ENABLE_FITS
#include <CCfits/CCfits>

using namespace CCfits;

#endif

void getDims(std::string fn
	     ,int *nx
	     ,int *ny){
#ifdef ENABLE_FITS
	try{
		std::auto_ptr<FITS> ff(new FITS (fn, Read));

		PHDU *h0=&ff->pHDU();

		*nx=h0->axis(0);
		*ny=h0->axis(1);
	}
	catch(FITS::CantOpen){
		std::cout << "can not open " << fn << std::endl;
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
void readImage(std::string fn
		,std::valarray<float> *convergence
		,std::valarray<float> *alpha1
		,std::valarray<float> *alpha2
		,std::valarray<float> *gamma1
		,std::valarray<float> *gamma2
	    ,LensHalo *LH){
#ifdef ENABLE_FITS
	int nx,ny;

	std:: cout << " reading MOKA file: " << fn << std:: endl;
	std::ostringstream checkfout;
	checkfout << fn << "_noisy.fits";	
	std:: string checkfilenameout = checkfout.str();
	std:: ifstream checkfileout;
	checkfileout.open(checkfilenameout.c_str());
	if(checkfileout.is_open()){
	  std:: cout << "  " << std:: endl;
	  // std:: cout << checkfilenameout << " exists I will STOP here " << std:: endl;
	  std:: cout << "     halo already processed! " << std:: endl;
	  std:: cout << "   I will measure the map proprieties only " << std:: endl;
	  // exit(1);
	  fn = fn + "_noisy.fits";	
	  std:: cout << " I am reading " << fn << std:: endl;
	  std:: cout << "  " << std:: endl;
	}

	std::auto_ptr<FITS> ff(new FITS (fn, Read));

	PHDU *h0=&ff->pHDU();

	nx=h0->axis(0);
	ny=h0->axis(1);

	h0->read(*convergence);

	h0->readKey ("SIDEL",LH->boxlarcsec);
	h0->readKey ("SIDEL2",LH->boxlMpc);
	h0->readKey ("ZLENS",LH->zl);
	h0->readKey ("ZSOURCE",LH->zs);
	h0->readKey ("OMEGA",LH->omegam);
	h0->readKey ("LAMBDA",LH->omegal);
	h0->readKey ("H",LH->h);
	h0->readKey ("W",LH->wq);
	h0->readKey ("MSTAR",LH->mstar);  
	h0->readKey ("MVIR",LH->m);  
	h0->readKey ("CONCENTRATION",LH->c);
	h0->readKey ("DL",LH->DL);
	h0->readKey ("DLS",LH->DLS);
	h0->readKey ("DS",LH->DS);

	ExtHDU &h1=ff->extension(1);
	h1.read(*alpha1);
	ExtHDU &h2=ff->extension(2);
	h2.read(*alpha2);
	ExtHDU &h3=ff->extension(3);
	h3.read(*gamma1);
	ExtHDU &h4=ff->extension(4);
	h4.read(*gamma2);

	std::cout << *h0 << h1 << h2 << h3  << h4 << std::endl;
#else
	std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
	exit(1);
#endif
}


/**
 * \brief write the fits file of the new MOKA map from the structure map
 */
void writeImage(std::string filename
		,std::valarray<float> convergence
		,std::valarray<float> gamma1
		,std::valarray<float> gamma2
		,std::valarray<float> gamma3
		,int nx
		,int ny
		,LensHalo *LH){
#ifdef ENABLE_FITS
	long naxis=2;
	long naxes[2]={nx,ny};

	std::auto_ptr<FITS> fout(0);

	try{
		fout.reset(new FITS(filename,FLOAT_IMG,naxis,naxes));
	}
	catch(FITS::CantCreate){
		exit(1);
	}

	std::vector<long> naxex(2);
	naxex[0]=nx;
	naxex[1]=ny;

	PHDU *phout=&fout->pHDU();

	phout->write( 1,nx*ny,convergence );

	phout->addKey ("SIDEL",LH->boxlarcsec,"arcsec");
	phout->addKey ("SIDEL2",LH->boxlMpc,"Mpc/h");
	phout->addKey ("ZLENS",LH->zl,"lens redshift");
	phout->addKey ("ZSOURCE",LH->zs, "source redshift");
	phout->addKey ("OMEGA",LH->omegam,"omega matter");
	phout->addKey ("LAMBDA",LH->omegal,"omega lamda");
	phout->addKey ("H",LH->h,"hubble/100");
	phout->addKey ("W",LH->wq,"dark energy equation of state parameter");
	phout->addKey ("MSTAR",LH->mstar,"stellar mass of the BCG in Msun/h");
	phout->addKey ("MVIR",LH->m,"virial mass of the halo in Msun/h");
	phout->addKey ("CONCENTRATION",LH->c,"NFW concentration");
	phout->addKey ("DL",LH->DL,"Mpc/h");
	phout->addKey ("DLS",LH->DLS,"Mpc/h");
	phout->addKey ("DS",LH->DS,"Mpc/h");


	ExtHDU *eh1=fout->addImage("gamma1", FLOAT_IMG, naxex);
	eh1->write(1,nx*ny,gamma1);
	ExtHDU *eh2=fout->addImage("gamma2", FLOAT_IMG, naxex);
	eh2->write(1,nx*ny,gamma2);
	ExtHDU *eh3=fout->addImage("gamma3", FLOAT_IMG, naxex);
	eh3->write(1,nx*ny,gamma3);

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
 * given a a set of grid points xci and yci and an interpixeld distance l return the id of the 
 * fof group of each cell point
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
      // consider as friends grid points less distant than 1.5 x l
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
  for(int ii=0;ii<np;ii++){
    int k = 0;
    for(int ji=0;ji<np;ji++){
      if(friends[ii+np*ji]!=0){
	k++;
	groupid[ji]=ii+1;
      }
    }
    if(k>0){
      ng++;
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
  /* Make a histogram of the data */
  std::vector< int > histogram(np,0);
  std::vector< int >::iterator it = groupid.begin();
  while(it != groupid.end()) histogram[*it++]++;
  int mode = std::max_element(histogram.begin(),histogram.end()) - histogram.begin();
  return mode;
}
