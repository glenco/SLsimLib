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

void LensHaloMOKA::getDims(){
#ifdef ENABLE_FITS
	try{
		std::auto_ptr<FITS> ff(new FITS (MOKA_input_file, Read));

		PHDU *h0=&ff->pHDU();

		map->nx=h0->axis(0);
		map->ny=h0->axis(1);
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

	h0->read(map->convergence);

	h0->readKey ("SIDEL",map->boxlarcsec);
	h0->readKey ("SIDEL2",map->boxlMpc);
	h0->readKey ("ZLENS",map->zlens);
	h0->readKey ("ZSOURCE",map->zsource);
	h0->readKey ("OMEGA",map->omegam);
	h0->readKey ("LAMBDA",map->omegal);
	h0->readKey ("H",map->h);
	h0->readKey ("W",map->wq);
	h0->readKey ("MSTAR",map->mstar);
	h0->readKey ("MVIR",map->m);
	h0->readKey ("CONCENTRATION",map->c);
	h0->readKey ("DL",map->DL);
	h0->readKey ("DLS",map->DLS);
	h0->readKey ("DS",map->DS);

	ExtHDU &h1=ff->extension(1);
	h1.read(map->alpha1);
	ExtHDU &h2=ff->extension(2);
	h2.read(map->alpha2);
	ExtHDU &h3=ff->extension(3);
	h3.read(map->gamma1);
	ExtHDU &h4=ff->extension(4);
	h4.read(map->gamma2);

	std::cout << *h0 << h1 << h2 << h3  << h4 << std::endl;
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
	phout->addKey ("DL",map->DL,"Mpc/h");
	phout->addKey ("DLS",map->DLS,"Mpc/h");
	phout->addKey ("DS",map->DS,"Mpc/h");


	ExtHDU *eh1=fout->addImage("gamma1", FLOAT_IMG, naxex);
	eh1->write(1,map->nx*map->ny,map->gamma1);
	ExtHDU *eh2=fout->addImage("gamma2", FLOAT_IMG, naxex);
	eh2->write(1,map->nx*map->ny,map->gamma2);
	ExtHDU *eh3=fout->addImage("gamma3", FLOAT_IMG, naxex);
	eh3->write(1,map->nx*map->ny,map->gamma3);

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
  /* Print out the frequencies of the values in v */
  // std::copy(histogram.begin(),histogram.end(),std::ostream_iterator< int >(std::cout, " "));
  // std::cout << std::endl;
  /* Find the mode */
  /*
  int mode = std::max_element(histogram.begin(),histogram.end()) - histogram.begin();
  return mode;
  */
}

