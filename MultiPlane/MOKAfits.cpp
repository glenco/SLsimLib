/*
 * fits.cpp
 *
 *  Created on: Jun 19, 2012
 *      Author: mpetkova
 */

#ifdef WITH_MOKA

#include <MOKAfits.h>
#include <CCfits/CCfits>

using namespace CCfits;

void getDims(std::string fn
	     ,int *nx
	     ,int *ny){

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
}

/*
 * reads in the fits file for the MOKA map and saves it in the structure map
 */
void readImage(std::string fn
		,std::valarray<float> *convergence
		,std::valarray<float> *alpha1
		,std::valarray<float> *alpha2
		,std::valarray<float> *gamma1
		,std::valarray<float> *gamma2
	       ,struct LensHalo *LH){ 
  /*
                ,double *boxl
		,double *boxlMpc
		,double *zlens
		,double *zsource
		,double *omegam
		,double *omegal
		,double *h
		,double *DL){
  */

	int nx,ny;

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
	/*
	h0->readKey ("SIDEL",*boxl);
	h0->readKey ("SIDEL2",*boxlMpc);
	h0->readKey ("ZLENS",*zlens);
	h0->readKey ("ZSOURCE",*zsource);
	h0->readKey ("OMEGA",*omegam);
	h0->readKey ("LAMBDA",*omegal);
	h0->readKey ("H",*h);
	*/
	ExtHDU &h1=ff->extension(1);
	h1.read(*alpha1);
	ExtHDU &h2=ff->extension(2);
	h2.read(*alpha2);
	ExtHDU &h3=ff->extension(3);
	h3.read(*gamma1);
	ExtHDU &h4=ff->extension(4);
	h4.read(*gamma2);

	std::cout << *h0 << h1 << h2 << h3  << h4 << std::endl;
}

/*
 * write the fits file of the new MOKA map from the structure map
 */
void writeImage(std::string filename
		,std::valarray<float> convergence
		,std::valarray<float> gamma1
		,std::valarray<float> gamma2
		,std::valarray<float> gamma3
		,int nx
		,int ny
		,struct LensHalo LH){ 
  /*
		,double boxl
		,double boxlMpc
		,double zlens
		,double zsource
		,double omegam
		,double omegal
		,double h
		,double DL){
  */

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

	phout->addKey ("SIDEL",LH.boxlarcsec,"arcsec");
	phout->addKey ("SIDEL2",LH.boxlMpc,"Mpc/h");
	phout->addKey ("ZLENS",LH.zl,"lens redshift");
	phout->addKey ("ZSOURCE",LH.zs, "source redshift");
	phout->addKey ("OMEGA",LH.omegam,"omega matter");
	phout->addKey ("LAMBDA",LH.omegal,"omega lamda");
	phout->addKey ("H",LH.h,"hubble/100");
	phout->addKey ("W",LH.wq,"dark energy equation of state parameter");
	phout->addKey ("MSTAR",LH.mstar,"stellar mass of the BCG in Msun/h");  
	phout->addKey ("MVIR",LH.m,"virial mass of the halo in Msun/h");  
	phout->addKey ("CONCENTRATION",LH.c,"NFW concentration");
	phout->addKey ("DL",LH.DL,"Mpc/h");
	phout->addKey ("DLS",LH.DLS,"Mpc/h");
	phout->addKey ("DS",LH.DS,"Mpc/h");
	

	ExtHDU *eh1=fout->addImage("gamma1", FLOAT_IMG, naxex);
	eh1->write(1,nx*ny,gamma1);
	ExtHDU *eh2=fout->addImage("gamma2", FLOAT_IMG, naxex);
	eh2->write(1,nx*ny,gamma2);
	ExtHDU *eh3=fout->addImage("gamma3", FLOAT_IMG, naxex);
	eh3->write(1,nx*ny,gamma3);

	std::cout << *phout << std::endl;

}

#endif
