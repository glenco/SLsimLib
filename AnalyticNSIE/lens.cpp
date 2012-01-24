/*
 * lens.c
 *
 *  Created on: Jan 12, 2012
 *      Author: mpetkova
 */

#include <slsimlib.h>

Lens::Lens(char filename[]){
	readParamfile(filename);

	set = true;
}

Lens::~Lens(){
}

int Lens::getNplanes(){
	return Nplanes;
}

void Lens::readParamfile(char *filename){
	  ifstream file_in(filename);
	  char label[20];
	  int i, type;
	  double tmp = 0;

	  std::cout << "reading from " << filename << std::endl;

	  if(!file_in){
	    std::cout << "Can't open file " << filename << std::endl;
	    exit(1);
	  }

	  // output file
	  file_in >> label >> outputfile;
	  std::cout << label << " " << outputfile << std::endl << std::endl;

	  // Nplanes file
	  file_in >> label >> Nplanes;
	  std::cout << label << " " << Nplanes << std::endl << std::endl;

	  // source information
	   std::cout << "**Source structure**" << std::endl;

	   file_in >> label >> type;
	   source_sb_type = (SBModel)type;
	   std::cout << label << " " << source_sb_type << std::endl;

	   if(source_sb_type == Uniform){
	 		  source_sb_func = uniform_SB;
	 		  std::cout << "uniform surface brightness source" << std::endl;
	   }else if(source_sb_type == Gaussian){
	 		  source_sb_func = gaussian_SB;
	 		  std::cout << "Gaussian surface brightness source" << std::endl;
	 		  file_in >> label >> source_gauss_r2;
	 		  std::cout << label << source_gauss_r2 << " Mpc" << std::endl;
	   }else{

	 	  std::cout << "BLR surface brightness source" << std::endl;
	 	  switch(source_sb_type){
	 	  case BLR_Disk:
	 		  source_sb_func = BLR_Disk_SB;
	 		  std::cout << "disk model" << std::endl;
	 		  break;
	 	  case BLR_Sph1:
	 		  source_sb_func = BLR_Sph1_SB;
	 		  std::cout << "spherical with circular orbits" << std::endl;
	 		  break;
	 	  case BLR_Sph2:
	 		  source_sb_func = BLR_Sph2_SB;
	 		  std::cout << "spherical with Gaussian velocities" << std::endl;
	 		  break;
	 	  default:
	 		  ERROR_MESSAGE();
	 		  std::cout << "ERROR: no submass internal profile chosen" << std::endl;
	 		  exit(1);
	 		  break;
	 	  }

	 	  file_in >> label >> source_BHmass;
	 	  std::cout << label << source_BHmass << " Msun" << std::endl;

	 	  file_in >> label >> source_gamma;
	 	  std::cout << label << source_gamma << std::endl;

	 	  file_in >> label >> source_inclination;
	 	  std::cout << label << source_inclination << " deg" << std::endl;

	 	  file_in >> label >> source_opening_angle;
	 	  std::cout << label << source_opening_angle << " deg" << std::endl;


	 	  source_inclination *= pi/180;
	 	  source_opening_angle *= pi/180;

	 	  file_in >> label >> source_r_in;
	 	  std::cout << label << source_r_in << " Mpc" << std::endl;

	 	  file_in >> label >> source_r_out;
	 	  std::cout << label << source_r_out << " Mpc" << std::endl;

	 	  file_in >> label >> source_nuo;
	 	  std::cout << label << source_nuo << " Hz" << std::endl;

	 	  file_in >> label >> source_fK;
	 	  std::cout << label << source_fK << " X V_Kepler" << std::endl;

	 	  source_monocrome = false;  // default value

	   }

	   // redshifts

	   file_in >> label >> zlens;
	   std::cout << label << " " << zlens << std::endl;

	   file_in >> label >> zsource;
	   std::cout << label << " " << zsource << std::endl;
}

void Lens::setInternal(CosmoHndl cosmo){
  MpcToAsec = 60*60*180 / pi / cosmo->angDist(0,zlens);
  std::cout << "Arcseconds/Mpc: " << MpcToAsec << std::endl;

  // find critical density
  Sigma_crit=cosmo->angDist(0,zsource)/cosmo->angDist(zlens,zsource)
		  /cosmo->angDist(0,zlens)/4/pi/Grav;

  to = (1+zlens)*cosmo->angDist(0,zsource)
		  /cosmo->angDist(zlens,zsource)/cosmo->angDist(0,zlens)
		  /8.39428142e-10;

  std::cout << "critical density is " << Sigma_crit << " Msun/Mpc^2   D_l = "
       << cosmo->angDist(0,zlens) << " Mpc D_s = " << cosmo->angDist(0,zsource) << " Mpc  to = " << to << " days/Mpc^2" << std::endl;

}


