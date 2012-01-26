/*
 * model.cpp
 *
 *  Created on: Jan 23, 2012
 *      Author: mpetkova
 */

#include <slsimlib.h>

Model::Model(char *paramfile, char *file_source, int flag){
	cosmo = new COSMOLOGY;
	source = new Source(file_source);

	switch(flag){
	case 0:
		lens = new AnaLens(paramfile);
		break;

	case 1:
		lens = new multiLens(paramfile);
		break;

	default:
		break;
	}

	setInternal();
}

Model::~Model(){
	delete lens;
	delete source;
	delete cosmo;
}

double Model::uniform_SB(double *y){
	return (double)( (y[0]*y[0] + y[1]*y[1]) < source->source_r*source->source_r );
}

double Model::gaussian_SB(double *y){
	return exp( -(y[0]*y[0] + y[1]*y[1])/source->source_gauss_r2 );
}

// surface brightness for models of the Broad Line Region
double Model::BLR_Disk_SB(double *y){
	return blr_surface_brightness_disk(y,this);
}

double Model::BLR_Sph1_SB(double *y){
	return blr_surface_brightness_spherical_circular_motions(sqrt(y[0]*y[0] + y[1]*y[1]),this);
}
double Model::BLR_Sph2_SB(double *y){
	return blr_surface_brightness_spherical_random_motions(sqrt(y[0]*y[0] + y[1]*y[1]),this);
}

void in_source(double *y_source,ListHndl sourcelist){
  return;
}

void Model::setInternal(){
	  lens->MpcToAsec = 60*60*180 / pi / cosmo->angDist(0,lens->zlens);

	  // in Mpc
	  lens->host_ro=4*pi*pow(lens->host_sigma/2.99792e5,2)*cosmo->angDist(0,lens->zlens)
			  *cosmo->angDist(lens->zlens,source->zsource)
			  /cosmo->angDist(0,source->zsource)/(1+lens->zlens);

	  // find critical density
	  lens->Sigma_crit=cosmo->angDist(0,source->zsource)/cosmo->angDist(lens->zlens,source->zsource)
			  /cosmo->angDist(0,lens->zlens)/4/pi/Grav;

	  lens->to = (1+lens->zlens)*cosmo->angDist(0,source->zsource)
			  /cosmo->angDist(lens->zlens,source->zsource)/cosmo->angDist(0,lens->zlens)
			  /8.39428142e-10;

	   if(source->source_sb_type == source->Uniform){
		   source_sb_func = &Model::uniform_SB;
	   }else if(source->source_sb_type == source->Gaussian){
		   source_sb_func = &Model::gaussian_SB;
	   }else if(source->source_sb_type ==  source->BLR_Disk){
		   source_sb_func = &Model::BLR_Disk_SB;
	   }else if(source->source_sb_type ==  source->BLR_Sph1){
		   source_sb_func = &Model::BLR_Sph1_SB;
	   }else if(source->source_sb_type == source->BLR_Sph2){
		   source_sb_func = &Model::BLR_Sph2_SB;
	   }
}

