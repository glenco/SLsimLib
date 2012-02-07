/*
 * lens.c
 *
 *  Created on: Jan 12, 2012
 *      Author: mpetkova
 */

#include <slsimlib.h>

using namespace std;

Lens::Lens(char filename[]){
	readParamfile(filename);

	set = true;
}

Lens::Lens(){
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

	  cout << "reading from " << filename << endl;

	  if(!file_in){
	    cout << "Can't open file " << filename << endl;
	    exit(1);
	  }

	  // output file
	  file_in >> label >> outputfile;
	  cout << label << " " << outputfile << endl << endl;

	  // Nplanes file
	  file_in >> label >> Nplanes;
	  cout << label << " " << Nplanes << endl << endl;

	   // redshifts
	   file_in >> label >> zlens;
	   cout << label << " " << zlens << endl;
}

Source::Source(){
}

SourceUniform::SourceUniform(char *filename) : Source(){
	readParamfile(filename);
}

SourceGaussian::SourceGaussian(char *filename) : Source(){
	readParamfile(filename);
}

SourceBLR::SourceBLR(char *filename) : Source(){
	readParamfile(filename);
}

SourceBLRDisk::SourceBLRDisk(char *filename) : SourceBLR(filename){

}

SourceBLRSph1::SourceBLRSph1(char *filename) : SourceBLR(filename){

}

SourceBLRSph2::SourceBLRSph2(char *filename) : SourceBLR(filename){

}

Source::~Source(){
}

SourceUniform::~SourceUniform(){
}

SourceGaussian::~SourceGaussian(){
}

SourceBLR::~SourceBLR(){
}

SourceBLRDisk::~SourceBLRDisk(){
}

SourceBLRSph1::~SourceBLRSph1(){
}

SourceBLRSph2::~SourceBLRSph2(){
}

void SourceUniform::readParamfile(char *filename){
	  ifstream file_in(filename);
	  char label[20];
	  int type;

	  cout << "reading from " << filename << endl;

	  if(!file_in){
	    cout << "Can't open file " << filename << endl;
	    exit(1);
	  }

	  // source information
	   cout << "**Source structure**" << endl;

	   file_in >> label >> type;
	   source_sb_type = (SBModel)type;
	   cout << label << " " << source_sb_type << endl;

	   // redshifts
	   file_in >> label >> zsource;
	   cout << label << " " << zsource << endl;
}

void SourceGaussian::readParamfile(char *filename){
	  ifstream file_in(filename);
	  char label[20];
	  int type;

	  cout << "reading from " << filename << endl;

	  if(!file_in){
	    cout << "Can't open file " << filename << endl;
	    exit(1);
	  }

	  // source information
	   cout << "**Source structure**" << endl;

	   file_in >> label >> type;
	   source_sb_type = (SBModel)type;
	   cout << label << " " << source_sb_type << endl;

	   file_in >> label >> source_gauss_r2;
	   cout << label << " " << source_gauss_r2 << " Mpc" << endl;

	   // redshifts
	   file_in >> label >> zsource;
	   cout << label << " " << zsource << endl;
}

void SourceBLR::readParamfile(char *filename){
	  ifstream file_in(filename);
	  char label[20];
	  int type;

	  cout << "reading from " << filename << endl;

	  if(!file_in){
	    cout << "Can't open file " << filename << endl;
	    exit(1);
	  }

	  // source information
	   cout << "**Source structure**" << endl;

	   file_in >> label >> type;
	   source_sb_type = (SBModel)type;
	   cout << label << " " << source_sb_type << endl;

	   cout << "BLR surface brightness source" << endl;

	   switch(source_sb_type){
	   case BLR_Disk:
		   cout << "disk model" << endl;
		   break;
	   case BLR_Sph1:
		   cout << "spherical with circular orbits" << endl;
		   break;
	   case BLR_Sph2:
		   cout << "spherical with Gaussian velocities" << endl;
		   break;
	   default:
		   ERROR_MESSAGE();
		   cout << "ERROR: no submass internal profile chosen" << endl;
		   exit(1);
		   break;
	   }

	   file_in >> label >> source_BHmass;
	   cout << label << " " << source_BHmass << " Msun" << endl;

	   file_in >> label >> source_gamma;
	   cout << label << " " << source_gamma << endl;

	   file_in >> label >> source_inclination;
	   cout << label << " " << source_inclination << " deg" << endl;

	   file_in >> label >> source_opening_angle;
	   cout << label << " " << source_opening_angle << " deg" << endl;

	   source_inclination *= pi/180;
	   source_opening_angle *= pi/180;

	   file_in >> label >> source_r_in;
	   cout << label << " " << source_r_in << " Mpc" << endl;

	   file_in >> label >> source_r_out;
	   cout << label << " " << source_r_out << " Mpc" << endl;

	   file_in >> label >> source_nuo;
	   cout << label << " " << source_nuo << " Hz" << endl;

	   file_in >> label >> source_fK;
	   cout << label << " " << source_fK << " X V_Kepler" << endl;

	   source_monocrome = false;  // default value

	   // redshifts
	   file_in >> label >> zsource;
	   cout << label << " " << zsource << endl;
}

double SourceUniform::source_sb_func(double *y){
	return (double)( (y[0]*y[0] + y[1]*y[1]) < source_r*source_r );
}

double SourceGaussian::source_sb_func(double *y){
	return exp( -(y[0]*y[0] + y[1]*y[1])/source_gauss_r2 );
}

// surface brightness for models of the Broad Line Region
double SourceBLRDisk::source_sb_func(double *y){
	return blr_surface_brightness_disk(y,this);
}

double SourceBLRSph1::source_sb_func(double *y){
	return blr_surface_brightness_spherical_circular_motions(sqrt(y[0]*y[0] + y[1]*y[1]),this);
}
double SourceBLRSph2::source_sb_func(double *y){
	return blr_surface_brightness_spherical_random_motions(sqrt(y[0]*y[0] + y[1]*y[1]),this);
}

void in_source(double *y_source,ListHndl sourcelist){
  return;
}
