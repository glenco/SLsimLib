/*
 * source.cpp
 *
 *  Created on: Feb 8, 2012
 *      Author: mpetkova
 */

#include <slsimlib.h>
#include <sstream>

using namespace std;

const string escape = "#";
char dummy[300];

Source::Source(){
}

SourceUniform::SourceUniform(string filename) : Source(){
	readParamfile(filename);
}

SourceGaussian::SourceGaussian(string filename) : Source(){
	readParamfile(filename);
}

SourceBLR::SourceBLR(string filename) : Source(){
	readParamfile(filename);
}

SourceBLRDisk::SourceBLRDisk(string filename) : SourceBLR(filename){

}

SourceBLRSph1::SourceBLRSph1(string filename) : SourceBLR(filename){

}

SourceBLRSph2::SourceBLRSph2(string filename) : SourceBLR(filename){

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

void SourceUniform::readParamfile(string filename){
	  string label, rlabel, rvalue;
	  stringstream ss;
	  double mydouble;
	  int flag;

	  label = "z_source";

	  cout << "uniform source: reading from " << filename << endl;

	  ifstream file_in(filename.c_str());
	  if(!file_in){
	    cout << "Can't open file " << filename << endl;
	    exit(1);
	  }

	  flag = 1;

	  // output file
	  while(!file_in.eof()){
		  file_in >> rlabel >> rvalue;
		  file_in.getline(dummy,100);

		  if(rlabel[0] == escape[0])
			  continue;

		  if(rlabel == label){

			  flag = 0;
			  ss << rvalue;

			  ss >> mydouble;
			  zsource = mydouble;

			  ss.clear();
			  ss.str(string());
		  }
	  }

	  file_in.close();


	  if(flag > 0){
		  ERROR_MESSAGE();
		  cout << "parameter " << label << " needs to be set!" << endl;
		  exit(0);
	  }

	  printSource();
}

void SourceGaussian::readParamfile(string filename){
	  string label[2], rlabel, rvalue;
	  stringstream ss;
	  void *addr[2];
	  int id[2];
	  int i;
	  double mydouble;
	  int flag;

	  addr[0] = &zsource;
	  id[0] = 0;
	  label[0] = "z_source";

	  addr[1] = &source_gauss_r2;
	  id[1] = 0;
	  label[1] = "gauss_r2";

	  cout << "gaussian source: reading from " << filename << endl;

	  ifstream file_in(filename.c_str());
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

		  flag = 1;

		  for(i = 0; i < 2; i++){

			  if(rlabel == label[i]){

				  flag = 0;
				  ss << rvalue;
				  ss >> mydouble;

				  *((double *)addr[i]) = mydouble;

				  ss.clear();
				  ss.str(string());
			  }

			  id[i] = -1;
		  }
	  }

	  file_in.close();

	  for(i = 0; i < 2; i++){
		  if(id[i] > 0){
			  ERROR_MESSAGE();
			  cout << "parameter " << label[i] << " needs to be set!" << endl;
			  exit(0);
		  }
	  }

	  printSource();
}

void SourceBLR::readParamfile(string filename){
	  string label[9], rlabel, rvalue;
	  stringstream ss;
	  void *addr[9];
	  int id[9];
	  int i, n;
	  float myfloat;
	  double mydouble;

	  n = 0;

	  addr[n] = &zsource;
	  id[n] = 1;
	  label[n++] = "z_source";

	  addr[n] = &source_BHmass;
	  id[n] = 0;
	  label[n++] = "BHmass";

	  addr[n] = &source_gamma;
	  id[n] = 0;
	  label[n++] = "gamma";

	  addr[n] = &source_inclination;
	  id[n] = 0;
	  label[n++] = "inclin";

	  addr[n] = &source_opening_angle;
	  id[n] = 0;
	  label[n++] = "opening_ang";

	  addr[n] = &source_r_in;
	  id[n] = 0;
	  label[n++] = "r_in";

	  addr[n] = &source_r_out;
	  id[n] = 0;
	  label[n++] = "r_out";

	  addr[n] = &source_nuo;
	  id[n] = 0;
	  label[n++] = "nuo";

	  addr[n] = &source_fK;
	  id[n] = 0;
	  label[n++] = "source_fk";

	  cout << "BLR source: reading from " << filename << endl;

	  ifstream file_in(filename.c_str());
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

		  for(i = 0; i < n; i++){
			  if(rlabel == label[i]){

				  ss << rvalue;

				  switch(id[i]){
				  case 0:
					  ss >> myfloat;
					  *((float *)addr[i]) = myfloat;
					  break;
				  case 1:
					  ss >> mydouble;
					  *((double *)addr[i]) = mydouble;
					  break;
				  }

				  ss.clear();
				  ss.str(string());

				  id[i] = -1;
			  }
		  }
	  }

	  file_in.close();

	  for(i = 0; i < n; i++){
		  if(id[i] > 0){
			  ERROR_MESSAGE();
			  cout << "parameter " << label[i] << " needs to be set!" << endl;
			  exit(0);
		  }
	  }

	  source_inclination *= pi/180;
	  source_opening_angle *= pi/180;
	  source_monocrome = false;

	  printSource();
}

void SourceUniform::printSource(){
	cout << endl << "**Source model**" << endl;

	cout << "z_source " << zsource << endl << endl;
}

void SourceGaussian::printSource(){
	cout << endl << "**Source model**" << endl;

	cout << "z_source " << zsource << endl;
	cout << "gauss_r2 " << source_gauss_r2 << endl << endl;
}

void SourceBLR::printSource(){
	cout << endl << "**Source model**" << endl;

	cout << "z_source " << zsource << endl;
	cout << "BHmass " << source_BHmass << endl;
	cout << "gamma " << source_gamma << endl;
	cout << "incl " << source_inclination*180/pi << endl;
	cout << "opening angl " << source_opening_angle*180/pi << endl;
	cout << "r_in " << source_r_in << endl;
	cout << "r_out " << source_r_out << endl;
	cout << "nuo " << source_nuo << endl;
	cout << "source_fK " << source_fK << endl << endl;
}

double SourceUniform::SurfaceBrightness(double *y){
	return (double)( (y[0]*y[0] + y[1]*y[1]) < source_r*source_r );
}

double SourceGaussian::SurfaceBrightness(double *y){
	return exp( -(y[0]*y[0] + y[1]*y[1])/source_gauss_r2 );
}

// surface brightness for models of the Broad Line Region
double SourceBLRDisk::SurfaceBrightness(double *y){
	return blr_surface_brightness_disk(y,this);
}

double SourceBLRSph1::SurfaceBrightness(double *y){
	return blr_surface_brightness_spherical_circular_motions(sqrt(y[0]*y[0] + y[1]*y[1]),this);
}
double SourceBLRSph2::SurfaceBrightness(double *y){
	return blr_surface_brightness_spherical_random_motions(sqrt(y[0]*y[0] + y[1]*y[1]),this);
}

void in_source(double *y_source,ListHndl sourcelist){
  return;
}

