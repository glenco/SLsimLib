/*
 * lens.c
 *
 *  Created on: Jan 12, 2012
 *      Author: mpetkova
 */

#include <slsimlib.h>
#include <sstream>

using namespace std;

Lens::Lens(string filename){
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

void Lens::readParamfile(string filename){
      const int MAXPARAM = 10;
	  string label[MAXPARAM], rlabel, rvalue;
	  void *addr[MAXPARAM];
	  int id[MAXPARAM];
	  stringstream ss;
	  int i ,n;
	  int myint;
	  double mydouble;
	  string mystring;

	  n = 0;

	  addr[n] = &outputfile;
	  id[n] = 2;
	  label[n++] = "outputfile";

	  addr[n] = &Nplanes;
	  id[n] = 2;
	  label[n++] = "Nplanes";

	  addr[n] = &zlens;
	  id[n] = 1;
	  label[n++] = "z_lens";

	  cout << "basic lens: reading from " << filename << endl;

	  ifstream file_in(filename.c_str());
	  if(!file_in){
	    cout << "Can't open file " << filename << endl;
	    exit(1);
	  }


	  // output file
	  while(!file_in.eof()){
		  file_in >> rlabel >> rvalue;

		  for(i = 0; i < n; i++){
			  if(rlabel == label[i]){

				  ss << rvalue;

				  switch(id[n]){
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
			  }
		  }
	  }


	  file_in.close();
}

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

	  label = "z_source";

	  cout << "uniform source: reading from " << filename << endl;

	  ifstream file_in(filename.c_str());
	  if(!file_in){
	    cout << "Can't open file " << filename << endl;
	    exit(1);
	  }

	  // output file
	  while(!file_in.eof()){
		  file_in >> rlabel >> rvalue;

		  if(rlabel == label){

			  ss << rvalue;

			  ss >> mydouble;
			  zsource = mydouble;

			  ss.clear();
			  ss.str(string());
		  }
	  }

	  file_in.close();

	  printSource();
}

void SourceGaussian::readParamfile(string filename){
	  string label[2], rlabel, rvalue;
	  stringstream ss;
	  void *addr[2];
	  int i;
	  double mydouble;

	  addr[0] = &zsource;
	  label[0] = "z_source";

	  addr[1] = &source_gauss_r2;
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

		  for(i = 0; i < 2; i++){

			  if(rlabel == label[i]){

				  ss << rvalue;
				  ss >> mydouble;

				  *((double *)addr[i]) = mydouble;

				  ss.clear();
				  ss.str(string());
			  }
		  }
	  }

	  file_in.close();

	  printSource();
}

void SourceBLR::readParamfile(string filename){
	  string label[9], rlabel, rvalue;
	  stringstream ss;
	  void *addr[9];
	  int i, n;
	  float myfloat;

	  n = 0;

	  addr[n] = &zsource;
	  label[n++] = "z_source";

	  addr[n] = &source_BHmass;
	  label[n++] = "BHmass";

	  addr[n] = &source_gamma;
	  label[n++] = "gamma";

	  addr[n] = &source_inclination;
	  label[n++] = "inclin";

	  addr[n] = &source_opening_angle;
	  label[n++] = "opening_ang";

	  addr[n] = &source_r_in;
	  label[n++] = "r_in";

	  addr[n] = &source_r_out;
	  label[n++] = "r_out";

	  addr[n] = &source_nuo;
	  label[n++] = "nuo";

	  addr[n] = &source_fK;
	  label[n++] = "source_sigma";

	  cout << "BLR source: reading from " << filename << endl;

	  ifstream file_in(filename.c_str());
	  if(!file_in){
	    cout << "Can't open file " << filename << endl;
	    exit(1);
	  }

	  // output file
	  while(!file_in.eof()){
		  file_in >> rlabel >> rvalue;

		  for(i = 0; i < 9; i++){

			  if(rlabel == label[i]){

				  ss << rvalue;
				  ss >> myfloat;

				  *((float *)addr[i]) = myfloat;

				  ss.clear();
				  ss.str(string());
			  }
		  }
	  }

	  file_in.close();

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
	cout << "incl " << source_inclination << endl;
	cout << "opening angl " << source_opening_angle << endl;
	cout << "r_in " << source_r_in << endl;
	cout << "r_out " << source_r_out << endl;
	cout << "nuo " << source_nuo << endl;
	cout << "source_sigma " << source_fK << endl << endl;
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
